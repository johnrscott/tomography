#!/usr/local/bin/python3
#############################################################
# Title: Linear adaptive test
#
# Date created: 19th June 2018
#
# Language:    Python 3
#
# Overview:    The program uses a simple adaptive strategy
#              to attempt to improve the linear estimate.
#              The first part of the data is used to choose
#              a measurement basis, which is then used for
#              all subsequent measurements.
#
# Details:     The script performs the following steps:
#
# Usage: python3 linear-adaptive.py
#
#############################################################

# Include
import importlib
import numpy as np
from numpy import sin, cos
from scipy.stats import unitary_group as ug
import scipy as sc
import simulation
importlib.reload(simulation)
import estimation
importlib.reload(estimation)
import stats
importlib.reload(stats)
import estimation
importlib.reload(estimation)
import cProfile
import pstats
from progress import *
from mpmath import chop
import meas

pr = cProfile.Profile()
pr.enable()

# ======= Test parameter ===============================
M = 2000  # Number of purity parameters x to try
x_start = 0 # Specify purity parameter x range
x_end = 1
N = 500  # Number of random density matrices per x value
S = 500  # Number of samples of each measurement to
         # simulate for each density matrix
S_1 = 100 # Number of samples of X, Y and Z to use to
         # perform the prelimiary estimate to adapt
         # the basis.
# ======================================================

# Seed random number generator
np.random.seed()

av_distances = np.zeros([M,3])
non_physical = np.zeros(M) # Proportion of non-physical estimates

# Preliminaries: compute the projectors

I = np.matrix([[1,0],[0,1]])
X = np.matrix([[0,1],[1,0]])
Y = np.matrix([[0,-1j],[1j,0]])
Z = np.matrix([[1,0],[0,-1]])

# This is a terrible function
proj_X, proj_Y, proj_Z, values_X, values_Y, values_Z = simulation.projectors(X,Y,Z)

# Open a file for writing
file = open("linear_adaptive_1_python.dat", "w")
file.write("Distances between estimated and original density matrices using various distances:\n\n")
file.write("Number of purity values tried = "+str(M)+"\n")
file.write("Number of density matrices per purity parameter = "+str(N)+"\n")
file.write("Total number of measurements = "+str(S)+"\n\n");
file.write("PURITY, \tOPERATOR, \tTRACE, \t\tFIDELITY, \tNON PHYSICAL\n")

# Define x -- put a loop here ----------------------------------------- LOOP for x between 0 and 1
#
# This loop runs through different values of the purity parameter x,
# and tests the ability of the linear estimator in each case
#
dist_op = np.zeros([N,1])
dist_trace = np.zeros([N,1])
dist_fid = np.zeros([N,1])

dp = 5

#for k,n in itertools.product(range(M),range(N)):
for k in range(M):
    non_physical_count = 0 # Temporary counter for non-physical estimates
    
    # Loop N times for each value of x ------------------ inner loop -- N trials for each x
    #
    # This loop generates N random density matrices for each fixed value of x
    # which used to simulate measurement data and run the estimator
    #
    for n in range(N):
        # Step 1: Prepare the density matrix
        #
        # The purity parameter x is picked between 0
        # and 1.
        #
        # Note: any time a numerical check is performed
        # and printed out, I've rounded the result to
        # make it more readable. Set the decimal places
        # to keep using the dp variable.
        #
        x = x_start + k * (x_end - x_start)/M
        #pr.enable()
        dens = simulation.random_density(x)
        #values_dens,vectors_dens = np.linalg.eig(dens)
        #pr.disable()

        # Step 2: Generate measurement data
        #
        # Generate data for X, Y and Z measurements.
        # Only generate S_1 samples for each measurement.
        # These will be used to estimate a preliminary
        # density matrix, that will be used to adapt the
        # basis.
        #
        X_data = simulation.simulate(dens,proj_X,values_X,S_1)
        Y_data = simulation.simulate(dens,proj_Y,values_Y,S_1)
        Z_data = simulation.simulate(dens,proj_Z,values_Z,S_1)

        # Step 3: Estimate the preliminary density matrix
        #
        # Compute linear estimator
        #
        # Then tr(pI) is computed by requiring that
        # the density matrix be normalised
        #
        dens_est_adapt = estimation.linear_estimate_XYZ(X_data, Y_data, Z_data)

        # Step 4: Adapt the measurement basis
        #
        # First, find the Bloch vector for the density matrix.
        # This is obtained using the Hilbert-Schmidt inner
        # product (A,B) = Tr(A^ B):
        #
        # V = 1/2 * [Tr(dens_est * X),
        #            Tr(dens_est * Z),
        #            Tr(dens_est * Z)]
        #
        # Then, rotate the vector V in an arbitrary direction by
        # 0.95531 rad (54 degrees).
        #
        # The new Bloch vector W will be one of the directions for the
        # new measurements. Then rotate this new vector around the first
        # Bloch vector by 120 degrees and then 240 degrees. These vectors
        # will be the other two measurement Bloch vectors.
        #
        # Use the Bloch vectors to obtain Hermitian matrices using the
        # expansion H = (1/2)(I + V.X) where V is the Bloch vector and
        # X is the vector of Pauli matrices.
        #
        V = np.array([np.trace(np.matmul(dens_est_adapt,X)),
                      np.trace(np.matmul(dens_est_adapt,Y)),
                      np.trace(np.matmul(dens_est_adapt,Z))]).transpose()
        #       
        #       M1, M2, M3 = meas.straddle(V)
        #
        #       # Get projectors
        #       proj_M1, proj_M2, proj_M3, values_M1, values_M2, values_M3 = simulation.projectors(M1,M2,M3)
        
        M1, M2, M3 = meas.straddle(V)
        
        # Get projectors
        proj_M1, proj_M2, proj_M3, values_M1, values_M2, values_M3 = simulation.projectors(M1,M2,M3)

        # Step 5: Simulate new measurements in the new basis
        #
        # This step simulates data in the new measurement
        # basis. These new measurements are then used to
        # obtain the final estimate of the density matrix.
        #
        # There are only S-S_1 measurements left at this
        # point, because S_1 of them are supposed to have been
        # used up in estimating the basis
        M1_data = simulation.simulate(dens,proj_M1,values_M1,S-S_1)
        M2_data = simulation.simulate(dens,proj_M2,values_M2,S-S_1)
        M3_data = simulation.simulate(dens,proj_M3,values_M3,S-S_1)

        # Step 6: Estimate density matrix
        #
        # Compute linear estimator using
        # the M1, M2, M3 measurements
        #
        dens_est = estimation.linear_estimate(M1_data, M2_data, M3_data,
                                              M1, M2, M3)

        #print("\nThe original density matrix was:\n")
        #print(dens)
        #print("\nThe estimate is:\n")
        #print(dens_est)
        #exit()
        # Step 4: Compute and the distances
        #
        # Compute distances between the estimated
        # and true density matrix using the
        # different distance fuctions.
        #
        dist_op[n] = stats.distance_op(dens, dens_est)
        dist_trace[n] = stats.distance_trace(dens, dens_est)
        dist_fid[n] = stats.distance_fid(dens, dens_est)

        # Count the number of non-physical matrices
        #
        eigenvalues = np.linalg.eigvals(dens_est)
        if eigenvalues[0] < 0 or eigenvalues[1] < 0:
            non_physical_count = non_physical_count + 1
        
    # Step 5: Average the distances 
    #
    # Average the distances for each value of x
    #
    av_distances[k,:] = [np.mean(dist_op), np.mean(dist_trace), np.mean(dist_fid)]
    non_physical[k] = non_physical_count/N
    p = (k+1)/M
    show_progress(pr,p)
    file.write("{0:.5f},\t{1:.5f},\t{2:.5f}, \t{3:.5f}, \t{4:.5f}\n".format(x, np.mean(dist_op),
                                                                            np.mean(dist_trace),
                                                                            np.mean(dist_fid),
                                                                            non_physical[k]))


pr.disable()
ps = pstats.Stats(pr)
total_time = ps.total_tt

file.write("\nTotal running time = "+str(np.around(total_time,3))+"s\n")    
file.close
