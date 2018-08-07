#############################################################
# Title: Measurement functions
#
# Date created: 26th June 2018
#
# Language: Python 3
#
# Overview: This file contains various functions for
#           manipulating measurement operators.
#
# Details:
#
# Usage: to be imported
#
#############################################################

import numpy as np
from numpy import sin, cos
import simulation

# Function: rotate one vector about another
#
# Rotate v about unit vector u
# by angle p
def axis_rot(u,v,p) :

    # Rotation of one vector about unit axis (u_x, u_y, u_z) by angle p
    r11 = cos(p) + u[0]**2 * (1 - cos(p))
    r22 = cos(p) + u[1]**2 * (1 - cos(p))
    r33 = cos(p) + u[2]**2 * (1 - cos(p))
    r12 = u[0] * u[1] * (1 - cos(p)) - u[2] * sin(p)
    r21 = u[0] * u[1] * (1 - cos(p)) + u[2] * sin(p)
    r13 = u[0] * u[2] * (1 - cos(p)) + u[1] * sin(p)
    r31 = u[0] * u[2] * (1 - cos(p)) - u[1] * sin(p)
    r23 = u[1] * u[2] * (1 - cos(p)) - u[0] * sin(p)
    r32 = u[1] * u[2] * (1 - cos(p)) + u[0] * sin(p)
    
    Ru = np.array([[r11,r12,r13],[r21,r22,r23],[r31,r32,r33]])
    #print('\n',Ru)
    I = np.identity(3)
    RuT = np.matrix.getH(np.asmatrix(Ru))
    #print('\n',RuT)
    #print(np.matmul(Ru, RuT))
    if((np.matmul(Ru, RuT) - I).max() > 1e-10) : print("Rotation matrix incorrect"); exit()
    
    # Perform rotation
    w = np.matmul(Ru, v)
    return w


# Function: straddle
#
# This function finds the measurement basis which
# 'straddles' a vector V on the Bloch sphere, meaning
# it is orthonormal and the vector sits equidistant
# from all the measurement operators.
# 
# First, rotate the vector V in an arbitrary direction by
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
def straddle(V) :

    # Normalise V
    V = V / np.linalg.norm(V)
    
    theta = np.arccos(1/np.sqrt(3))

    Y = np.array([1,1,0])
    Z = Y - np.dot(V,Y) * V
    axis = Z / np.linalg.norm(Z)
    W1 = axis_rot(axis, V, theta)

    # Rotate by 120 degrees
    p = (2 * np.pi)/3

    # Obtain the other measurement axes
    W2 = axis_rot(V, W1, p)
    W3 = axis_rot(V, W1, 2*p)

    # Check inner products between basis Bloch vectors
    # They should be orthonormal
    #print("\nW1.W2", np.dot(W1,W2))
    #print("W1.W3", np.dot(W1,W3))
    #print("W2.W3", np.dot(W2,W3))

    # Generate the measurement matrices
    I = np.matrix([[1,0],[0,1]])
    X = np.matrix([[0,1],[1,0]])
    Y = np.matrix([[0,-1j],[1j,0]])
    Z = np.matrix([[1,0],[0,-1]])
    M1 = W1[0]*X + W1[1]*Y + W1[2]*Z
    M2 = W2[0]*X + W2[1]*Y + W2[2]*Z
    M3 = W3[0]*X + W3[1]*Y + W3[2]*Z


    #print(M1,"\n\n")
    #print(M2,"\n\n")
    #print(M3,"\n\n")
    #exit()
    
    # Return measurements
    return M1, M2, M3
