/************************************************************
* Title: Testing the linear adaptive estimator
*
* Date created: 14th July 2018
*
* Language:    Cpp
*
* Overview:    The program uses a simple adaptive strategy
*              to attempt to improve the linear estimate.
*              The first part of the data is used to choose
*              a measurement basis, which is then used for
*              all subsequent measurements.
*
* Details:     The steps are similar to the other tests:
*
*              1) Pick a random density matrix for 2 qubits
*
*              2) Fix X, Y Z measurement operators
*
*              3) Compute the density matrix using the
*                 linear estimator using a certain (small)
*                 number of X, Y, Z measurements. This will
*                 serve as a prelimiary estimate for 
*                 refining the basis.
*
*              4) Pick a new basis which is centred on the
*                 density matrix estimate
* 
*              5) Perform all the other measurements in this
*                 new basis
*
*              6) Compute the density matrix using the
*                 linear estimator
*
*              6) Compute the distance between the
*                 estimate and the true density matrix
*
*              7) Perform averages of each distance at each
*                 value of x. 
*
*              8) Repeat for different values of x and repeat 
*                 each value of x a large number of times. 
*                 Store all the distances.
*
* Compilation: make 
*
*************************************************************/

#include "linear-adaptive.h"

int main() {

  // Start the clock!
  std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

  
  // Step 0: Seed the random number generators
  //
  // The random engine is initialised outside the main program
  // This is so that the engine is only initialised once,
  // as opposed to every time the simulate function is
  // called. That way each set of measurements is independent
  // of the others. In fact, it is critical that it is passed
  // by reference -- otherwise the object gets copied and
  // you get the same problem. If the generator is passed
  // by reference, the same object is used by all the
  // random number generators.
  //
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  srand(seed);
  std::mt19937 gen;
  gen.seed(seed);

  // ======= Test parameter ===============================
  int M = 2000;  // Number of purity parameters x to try
  double x_start = 0; // Specify purity parameter x range
  double x_end = 1;
  int N = 500;  // Number of random density matrices per x value
  int S = 500;  // Number of samples of each measurement to
                // simulate for each density matrix
  int S_1 = 100; // Number of samples of X, Y and Z to use to
                 // perform the prelimiary estimate to adapt
                 // the basis.
  // ======================================================

  double non_physical[M];

  // Get an output file ready
  std::ofstream file;
  file.open("linear_test_1_cpp.dat");
  file << "Distances between estimated and "
       << "original density matrices using various distances."
       << "Go to the end of the file for the running time."
       << std::endl << std::endl;

  // Write the simulation parameters
  file << "Number of purity values tried = "
       << M << std::endl
       << "Number of density matrices per purity parameter = "
       << N << std::endl
       << "Total number of measurements for each of X, Y and Z = "
       << S << std::endl
       << std::endl;
  
  file << "PURITY, \tOPERATOR, \tTRACE, \t\tFIDELITY, \tNON PHYSICAL";
  // Set precision
  file << std::fixed << std::setprecision(5) << std::endl;

  // Preliminaries: define measurement operators
  MatrixXc I(2,2); I << 1, 0, 0, 1;
  MatrixXc X(2,2); X << 0, 1, 1, 0;
  MatrixXc Y(2,2); Y << 0, std::complex<double>(0,-1),
		     std::complex<double>(0,1), 0;
  MatrixXc Z(2,2); Z << 1, 0, 0, -1;

  // Compute the projectors
  MatrixXc proj_X[2];
  MatrixXc proj_Y[2];
  MatrixXc proj_Z[2];
  double outcomes_X[2];
  double outcomes_Y[2];
  double outcomes_Z[2];
  
  make_projector(X, proj_X, outcomes_X);
  make_projector(Y, proj_Y, outcomes_Y);
  make_projector(Z, proj_Z, outcomes_Z);
  
  // Define x -- put a loop here ------------------- LOOP for x between 0 and 1
  //
  // This loop runs through different values of the purity parameter x,
  // and tests the ability of the linear estimator in each case
  //

  // Variables to store the estimation error distances
  double dist_op[N];
  double dist_trace[N];
  double dist_fid[N];

  //int dp = 5; // Decimal places for printing
  double x = 0; // Purity parameter
  
  for(int k=0; k<M; k++) {
#ifdef DEBUG
    std::cout << "======================= "
	      << "START OF AN X LOOP"
	      << "======================="
	      << std::endl;
#endif
    // Temporary counter for non-physical estimates
    double non_physical_count = 0;
    
    // Loop N times for each value of x ------ inner loop -- N trials for each x
    //
    // This loop generates N random density matrices for each fixed value of x
    // which used to simulate measurement data and run the estimator
    //
    for(int n=0; n<N; n++) {
#ifdef DEBUG
      std::cout << std::endl
		<< "++++++++ "
		<< "FIXED DENSITY MATRIX "
		<< "++++++++"
		<< std::endl << std::endl;
#endif
      // Step 1: Prepare the density matrix
      //
      // The purity parameter x is picked between 0
      // and 1.
      //
      // Note: any time a numerical check is performed
      // and printed out, I've rounded the result to
      // make it more readable. Set the decimal places
      // to keep using the dp variable.
      //
      x = x_start + k * (x_end - x_start)/M;
      MatrixXc dens = random_density(x, gen);
      //MatrixXc dens(2,2); dens << 1,0,0,0;
      
      // Step 2: Generate measurement data
      //
      // Generate data for X, Y and Z measurements. 
      //
      double X_data[S_1]; // To contain the (real) measurement values
      double Y_data[S_1];
      double Z_data[S_1];
      simulate(dens, proj_X, outcomes_X, S_1, X_data, gen);
      simulate(dens, proj_Y, outcomes_Y, S_1, Y_data, gen);
      simulate(dens, proj_Z, outcomes_Z, S_1, Z_data, gen);
      
      // Step 3: Estimate density matrix
      //
      // Compute linear estimator
      //
      // Then tr(pI) is computed by requiring that
      // the density matrix be normalised
      //
      // This estimate is used to refine the
      // measurement basis.
      //
      MatrixXc dens_est_adapt = linear_estimate_XYZ(X_data, Y_data, Z_data, S_1);

      // Step 4: Get new basis
      //
      // First, find the Bloch vector for the density matrix.
      // This is obtained using the Hilbert-Schmidt inner
      // product (A,B) = Tr(A^ B):
      //
      // V = 1/2 * [Tr(dens_est * X),
      //            Tr(dens_est * Z),
      //            Tr(dens_est * Z)]
      //
      // Then, rotate the vector V in an arbitrary direction by
      // 0.95531 rad (54 degrees).
      //
      // The new Bloch vector W will be one of the directions for the
      // new measurements. Then rotate this new vector around the first
      // Bloch vector by 120 degrees and then 240 degrees. These vectors
      // will be the other two measurement Bloch vectors.
      //
      // Use the Bloch vectors to obtain Hermitian matrices using the
      // expansion H = (1/2)(I + V.X) where V is the Bloch vector and
      // X is the vector of Pauli matrices.
      //
      Vector3d V; V << (dens_est_adapt * X).trace().real(),
		       (dens_est_adapt * Y).trace().real(),
		       (dens_est_adapt * Z).trace().real();

      // Array of three 2x2 matrices to store operators
      MatrixXc M1(2,2), M2(2,2), M3(2,2);
      straddle(V,M1,M2,M3);

      // Compute the projectors and outcomes for M1, M2 and M3
      MatrixXc proj_M1[2];
      MatrixXc proj_M2[2];
      MatrixXc proj_M3[2];
      double outcomes_M1[2];
      double outcomes_M2[2];
      double outcomes_M3[2];
      make_projector(M1, proj_M1, outcomes_M1);
      make_projector(M2, proj_M2, outcomes_M2);
      make_projector(M3, proj_M3, outcomes_M3);

      // Step 5: Generate measurement data again
      //
      // Generate data for the M1, M2 and M3
      // measurement operators
      //
      double X_data[S-S_1]; // To contain the (real) measurement values
      double Y_data[S-S_1];
      double Z_data[S-S_1];
      simulate(dens, proj_X, outcomes_X, S-S_1, X_data, gen);
      simulate(dens, proj_Y, outcomes_Y, S-S_1, Y_data, gen);
      simulate(dens, proj_Z, outcomes_Z, S-S_1, Z_data, gen);

      abort()
      
      // Re
      //
      MatrixXc dens_est = linear_estimate_XYZ(X_data, Y_data, Z_data, S);
      
      // Step 4: Compute and the distances
      //
      // Compute distances between the estimated
      // and true density matrix using the
      // different distance fuctions.
      //
      dist_op[n] = distance_op(dens_est, dens);
      dist_trace[n] = distance_trace(dens_est, dens);
      dist_fid[n] = distance_fid_2(dens_est, dens);

      // Count the number of non-physical matrices
      //
      Eigen::SelfAdjointEigenSolver<MatrixXc> eigenD(dens_est);
      if(eigenD.info() != Eigen::Success) abort();

#ifdef DEBUG
#ifdef DEBUG_PRINT_ESTIMATES_EIGENVALUES
      std::cout << "The estimated eigenvalues are "
		<< eigenD.eigenvalues()
		<< std::endl;
#endif
#endif
      if ((eigenD.eigenvalues()[0] < 0) || (eigenD.eigenvalues()[1] < 0)) {
	non_physical_count = non_physical_count + 1;
      }
      
    } // end of inner loop (fixed purity, random density matrices)
    
    // Step 5: Average the distances 
    //
    // Average the distances for each value of x. There are N density
    // matrices for each value of X, and consequently N distances to
    // compute. Therefore the mean is computed over N points.
    //
    double mean_op = mean(dist_op, N);
    double mean_trace = mean(dist_trace, N);
    double mean_fid = mean(dist_fid, N);
    non_physical[k] = non_physical_count/N;

#ifdef SHOW_PROGRESS
    // Show progress
    double p = static_cast<double>(k+1)/M;
    show_progress(start,p);
#endif
    
    
#ifdef DEBUG
#ifdef DEBUG_PRINT_DISTANCE_AVERAGES
    std::cout << std::endl
	      << "The average operator distance is: " << mean_op << std::endl    
	      << "The average trace distance is: " << mean_trace << std::endl
	      << "The average fidelity distance is: " << mean_fid << std::endl;
#endif
#endif
    
    // Step 6: Write the results to a file
    //
    file << x << ",\t"
	 << mean_op << ",\t"
	 << mean_trace << ",\t"
	 << mean_fid << ",\t"
	 << non_physical[k]
	 << std::endl;

#ifdef DEBUG
    std::cout << std::endl
	      << "++++++++ "
	      << "END OF INNER LOOP "
	      << "++++++++"
	      << std::endl << std::endl;
#endif
    
  } // end of outer loop (looping through different purities)

  // Get stop time
  auto end = std::chrono::steady_clock::now();
  auto time = end - start;

  auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(time).count();
  
  // Store the running time
  file << std::endl
       << "Total running time = "
       << dur << "ms" << std::endl;
    
  file.close();
}
