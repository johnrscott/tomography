/************************************************************
* Title: Density matrix estimation
*
* Date created: 14th July 2018
*
* Language: Cpp
*
* Overview:    
*
* Details:
*
* Usage: to be imported
*
************************************************************/

#include "estimation.h"

// Function: linear_estimate_XYZ(X_data, Y_data, Z_data)
//
// The function takes in X, Y and Z
// measurement data and produces an
// estimate of the density matrix
// by computing the means of the X,
// Y and z data. The estimated density
// matrix is given by
//
// dens_est = 1/2 (mean_X // X +
//                 mean_Y // Y +
//                 mean_Z // Z +
//                 I) 
//
// This function returns a density
// matrix that satisfies the trace = 1
// constraint. This is achieved by
// leaving the multiple of I as a
// free parameter to be set after
// the other parameters have been
// estimated.
//
// The matrix is also Hermitian by
// construction. The only property
// not necessarily satisfied by the
// density matrix is positivity.
//
// S contains the number of samples
// in each of X, Y and Z
//
//
MatrixXc linear_estimate_XYZ(double X_data[],
			     double Y_data[],
			     double Z_data[],
			     int S) {
  MatrixXc I(2,2); I << 1, 0, 0, 1;
  MatrixXc X(2,2); X << 0, 1, 1, 0;
  MatrixXc Y(2,2); Y << 0, std::complex<double>(0,-1),
		     std::complex<double>(0,1), 0;
  MatrixXc Z(2,2); Z << 1, 0, 0, -1;

  // Compute means
  double mean_X = mean(X_data, S);
  double mean_Y = mean(Y_data, S);
  double mean_Z = mean(Z_data, S);
  
#ifdef DEBUG
#ifdef DEBUG_PRINT_MEANS
  std::cout << std::endl
	    << "X mean: " << mean_X << std::endl
	    << "Y mean: " << mean_Y << std::endl
	    << "Z mean: " << mean_Z << std::endl
	    << std::endl;
  
#endif
#endif
  
  // Reconstruct density matrix
  MatrixXc dens_est = (mean_X * X + mean_Y * Y + mean_Z * Z + I)/2;

#ifdef DEBUG
#ifdef DEBUG_PRINT_ESTIMATE
  std::cout << "The estimated density matrix is"
	    << std::endl << std::endl
	    << dens_est
	    << std::endl << std::endl;
#endif
#endif
  
  return dens_est;
}

// The objective function d is the distance between two matrices
// parametrised by X and Y.
//
// The last parameter void * is expected to be MatrixXc *
double d(const std::vector<double> & x, std::vector<double> & grad, void * f_data) {
  // x specifies dens = T T^
  MatrixXc T(2,2); T << x[0], 0, std::complex<double>(x[1],x[2]), x[3];
  MatrixXc dens_1 = T * T.adjoint();
  MatrixXc dens_2 = * static_cast<MatrixXc * >(f_data);
  //std::cout << "The linear_estimate is:\n\n" << dens_2 << std::endl;
  double distance = distance_l2norm(dens_1, dens_2);
  return distance;
}

// Constrain data for the trace 1 condition
double cons(const std::vector<double> & x, std::vector<double> & grad, void * f_data) {
  double value = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] - 1;
  return value;
}

// Function: enm_estimate_XYZ(X_data, Y_data, Z_data)
//
// This function estimates the density matrix using
// the extended norm minimisation method.
//
// This estimation method first uses the linear
// estimator to get an initial estimate of the
// density matrix, and then using a minimisation
// procedure to find the closest physical density
// matrix to this estimate using the trace norm.
//
MatrixXc enm_estimate_XYZ(double X_data[],
			  double Y_data[],
			  double Z_data[],
			  int S) {  
  // Estimate the density matrix with the linear estimator
  MatrixXc dens_lin = linear_estimate_XYZ(X_data, Y_data, Z_data, S);
  //std::cout << "The linear_estimate is:\n\n" << dens_lin << std::endl;

  // Something that would speed up the algorithm would be
  // to check whether the density matrix is physical here,
  // then skip the optimisation procedure if it is.
  //Eigen::SelfAdjointEigenSolver<MatrixXc> eigenD(dens_lin);
  //if(eigenD.info() != Eigen::Success) abort();
  //if((eigenD.eigenvalues()[0] < 1) && (eigenD.eigenvalues()[0] > 0)) {
  //  return dens_lin; // Return dens_est -- no optimisation to do
  //}

  // Declare variables outside try
  MatrixXc T(2,2); 
  MatrixXc dens_enm;
  
  try {
    // Create an nlop object
    // For some reason SLSQP appears not to work
    nlopt::opt opt(nlopt::/*LD_SLSQP*/LN_NELDERMEAD, 4); // 4 optimisation parameters
    // Set objective function
    opt.set_min_objective(d, &dens_lin);
    // Add trace 1 constraint
    opt.add_equality_constraint(cons, &dens_lin);
    opt.set_ftol_rel(1e-8);
    double ftol_rel = opt.get_ftol_rel();
    //std::cout << ftol_rel;
    //abort();
    std::vector<double> x{0,0,0,0};
  
    // Variable to contain the minimum of the function
    double minf;
    
    // Do the optimisation
    nlopt::result result = opt.optimize(x,minf);
    // Reconstruct the density matrix
    T << x[0], 0, std::complex<double>(x[1],x[2]), x[3];
    dens_enm = T * T.adjoint();
    //#ifdef DEBUG
    //#ifdef DEBUG_PRINT_ENM_OUTPUT
    std::cout << "Linear estimator: " << std::endl
	      << dens_lin << std::endl;
    Eigen::SelfAdjointEigenSolver<MatrixXc> eig1(dens_lin);
    std::cout << "Eig1: " << eig1.eigenvalues()[0] << std::endl; 
    std::cout << "Eig2: " << eig1.eigenvalues()[1] << std::endl;
    
    std::cout << "ENM estimator: " << std::endl
	      << dens_enm << std::endl;
    Eigen::SelfAdjointEigenSolver<MatrixXc> eig2(dens_enm);
    std::cout << "Eig1: " << eig2.eigenvalues()[0] << std::endl; 
    std::cout << "Eig2: " << eig2.eigenvalues()[1] << std::endl;
    abort();
    //#endif
    //#endif
  }
  catch (std::exception & e) {
    std::cout << "nlopt failed: " << e.what() << std::endl;
    abort();
  }
    
  return dens_enm;
    
}

// The likelihood function is used by the maximum likelihood
// method to find the most likely density matrix for a
// given set of expectation values
//
// The last parameter void * contains the expectatation values
//
double L(const std::vector<double> & x, std::vector<double> & grad, void * f_data) {
  // x specifies dens = T T^
  MatrixXc T(2,2); T << x[0], 0, std::complex<double>(x[1],x[2]), x[3];
  MatrixXc dens_1 = T * T.adjoint();
  MatrixXc dens_2 = * static_cast<MatrixXc * >(f_data);
  //std::cout << "The linear_estimate is:\n\n" << dens_2 << std::endl;
  double distance = distance_l2norm(dens_1, dens_2);
  return distance;
}

// Function: ml_estimate_XYZ(X_data, Y_data, Z_data)
//
// This function estimates the density matrix using
// the maximum likelihood method
//
// The maximum likelihood method using a likelihood
// function to determine
MatrixXc ml_estimate_XYZ(double X_data[],
			  double Y_data[],
			  double Z_data[],
			  int S) {  
  // Estimate the density matrix with the linear estimator
  MatrixXc dens_lin = linear_estimate_XYZ(X_data, Y_data, Z_data, S);
  //std::cout << "The linear_estimate is:\n\n" << dens_lin << std::endl;

  // Something that would speed up the algorithm would be
  // to check whether the density matrix is physical here,
  // then skip the optimisation procedure if it is.
  //Eigen::SelfAdjointEigenSolver<MatrixXc> eigenD(dens_lin);
  //if(eigenD.info() != Eigen::Success) abort();
  //if((eigenD.eigenvalues()[0] < 1) && (eigenD.eigenvalues()[0] > 0)) {
  //  return dens_lin; // Return dens_est -- no optimisation to do
  //}

  // Create an nlop object
  // For some reason SLSQP appears not to work
  nlopt::opt opt(nlopt::/*LD_SLSQP*/LN_NELDERMEAD, 4); // 4 optimisation parameters
  // Set objective function
  opt.set_min_objective(d, &dens_lin);
  opt.set_ftol_rel(1e-8);
  double ftol_rel = opt.get_ftol_rel();
  //std::cout << ftol_rel;
  //abort();
  std::vector<double> x{0,0,0,0};

  // Variable to contain the minimum of the function
  double minf;

  // Declare variables outside try
  MatrixXc T(2,2); 
  MatrixXc dens_enm;
  
  // Do the optimisation
  try {
    nlopt::result result = opt.optimize(x,minf);
    // Reconstruct the density matrix
    T << x[0], 0, std::complex<double>(x[1],x[2]), x[3];
    dens_enm = T * T.adjoint();
    #ifdef DEBUG
    #ifdef DEBUG_PRINT_ENM_OUTPUT
    std::cout << "Linear estimator: " << std::endl
	      << dens_lin << std::endl;
    Eigen::SelfAdjointEigenSolver<MatrixXc> eig1(dens_lin);
    std::cout << "Eig1: " << eig1.eigenvalues()[0] << std::endl; 
    std::cout << "Eig2: " << eig1.eigenvalues()[1] << std::endl;
    
    std::cout << "ENM estimator: " << std::endl
	      << dens_enm << std::endl;
    Eigen::SelfAdjointEigenSolver<MatrixXc> eig2(dens_enm);
    std::cout << "Eig1: " << eig2.eigenvalues()[0] << std::endl; 
    std::cout << "Eig2: " << eig2.eigenvalues()[1] << std::endl;
    abort();
    #endif
    #endif
  }
  catch (std::exception & e) {
    std::cout << "nlopt failed: " << e.what() << std::endl;
    abort();
  }
    
  return dens_enm;
    
}
