/************************************************************
* Title: Summary statistics
*
* Date created: 19th June 2018
*
* Language:  Cpp
*
* Overview:    
*
* Details: source
*
* Usage:
*
************************************************************/

#include "stats.h"

// Distance using the trace norm
//
// The trace distance is defined
// in '2013 Sugiyama et al. - Precision
// guaranteed quantum tomography' as
//
//  1/2 * tr |A - B|
//
// where the 
double distance_trace(MatrixXc A, MatrixXc B) {
  MatrixXc C = (A - B).cwiseAbs();
  double distance  = std::abs(0.5 * C.trace());
  
#ifdef DEBUG
#ifdef DEBUG_PRINT_DISTANCES
  std::cout << "The operator distance is: " << distance << std::endl;
#endif
#endif
  
  return distance;
}

// Distance using the l2 norm
//
// The l2 distance is the square root of the
// sum of the squares of the difference of
// the matrices
//
double distance_l2norm(MatrixXc A, MatrixXc B) {
  //MatrixXc C = (A - B).norm();
  double distance  = (A - B).norm();
  
#ifdef DEBUG
#ifdef DEBUG_PRINT_DISTANCES
  std::cout << "The l2 distance is: " << distance << std::endl;
#endif
#endif
  
  return distance;
}


// Distance using the Hilbert-Schmidt norm
//
// The Hilbert-Schmidt distance is defined
// in '2013 Sugiyama et al. - Precision
// guaranteed quantum tomography' as
//
//  1/sqrt(2) * tr[ (A - B)^2 ] ^ (1/2) 
//
double distance_hs(MatrixXc A, MatrixXc B){
  MatrixXc C = (A - B) * (A - B);
  double distance  = std::abs(1/std::sqrt(2) * C.trace());

#ifdef DEBUG
#ifdef DEBUG_PRINT_DISTANCES
  std::cout << "The trace distance is: " << distance << std::endl;
#endif
#endif
  
  return distance;
}

// chopit
//
//
int chopit(MatrixXc & input) {
  int rows = input.rows();
  int cols = input.cols();
  for(int r=0; r<rows; r++) {
    for(int c=0; c<cols; c++) {
      if(std::abs(input(r,c).real()) < 1e-15) input(r,c) = std::complex<double>(0,input(r,c).imag());
      if(std::abs(input(r,c).imag()) < 1e-15) input(r,c) = std::complex<double>(input(r,c).real(),0);
    }
  }
    return 0;
}

// Infidelity
//
// Compute the infidelity of two states using
//
//  F(A,B) = 1 - tr[ sqrt(sqrt(A) B sqrt(A)) ]^2
//
double infidelity(const MatrixXc A, const MatrixXc B) {
  Eigen::SelfAdjointEigenSolver<MatrixXc> eigen1(A);
  //MatrixXc sqrtA = A.sqrt();
  MatrixXc sqrtA = eigen1.operatorSqrt();
  //std::cout << "A is: \n\n" << sqrtA << "\n\n";
  //std::cout << "B is: \n\n" << B << "\n\n";
  Eigen::SelfAdjointEigenSolver<MatrixXc> eigen2(sqrtA * B * sqrtA);
  //MatrixXc thingy = sqrtA * B * sqrtA;
  MatrixXc D = eigen2.operatorSqrt();
  //chopit(D);
  //std::cout << "Next is: \n\n" << D << "\n\n";
  //Eigen::SelfAdjointEigenSolver<MatrixXc> eig1(D);
  //std::cout << "Eig1: " << eig1.eigenvalues()[0] << std::endl; 
  //std::cout << "Eig2: " << eig1.eigenvalues()[1] << std::endl;
  //abort();
  //MatrixXc C = sqrtA * B * sqrtA;
  double infidelity = 1 - std::pow(std::real(D.trace()),2);
  return infidelity;
}

// Infidelity2
//
// Compute the infidelity of two states using
//
//  F(A,B) = 1 - tr[ sqrt(sqrt(A) B sqrt(A)) ]^2
//
// Computed using diagonalisation tricks. Note that the
// input must be a positive matrix (or very close) for the
// results to be valid
double infidelity2(const MatrixXc A, const MatrixXc B) {
  // Diagonalise A. Put eigenvectors as columns of P.
  Eigen::SelfAdjointEigenSolver<MatrixXc> eigenA(A);
  if(eigenA.info() != Eigen::Success) abort();
  MatrixXc P = eigenA.eigenvectors();
  // Change the basis of B: let T = Inv(P) * B * P
  MatrixXc T = P.inverse() * B * P;
  // Get A^(1/2) in the new basis: S, with sqare-rooted eigenvalues of A on the diagonal
  // Assume that the eigenvalues are positive.
  MatrixXc S(2,2); S << std::sqrt(std::abs(eigenA.eigenvalues()[0])), 0,
  		        0, std::sqrt(std::abs(eigenA.eigenvalues()[1]));
  // Get the product R = S * T * S
  MatrixXc R = S * T * S;
  // Diagonalise R
  Eigen::SelfAdjointEigenSolver<MatrixXc> eigenR(R);
  if(eigenR.info() != Eigen::Success) abort();
  // The fidelity is square of the sum of the square roots of the eigenvalues
  double fid = std::pow(std::sqrt(eigenR.eigenvalues()[0]) + std::sqrt(eigenR.eigenvalues()[1]),2); 
  double infid = 1 - fid;
  return infid;
}

// Mean calcluator
//
// Compute the mean of an array
//
double mean(double array[], int N) {
  
  double tmp(0);
  for(int n=0; n<N; n++) tmp += array[n];
  double mean = tmp/N;
  return mean;

}
