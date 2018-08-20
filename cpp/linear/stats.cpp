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

// Infidelity
//
// Compute the infidelity of two states using
//
//  F(A,B) = 1 - tr[ sqrt(A) B sqrt(A) ]
//
double infidelity(const MatrixXc A, const MatrixXc B) {
  Eigen::SelfAdjointEigenSolver<MatrixXc> eigen1(A);
  MatrixXc sqrtA = eigen1.operatorSqrt();
  Eigen::SelfAdjointEigenSolver<MatrixXc> eigen2(sqrtA * B * sqrtA);
  MatrixXc D = eigen2.operatorSqrt();
  //MatrixXc C = sqrtA * B * sqrtA;
  double infidelity = 1 - std::pow(std::real(D.trace()),2);
  return infidelity;
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
