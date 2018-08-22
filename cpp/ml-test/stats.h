/************************************************************
* Title: Summary statistics
*
* Date created: 14th July 2018
*
* Language: Cpp
*
* Overview:    
*
* Details: Header file
*
* Usage: 
*
************************************************************/

#define DEBUG_PRINT_DISTANCES

#include "iostream"
#include "Eigen/Dense"
#include "Eigen/SVD"

typedef Eigen::Matrix<std::complex<double>,
		      Eigen::Dynamic,
		      Eigen::Dynamic> MatrixXc;

// Distance using the operator norm
// Here, this is computed by finding
// the largest singular value, but
// I think it might be more complicated
// than that (as in, there are more
// variants of the norm to consider)
// 
double distance_op(MatrixXc A, MatrixXc B);


// Distance using the Frobenius norm
//
// This uses a library function from
// Eigen which computes the Frobenius
// (aka Hilbert-Schmidt) norm by
// default
//
// The Hilbert-Schmidt norm is defined
// in '2013 Sugiyama et al. - Precision guaranteed
// quantum tomography' as
//
//  1/sqrt(2) * tr[ (A - b)^2 ] ^ (1/2) 
//
double distance_hs(MatrixXc A, MatrixXc B);

// Distance using the trace norm
//
// The trace distance is defined
// in '2013 Sugiyama et al. - Precision
// guaranteed quantum tomography' as
//
//  1/2 * tr |A - B|
//
// where the 
double distance_trace(MatrixXc A, MatrixXc B);

// Distance using the l2 norm
//
// The l2 distance is the square root of the
// sum of the squares of the difference of
// the matrices
//
double distance_l2norm(MatrixXc A, MatrixXc B);

// Fidelity distance
//
// The fidelity distance d is
// defined as follows:
//
//  F(A,B) = tr[ sqrt(sqrt(A) B sqrt(A)) ]
//  d(A.B) = arccos[F(A,B)]
//
// This is computed using a trick that
// involves diagonalising A to ease the
// sqrt function. This originated as a
// python optimisation -- it might not
// be necessary here.
//
// distance_fid_2 implements the same
// operation in a different method
double distance_fid(const MatrixXc A, const MatrixXc B);
double distance_fid_2(const MatrixXc A, const MatrixXc B);

// Infidelity
//
// Compute the infidelity of two states using
//
//  F(A,B) = 1 - tr[ sqrt(A) B sqrt(A) ]
//
double infidelity(const MatrixXc A, const MatrixXc B);

// Infidelity2
//
// Compute the infidelity of two states using
//
//  F(A,B) = 1 - tr[ sqrt(sqrt(A) B sqrt(A)) ]^2
//
// Computed using diagonalisation tricks
double infidelity2(const MatrixXc A, const MatrixXc B);

// Mean calcluator
//
// Compute the mean of an array
//
double mean(double array[], int N);
