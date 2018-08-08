/************************************************************
* Title: Functions for measurement operators
*
* Date created: 14th July 2018
*
* Language: Cpp
*
* Overview: Function for generating straddled basis
*
* Details: header
*
* Usage: Include in the main program
*
************************************************************/


#define DEBUG
#define DEBUG_PRINT_MEASUREMENT_OPS
#define DEBUG_PRINT_INNER_PRODS
#include <iostream>
#include "Eigen/Dense"
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>

typedef Eigen::Matrix<double,
		      Eigen::Dynamic,
		      Eigen::Dynamic> MatrixXd;

typedef Eigen::Matrix<double,3,1> Vector3d;


typedef Eigen::Matrix<std::complex<double>,
		      Eigen::Dynamic,
		      Eigen::Dynamic> MatrixXc;

// Function: straddle
//
// This function finds the measurement basis which
// 'straddles' a vector V on the Bloch sphere, meaning
// it is orthonormal and the vector sits equidistant
// from all the measurement operators.
// 
// First, rotate the vector V in an arbitrary direction by
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
// The three arrays are passed in by reference
//
int straddle(Vector3d V, MatrixXc M1, MatrixXc M2, MatrixXc M3);
