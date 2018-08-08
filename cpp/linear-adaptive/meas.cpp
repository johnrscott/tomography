/************************************************************
* Title: Measurement functions
*
* Date created: 7th August 2018
*
* Language: Cpp
*
* Overview: This file contains various functions for
*           manipulating measurement operators.
*
* Details:
*
* Usage: 
*
************************************************************/

#include "meas.h"

// Function: rotate one vector about another
//
// Rotate v about unit vector u
// by angle p
Vector3d axis_rot(Vector3d u, Vector3d v, double p){

  // Rotation of one vector about unit axis (u_x, u_y, u_z) by angle p
  MatrixXd R(3,3);
  
  double r11 = std::cos(p) + std::pow(u[0],2) * (1 - std::cos(p));
  double r22 = std::cos(p) + std::pow(u[1],2) * (1 - std::cos(p));
  double r33 = std::cos(p) + std::pow(u[2],2) * (1 - std::cos(p));
  double r12 = u[0] * u[1] * (1 - std::cos(p)) - u[2] * std::sin(p);
  double r21 = u[0] * u[1] * (1 - std::cos(p)) + u[2] * std::sin(p);
  double r13 = u[0] * u[2] * (1 - std::cos(p)) + u[1] * std::sin(p);
  double r31 = u[0] * u[2] * (1 - std::cos(p)) - u[1] * std::sin(p);
  double r23 = u[1] * u[2] * (1 - std::cos(p)) - u[0] * std::sin(p);
  double r32 = u[1] * u[2] * (1 - std::cos(p)) + u[0] * std::sin(p);
    
  R << r11,r12,r13,
       r21,r22,r23,
       r31,r32,r33;
  
  //#print('\n',Ru)
  MatrixXc I(3,3); I << 1, 0, 0,
		        0, 1, 0,
		        0, 0, 1;
  //RuT = R.adjoint();
  //#print('\n',RuT)
  //#print(np.matmul(Ru, RuT))
  //if((np.matmul(Ru, RuT) - I).max() > 1e-10) : print("Rotation matrix incorrect"); exit()
    
  // Perform rotation
  Vector3d w = R * v;
  return w;
}    

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

int straddle(Vector3d V, MatrixXc M1, MatrixXc M2, MatrixXc M3) {

  // Normalise V
  V.normalize();
    
  double theta = std::acos(1/std::sqrt(3));

  // Pick an arbitrary axis and use Gram-Scmidt to orthogonalize
  Vector3d y; y << 1,1,0;
  Vector3d z = y - V.dot(y) * V;
  z.normalize();
  Vector3d W1 = axis_rot(z, V, theta);

  // Rotate by 120 degrees
  double p = (2 * M_PI)/3;

  // Obtain the other measurement axes
  Vector3d W2 = axis_rot(V, W1, p);
  Vector3d W3 = axis_rot(V, W1, 2*p);

  // Check inner products between basis Bloch vectors
  // They should be orthonormal
  //print("\nW1.W2", np.dot(W1,W2))
  //print("W1.W3", np.dot(W1,W3))
  //print("W2.W3", np.dot(W2,W3))

  // Generate the measurement matrices
  MatrixXc I(2,2); I << 1, 0, 0, 1;
  MatrixXc X(2,2); X << 0, 1, 1, 0;
  MatrixXc Y(2,2); Y << 0, std::complex<double>(0,-1),
		     std::complex<double>(0,1), 0;
  MatrixXc Z(2,2); Z << 1, 0, 0, -1;

  M1 = W1[0]*X + W1[1]*Y + W1[2]*Z;
  M2 = W2[0]*X + W2[1]*Y + W2[2]*Z;
  M3 = W3[0]*X + W3[1]*Y + W3[2]*Z;

#ifdef DEBUG
#ifdef DEBUG_PRINT_MEASUREMENT_OPS
  std::cout << "===============================" << std::endl;
  std::cout << "The measurement operators are:" << std::endl;
  std::cout << std::endl;
  std::cout << "M1:" << std::endl << M1 << std::endl << std::endl;
  std::cout << "M2:" << std::endl << M2 << std::endl << std::endl;
  std::cout << "M3:" << std::endl << M3 << std::endl << std::endl;
  std::cout << std::endl;
#endif
#ifdef DEBUG_PRINT_INNER_PRODS
  std::cout << "===============================" << std::endl;
  std::cout << "The inner products between the " << std::endl;
  std::cout << "measurement operators are:" << std::endl;
  std::cout << std::endl;
  std::cout << "M1 . M2 = " << (M1 * M2).trace() << std::endl << std::endl;
  std::cout << "M1 . M3 = " << (M1 * M3).trace() << std::endl << std::endl;
  std::cout << "M2 . M3 = " << (M2 * M3).trace() << std::endl << std::endl;
  std::cout << std::endl;
#endif
#endif

  abort();
  
  return 0;
  
}
