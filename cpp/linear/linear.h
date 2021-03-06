/************************************************************
* Title: ENM Test header file
*
* Date created: 13th August 2018
*
* Language: Cpp
*
* Overview: 
*
* Details: Header
*
* Usage: Include the header file
*
************************************************************/

#define SHOW_PROGRESS

#include <iostream>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <random> 
#include <fstream>
#include <iomanip>
#include <chrono>
#include "Eigen/Dense"
//#include "simulation.h"
#include "estimation.h"
#include "stats.h"
//#include "proj.h"
#include "progress.h"
#include "qpp.h"
#include <tuple>

typedef Eigen::Matrix<std::complex<double>,
		      Eigen::Dynamic,
		      Eigen::Dynamic> MatrixXc;
