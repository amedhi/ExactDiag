/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-16 22:12:57
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-27 23:49:33
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
//#define EIGEN_USE_MKL_ALL
#include <Eigen/Core>

#ifndef MATRIX_H
#define MATRIX_H

using Vector3i = Eigen::Vector3i;
using Vector3d = Eigen::Vector3d;
using RealVector = Eigen::VectorXd;
using intMatrix = Eigen::MatrixXi;
using RealMatrix = Eigen::MatrixXd;
using ComplexVector = Eigen::VectorXcd;
using ComplexMatrix = Eigen::MatrixXcd;
using cmplVector = Eigen::VectorXcd;
using realArray1D = Eigen::ArrayXd; 
using cmplArray1D = Eigen::ArrayXcd; 
using cmplArray2D = Eigen::Array<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic>;

#endif
