/*

   calin/calin_global_definitions.hpp -- Stephen Fegan -- 2015-04-16

   Definitions of types and macros that apply to the full package

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#pragma once

#include <vector>
#include <Eigen/Core>

namespace calin {

#ifndef SWIG
#define CALIN_TYPEALIAS(A,B...) using A = B
#else
#define CALIN_TYPEALIAS(A,B...) typedef B A
#endif

//#define CALIN_USE_EIGEN_REF

CALIN_TYPEALIAS(VecRef, Eigen::VectorXd&);
CALIN_TYPEALIAS(MatRef, Eigen::MatrixXd&);
CALIN_TYPEALIAS(ConstVecRef, const Eigen::VectorXd&);
CALIN_TYPEALIAS(ConstMatRef, const Eigen::MatrixXd&);

CALIN_TYPEALIAS(IntVecRef, Eigen::VectorXi&);
CALIN_TYPEALIAS(IntMatRef, Eigen::MatrixXi&);
CALIN_TYPEALIAS(ConstIntVecRef, const Eigen::VectorXi&);
CALIN_TYPEALIAS(ConstIntMatRef, const Eigen::MatrixXi&);

inline std::vector<double> eigen_to_stdvec(const Eigen::VectorXd& x)
{
  return std::vector<double>(x.data(), x.data()+x.size());
}

inline Eigen::VectorXd std_to_eigenvec(const std::vector<double> &x)
{
  return Eigen::Map<const Eigen::VectorXd>(&x.front(), x.size());
}

inline std::vector<int> eigen_to_stdvec(const Eigen::VectorXi& x)
{
  return std::vector<int>(x.data(), x.data()+x.size());
}

inline Eigen::VectorXi std_to_eigenvec(const std::vector<int> &x)
{
  return Eigen::Map<const Eigen::VectorXi>(&x.front(), x.size());
}

}; // namespace calin
