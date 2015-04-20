/* 

   calin/package_wide_definitions.hpp -- Stephen Fegan -- 2015-04-16

   Definitions of types and macros that apply to the full package

*/

#pragma once

#include <Eigen/Core>

namespace calin { 

#ifndef SWIG
#define CALIN_TYPEALIAS(A,B) using A = B
#else
#define CALIN_TYPEALIAS(A,B) typedef B A
#endif

//#define CALIN_USE_EIGEN_REF

#ifdef CALIN_USE_EIGEN_REF
CALIN_TYPEALIAS(VecRef, Eigen::Ref<Eigen::VectorXd>);
CALIN_TYPEALIAS(ConstVecRef, const Eigen::Ref<const Eigen::VectorXd>&);
CALIN_TYPEALIAS(MatRef, Eigen::Ref<Eigen::MatrixXd>);
CALIN_TYPEALIAS(ConstMatRef, const Eigen::Ref<const Eigen::MatrixXd>&);
#else
CALIN_TYPEALIAS(VecRef, Eigen::VectorXd&);
CALIN_TYPEALIAS(MatRef, Eigen::MatrixXd&);
CALIN_TYPEALIAS(ConstVecRef, const Eigen::VectorXd&);
CALIN_TYPEALIAS(ConstMatRef, const Eigen::MatrixXd&);
#endif

}; // namespace calin
