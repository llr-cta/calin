/*

   calin/unit_tests/math/test_pmt_ses_models.cpp -- Stephen Fegan -- 2017-04-24

   Unit tests for pmt_ses_models classes

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>
#include <vector>

#include "Eigen/Dense"
#include "calib/pmt_ses_models.hpp"

using namespace calin::calib::pmt_ses_models;
using namespace calin::math::function;

TEST(TestTwoGaussianSES, GradientCheck) {
  Eigen::VectorXd p(4);
  p << 0.35, 3.5, 6.0, 1.2; // frac, sigma-lo, mean-hi, sigma-hi
  TwoGaussianSES pdf;
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(double x=0; x<20.0; x+=0.1)
  {
    double good;
    EXPECT_TRUE(gradient_check(pdf, x, 1e-6, good));
    EXPECT_LE(good, 0.5);
  }
}

TEST(TestTwoGaussianSES, HessianCheck) {
  Eigen::VectorXd p(4);
  p << 0.35, 3.5, 6.0, 1.2; // frac, sigma-lo, mean-hi, sigma-hi
  TwoGaussianSES pdf;
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(double x=0; x<20.0; x+=0.1)
  {
    double good;
    EXPECT_TRUE(hessian_check(pdf, x, 1e-6, good));
    EXPECT_LE(good, 0.5);
  }
}

TEST(TestTwoGaussianSES, ParameterGradientCheck) {
  Eigen::VectorXd p(4);
  p << 0.35, 3.5, 6.0, 1.2; // frac, sigma-lo, mean-hi, sigma-hi
  Eigen::VectorXd dp(4);
  dp << 1e-4, 1e-4, 1e-4, 1e-4;
  TwoGaussianSES pdf;
  for(double x=0; x<20.0; x+=0.1)
  {
    Eigen::VectorXd good(4);
    EXPECT_TRUE(gradient_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
    EXPECT_LE(good(2), 0.5);
    EXPECT_LE(good(3), 0.5);
  }
}

TEST(TestTwoGaussianSES, ParameterHessianCheck) {
  Eigen::VectorXd p(4);
  p << 0.35, 3.5, 6.0, 1.2; // frac, sigma-lo, mean-hi, sigma-hi
  Eigen::VectorXd dp(4);
  dp << 1e-4, 1e-4, 1e-4, 1e-4;
  TwoGaussianSES pdf;
  for(double x=0; x<20.0; x+=0.1)
  {
    Eigen::MatrixXd good(4,4);
    EXPECT_TRUE(hessian_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0,0), 0.5);
    EXPECT_LE(good(0,1), 0.5);
    EXPECT_LE(good(0,2), 0.5);
    EXPECT_LE(good(0,3), 0.5);
    EXPECT_LE(good(1,0), 0.5);
    EXPECT_LE(good(1,1), 0.5);
    EXPECT_LE(good(1,2), 0.5);
    EXPECT_LE(good(1,3), 0.5);
    EXPECT_LE(good(2,0), 0.5);
    EXPECT_LE(good(2,1), 0.5);
    EXPECT_LE(good(2,2), 0.5);
    EXPECT_LE(good(2,3), 0.5);
    EXPECT_LE(good(3,0), 0.5);
    EXPECT_LE(good(3,1), 0.5);
    EXPECT_LE(good(3,2), 0.5);
    EXPECT_LE(good(3,3), 0.5);
  }
}

TEST(TestReparameterizedTwoGaussianSES, VerifyGainAndResolution) {
  Eigen::VectorXd p(4);
  p << 0.176559, 0.362319, 100, 0.451781; // frac, beta_lo, gain, beta
  ReparameterizedTwoGaussianSES pdf;
  pdf.set_parameter_values(p);
  EXPECT_NEAR(p[2], pdf.underlying_ses()->gain(), 1e-5);
  EXPECT_NEAR(p[3], pdf.underlying_ses()->resolution(), 1e-5);

  p << 0.176559, 0.362319, 30, 0.451781; // frac, beta_lo, gain, beta
  pdf.set_parameter_values(p);
  EXPECT_NEAR(p[2], pdf.underlying_ses()->gain(), 1e-5);
  EXPECT_NEAR(p[3], pdf.underlying_ses()->resolution(), 1e-5);

  p << 0.176559, 0.362319, 100, 0.42; // frac, beta_lo, gain, beta
  pdf.set_parameter_values(p);
  EXPECT_NEAR(p[2], pdf.underlying_ses()->gain(), 1e-5);
  EXPECT_NEAR(p[3], pdf.underlying_ses()->resolution(), 1e-5);
}

TEST(TestReparameterizedTwoGaussianSES, ParameterGradientCheck) {
  Eigen::VectorXd p(4);
  p << 0.176559, 0.362319, 100, 0.451781; // frac, beta_lo, gain, beta
  Eigen::VectorXd dp(4);
  dp << 1e-5, 1e-5, 1e-3, 1e-5;
  ReparameterizedTwoGaussianSES pdf;
  for(double x=0; x<200.0; x+=1.0)
  {
    Eigen::VectorXd good(4);
    EXPECT_TRUE(gradient_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
    EXPECT_LE(good(2), 0.5);
    EXPECT_LE(good(3), 0.5);
  }
}

TEST(TestReparameterizedTwoGaussianSES, ParameterHessianCheck) {
  Eigen::VectorXd p(4);
  p << 0.176559, 0.362319, 100, 0.451781; // frac, beta_lo, gain, beta
  Eigen::VectorXd dp(4);
  dp << 1e-5, 1e-5, 1e-3, 1e-5;
  ReparameterizedTwoGaussianSES pdf;
  for(double x=1.0; x<200.0; x+=1.0)
  {
    Eigen::MatrixXd good(4,4);
    EXPECT_TRUE(hessian_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0,0), 0.5);
    EXPECT_LE(good(0,1), 0.5);
    EXPECT_LE(good(0,2), 0.5);
    EXPECT_LE(good(0,3), 0.5);
    EXPECT_LE(good(1,0), 0.5);
    EXPECT_LE(good(1,1), 0.5);
    EXPECT_LE(good(1,2), 0.5);
    EXPECT_LE(good(1,3), 0.5);
    EXPECT_LE(good(2,0), 0.5);
    EXPECT_LE(good(2,1), 0.5);
    EXPECT_LE(good(2,2), 0.5);
    EXPECT_LE(good(2,3), 0.5);
    EXPECT_LE(good(3,0), 0.5);
    EXPECT_LE(good(3,1), 0.5);
    EXPECT_LE(good(3,2), 0.5);
    EXPECT_LE(good(3,3), 0.5);
  }
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
