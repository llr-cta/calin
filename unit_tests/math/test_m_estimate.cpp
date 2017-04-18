/*

   calin/unit_tests/math/test_m_estimate.cpp -- Stephen Fegan -- 2017-04-18

   Unit tests for m_estimate

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>
#include <vector>

#include "Eigen/Dense"
#include "math/m_estimate.hpp"

using namespace calin::math::m_estimate;

namespace {
  constexpr double maxgood = 1.0;
}

TEST(TestHyperbolicLikelihoodRhoFunction, GradientCheck) {
  for(double C = 0; C < 10; C += 1) {
    for(double D = 0.1; D < 2; D += 0.1) {
      HyperbolicLikelihoodRhoFunction rho(C,D);
      Eigen::VectorXd x(1);
      Eigen::VectorXd dx(1);
      dx << 0.001;
      for(x(0) = -10.0; x(0)<20.0; x(0)+=0.1) {
        Eigen::VectorXd good(1);
        EXPECT_TRUE(gradient_check(rho, x, dx, good, maxgood));
        EXPECT_LE(good(0), maxgood);
      }
    }
  }
}

TEST(TestHyperbolicLikelihoodRhoFunction, HessianCheck) {
  for(double C = 0; C < 10; C += 1) {
    for(double D = 0.1; D < 2; D += 0.1) {
      HyperbolicLikelihoodRhoFunction rho(C,D);
      Eigen::VectorXd x(1);
      Eigen::VectorXd dx(1);
      dx << 0.001;
      for(x(0) = -10.0; x(0)<20.0; x(0)+=0.1) {
        Eigen::MatrixXd good(1,1);
        EXPECT_TRUE(hessian_check(rho, x, dx, good, maxgood));
        EXPECT_LE(good(0,0), maxgood);
      }
    }
  }
}

#if 0
TEST(TestLikelihood, HessianCheck) {
  SimpleHist hist(1.0);
  RNG rng(12345);
  for(unsigned i=0; i<10000; i++)
    hist.insert(rng.normal()*10.0 + 100.0);
  BinnedGaussianPDF gaussian(1.0);
  IID1DDataLikelihoodFunction cost_fn(&gaussian, hist);
  Eigen::VectorXd x(2);
  Eigen::VectorXd dx(2);
  dx << 0.00001, 0.000001;
  for(x(0) = 80.0; x(0)<120.0; x(0)+=1.0) {
    for(x(1) = 8.0; x(1)<12.0; x(1)+=0.1) {
      Eigen::MatrixXd good(2,2);
      EXPECT_TRUE(hessian_check(cost_fn, x, dx, good, maxgood));
      EXPECT_LE(good(0,0), maxgood);
      EXPECT_LE(good(0,1), maxgood);
      EXPECT_LE(good(1,0), maxgood);
      EXPECT_LE(good(1,1), maxgood);
    }
  }
}
#endif

TEST(TestModifiedHyperbolicLikelihoodRhoFunction, GradientCheck) {
  for(double C = 1; C < 10; C += 1) {
    for(double D = 0.1; D < std::min(C,2.0); D += 0.1) {
      ModifiedHyperbolicLikelihoodRhoFunction rho(C,D);
      Eigen::VectorXd x(1);
      Eigen::VectorXd dx(1);
      dx << 0.001;
      for(x(0) = -10.0; x(0)<20.0; x(0)+=0.1) {
        //std::cout << C << ' ' << D << ' ' << x(0) << '\n';
        Eigen::VectorXd good(1);
        EXPECT_TRUE(gradient_check(rho, x, dx, good, maxgood));
        EXPECT_LE(good(0), maxgood);
      }
    }
  }
}

TEST(TestModifiedHyperbolicLikelihoodRhoFunction, HessianCheck) {
  for(double C = 1; C < 10; C += 1) {
    for(double D = 0.1; D < std::min(C,2.0); D += 0.1) {
      ModifiedHyperbolicLikelihoodRhoFunction rho(C,D);
      Eigen::VectorXd x(1);
      Eigen::VectorXd dx(1);
      dx << 0.001;
      for(x(0) = -10.0; x(0)<20.0; x(0)+=0.1) {
        //std::cout << C << ' ' << D << ' ' << x(0) << '\n';
        Eigen::MatrixXd good(1,1);
        EXPECT_TRUE(hessian_check(rho, x, dx, good, maxgood));
        EXPECT_LE(good(0,0), maxgood);
      }
    }
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
