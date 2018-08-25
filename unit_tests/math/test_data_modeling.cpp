/*

   calin/unit_tests/math/test_data_modeling.cpp -- Stephen Fegan -- 2017-04-18

   Unit tests for data_modeling

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

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
#include "math/histogram.hpp"
#include "math/rng.hpp"
#include "math/pdf_1d.hpp"
#include "math/data_modeling.hpp"
#include "math/m_estimate.hpp"

using namespace calin::math::data_modeling;
using namespace calin::math::rng;
using namespace calin::math::pdf_1d;
using namespace calin::math::histogram;
using namespace calin::math::m_estimate;

namespace {
  constexpr double maxgood = 1.0;
}

TEST(TestLikelihood, GradientCheck) {
  SimpleHist hist(1.0);
  RNG rng(RNG::std_test_seed);
  for(unsigned i=0; i<10000; i++)
    hist.insert(rng.normal()*10.0 + 100.0);
  BinnedGaussianPDF gaussian(1.0);
  IID1DDataLikelihoodFunction cost_fn(&gaussian, hist);
  Eigen::VectorXd x(2);
  Eigen::VectorXd dx(2);
  dx << 0.01, 0.001;
  for(x(0) = 80.0; x(0)<120.0; x(0)+=1.0) {
    for(x(1) = 8.0; x(1)<12.0; x(1)+=0.1) {
      Eigen::VectorXd good(2);
      EXPECT_TRUE(gradient_check(cost_fn, x, dx, good, maxgood));
      EXPECT_LE(good(0), maxgood);
      EXPECT_LE(good(1), maxgood);
    }
  }
}

TEST(TestLikelihood, HessianCheck) {
  SimpleHist hist(1.0);
  RNG rng(RNG::std_test_seed);
  for(unsigned i=0; i<10000; i++)
    hist.insert(rng.normal()*10.0 + 100.0);
  BinnedGaussianPDF gaussian(1.0);
  IID1DDataLikelihoodFunction cost_fn(&gaussian, hist);
  Eigen::VectorXd x(2);
  Eigen::VectorXd dx(2);
  dx << 0.01, 0.001;
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

TEST(TestChi2, GradientCheck) {
  SimpleHist hist(1.0);
  RNG rng(RNG::std_test_seed);
  for(unsigned i=0; i<10000; i++)
    hist.insert(rng.normal()*10.0 + 100.0);
  BinnedGaussianPDF gaussian(1.0);
  IID1DDataChi2Function cost_fn(&gaussian, hist);
  Eigen::VectorXd x(3);
  Eigen::VectorXd dx(3);
  dx << 0.1, 0.00001, 0.000001;
  for(x(1) = 80.0; x(1)<120.0; x(1)+=1.0) {
    for(x(2) = 8.0; x(2)<12.0; x(2)+=0.1) {
      x(0) = 10000 + (x(1)-100.0)*10.0;
      Eigen::VectorXd good(2);
      EXPECT_TRUE(gradient_check(cost_fn, x, dx, good, maxgood));
      EXPECT_LE(good(0), maxgood);
      EXPECT_LE(good(1), maxgood);
      EXPECT_LE(good(2), maxgood);
    }
  }
}

TEST(TestChi2, HessianCheck) {
  SimpleHist hist(1.0);
  RNG rng(RNG::std_test_seed);
  for(unsigned i=0; i<10000; i++)
    hist.insert(rng.normal()*10.0 + 100.0);
  BinnedGaussianPDF gaussian(1.0);
  IID1DDataChi2Function cost_fn(&gaussian, hist);
  Eigen::VectorXd x(3);
  Eigen::VectorXd dx(3);
  dx << 0.1, 0.00001, 0.000001;
  for(x(1) = 80.0; x(1)<120.0; x(1)+=1.0) {
    for(x(2) = 8.0; x(2)<12.0; x(2)+=0.1) {
      x(0) = 10000 + (x(1)-100.0)*10.0;
      Eigen::MatrixXd good(3,3);
      EXPECT_TRUE(hessian_check(cost_fn, x, dx, good, 2.0));
      EXPECT_LE(good(0,0), maxgood);
      EXPECT_LE(good(0,1), maxgood);
      EXPECT_LE(good(0,2), maxgood);
      EXPECT_LE(good(1,0), maxgood);
      EXPECT_LE(good(1,1), maxgood);
      EXPECT_LE(good(1,2), maxgood);
      EXPECT_LE(good(2,0), maxgood);
      EXPECT_LE(good(2,1), maxgood);
      EXPECT_LE(good(2,2), maxgood);
    }
  }
}

TEST(TestMEstimateWithNULL, GradientCheck) {
  SimpleHist hist(1.0);
  RNG rng(RNG::std_test_seed);
  for(unsigned i=0; i<10000; i++)
    hist.insert(rng.normal()*10.0 + 100.0);
  BinnedGaussianPDF gaussian(1.0);
  NullLikelihoodRhoFunction rho;
  IID1DDataMEstimateLikelihoodFunction cost_fn(&gaussian, &rho, hist);
  Eigen::VectorXd x(2);
  Eigen::VectorXd dx(2);
  dx << 0.01, 0.001;
  for(x(0) = 80.0; x(0)<120.0; x(0)+=1.0) {
    for(x(1) = 8.0; x(1)<12.0; x(1)+=0.1) {
      Eigen::VectorXd good(2);
      EXPECT_TRUE(gradient_check(cost_fn, x, dx, good, maxgood));
      EXPECT_LE(good(0), maxgood);
      EXPECT_LE(good(1), maxgood);
    }
  }
}

TEST(TestMEstimateWithNULL, HessianCheck) {
  SimpleHist hist(1.0);
  RNG rng(RNG::std_test_seed);
  for(unsigned i=0; i<10000; i++)
    hist.insert(rng.normal()*10.0 + 100.0);
  BinnedGaussianPDF gaussian(1.0);
  NullLikelihoodRhoFunction rho;
  IID1DDataMEstimateLikelihoodFunction cost_fn(&gaussian, &rho, hist);
  Eigen::VectorXd x(2);
  Eigen::VectorXd dx(2);
  dx << 0.01, 0.001;
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

TEST(TestMEstimateWithHyperbolic, GradientCheck) {
  SimpleHist hist(1.0);
  RNG rng(RNG::std_test_seed);
  for(unsigned i=0; i<10000; i++)
    hist.insert(rng.normal()*10.0 + 100.0);
  BinnedGaussianPDF gaussian(1.0);
  for(double C = 0; C < 10; C += 2) {
    for(double D = 0.1; D < 2; D += 0.2) {
      HyperbolicLikelihoodRhoFunction rho(C,D);
      IID1DDataMEstimateLikelihoodFunction cost_fn(&gaussian, &rho, hist);
      Eigen::VectorXd x(2);
      Eigen::VectorXd dx(2);
      dx << 0.01, 0.001;
      for(x(0) = 80.0; x(0)<120.0; x(0)+=2.0) {
        for(x(1) = 8.0; x(1)<12.0; x(1)+=0.2) {
          Eigen::VectorXd good(2);
          EXPECT_TRUE(gradient_check(cost_fn, x, dx, good, maxgood));
          EXPECT_LE(good(0), maxgood);
          EXPECT_LE(good(1), maxgood);
        }
      }
    }
  }
}

TEST(TestMEstimateWithHyperbolic, HessianCheck) {
  SimpleHist hist(1.0);
  RNG rng(RNG::std_test_seed);
  for(unsigned i=0; i<10000; i++)
    hist.insert(rng.normal()*10.0 + 100.0);
  BinnedGaussianPDF gaussian(1.0);
  for(double C = 0; C < 10; C += 2) {
    for(double D = 0.1; D < 2; D += 0.2) {
      HyperbolicLikelihoodRhoFunction rho(C,D);
      IID1DDataMEstimateLikelihoodFunction cost_fn(&gaussian, &rho, hist);
      Eigen::VectorXd x(2);
      Eigen::VectorXd dx(2);
      dx << 0.01, 0.001;
      for(x(0) = 80.0; x(0)<120.0; x(0)+=2.0) {
        for(x(1) = 8.0; x(1)<12.0; x(1)+=0.2) {
          Eigen::MatrixXd good(2,2);
          EXPECT_TRUE(hessian_check(cost_fn, x, dx, good, maxgood));
          EXPECT_LE(good(0,0), maxgood);
          EXPECT_LE(good(0,1), maxgood);
          EXPECT_LE(good(1,0), maxgood);
          EXPECT_LE(good(1,1), maxgood);
        }
      }
    }
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
