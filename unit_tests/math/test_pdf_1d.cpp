/*

   calin/unit_tests/math/test_pdf_1d.cpp -- Stephen Fegan -- 2015-04-03

   Unit tests for pdf_1d classes

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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
#include "math/pdf_1d.hpp"
#include "math/log_quadratic_spline_pdf_1d.hpp"

using namespace calin::math;
using namespace calin::math::function;
using namespace calin::math::pdf_1d;

TEST(TestGaussianPDF, GradientCheck) {
  Eigen::VectorXd p(2);
  p << 3, 2; // mean, sigma
  GaussianPDF pdf;
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(double x = 0; x<10.0; x+=0.1)
  {
    double good;
    EXPECT_TRUE(gradient_check(pdf, x, 1e-6, good));
    EXPECT_LE(good, 0.5);
  }
}

TEST(TestGaussianPDF, HessianCheck) {
  Eigen::VectorXd p(2);
  p << 3, 2; // mean, sigma
  GaussianPDF pdf;
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(double x = 0; x<10.0; x+=0.1)
  {
    double good;
    EXPECT_TRUE(hessian_check(pdf, x, 1e-6, good));
    EXPECT_LE(good, 0.5);
  }
}

TEST(TestGaussianPDF, ParameterGradientCheck) {
  Eigen::VectorXd p(2);
  p << 3, 2; // mean, sigma
  Eigen::VectorXd dp(2);
  dp << 1e-6, 1e-6;
  GaussianPDF pdf;
  for(double x = 0; x<10.0; x+=0.1)
  {
    Eigen::VectorXd good(2);
    EXPECT_TRUE(gradient_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
  }
}

TEST(TestGaussianPDF, ParameterHessianCheck) {
  Eigen::VectorXd p(2);
  p << 3, 2; // mean, sigma
  Eigen::VectorXd dp(2);
  dp << 1e-6, 1e-6;
  GaussianPDF pdf;
  for(double x = 0; x<10.0; x+=0.1)
  {
    Eigen::MatrixXd good(2,2);
    EXPECT_TRUE(hessian_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0,0), 0.5);
    EXPECT_LE(good(0,1), 0.5);
    EXPECT_LE(good(1,0), 0.5);
    EXPECT_LE(good(1,1), 0.5);
  }
}

TEST(TestBinnedGaussianPDF, NormalisedValues) {
  Eigen::VectorXd p(2);
  p << 3, 0.25; // mean, sigma
  Eigen::VectorXd dp(2);
  dp << 1e-6, 1e-6;
  double dx = 0.5;
  BinnedGaussianPDF pdf(dx);
  pdf.set_parameter_values(p);
  double sum = 0;
  for(double x=-10; x<13; x+=dx)sum += pdf.value_1d(x);
  EXPECT_NEAR(sum*dx, 1.0, dx);
}

TEST(TestBinnedGaussianPDF, GradientCheck) {
  Eigen::VectorXd p(2);
  p << 3, 0.25; // mean, sigma
  BinnedGaussianPDF pdf(0.5);
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(double x = 0; x<10.0; x+=0.1)
  {
    double good;
    EXPECT_TRUE(gradient_check(pdf, x, 1e-6, good));
    EXPECT_LE(good, 0.5);
  }
}

TEST(TestBinnedGaussianPDF, HessianCheck) {
  Eigen::VectorXd p(2);
  p << 3, 0.25; // mean, sigma
  BinnedGaussianPDF pdf(0.5);
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(double x = 0; x<10.0; x+=0.1)
  {
    double good;
    EXPECT_TRUE(hessian_check(pdf, x, 1e-6, good));
    EXPECT_LE(good, 0.5);
  }
}

TEST(TestBinnedGaussianPDF, ParameterGradientCheck) {
  Eigen::VectorXd p(2);
  p << 3, 0.25; // mean, sigma
  Eigen::VectorXd dp(2);
  dp << 1e-6, 1e-6;
  BinnedGaussianPDF pdf(0.5);
  for(double x = 0; x<10.0; x+=0.1)
  {
    Eigen::VectorXd good(2);
    EXPECT_TRUE(gradient_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
  }
}

TEST(TestBinnedGaussianPDF, ParameterHessianCheck) {
  Eigen::VectorXd p(2);
  p << 3, 0.25; // mean, sigma
  Eigen::VectorXd dp(2);
  dp << 1e-6, 1e-6;
  BinnedGaussianPDF pdf(0.5);
  for(double x = 0; x<10.0; x+=0.1)
  {
    Eigen::MatrixXd good(2,2);
    EXPECT_TRUE(hessian_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0,0), 0.5);
    EXPECT_LE(good(0,1), 0.5);
    EXPECT_LE(good(1,0), 0.5);
    EXPECT_LE(good(1,1), 0.5);
  }
}

TEST(TestLimitedGaussianPDF, GradientCheck) {
  Eigen::VectorXd p(2);
  p << 3, 2; // mean, sigma
  LimitedGaussianPDF pdf(1.0,9.0);
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(double x = 0; x<10.0; x+=0.1)
  {
    double good;
    EXPECT_TRUE(gradient_check(pdf, x, 1e-6, good));
    EXPECT_LE(good, 0.5);
  }
}

TEST(TestLimitedGaussianPDF, HessianCheck) {
  Eigen::VectorXd p(2);
  p << 3, 2; // mean, sigma
  LimitedGaussianPDF pdf(1.0,9.0);
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(double x = 0; x<10.0; x+=0.1)
  {
    double good;
    EXPECT_TRUE(hessian_check(pdf, x, 1e-6, good));
    EXPECT_LE(good, 0.5);
  }
}

TEST(TestLimitedGaussianPDF, ParameterGradientCheck) {
  Eigen::VectorXd p(2);
  p << 3, 2; // mean, sigma
  Eigen::VectorXd dp(2);
  dp << 1e-6, 1e-6;
  LimitedGaussianPDF pdf(1.0,9.0);
  for(double x = 0; x<10.0; x+=0.1)
  {
    Eigen::VectorXd good(2);
    EXPECT_TRUE(gradient_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
  }
}

TEST(TestLimitedGaussianPDF, ParameterHessianCheck) {
  Eigen::VectorXd p(2);
  p << 3, 2; // mean, sigma
  Eigen::VectorXd dp(2);
  dp << 1e-6, 1e-6;
  LimitedGaussianPDF pdf(1.0,9.0);
  for(double x = 0; x<10.0; x+=0.1)
  {
    Eigen::MatrixXd good(2,2);
    EXPECT_TRUE(hessian_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0,0), 0.5);
    EXPECT_LE(good(0,1), 0.5);
    EXPECT_LE(good(1,0), 0.5);
    EXPECT_LE(good(1,1), 0.5);
  }
}

TEST(TestLimitedGaussianPDF_DX, GradientCheck) {
  Eigen::VectorXd p(2);
  p << 3, 2; // mean, sigma
  LimitedGaussianPDF pdf(1.0,9.0,0.5);
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(double x = 0; x<10.0; x+=0.1)
  {
    double good;
    EXPECT_TRUE(gradient_check(pdf, x, 1e-4, good));
    EXPECT_LE(good, 0.5);
  }
}

TEST(TestLimitedGaussianPDF_DX, HessianCheck) {
  Eigen::VectorXd p(2);
  p << 3, 2; // mean, sigma
  LimitedGaussianPDF pdf(1.0,9.0,0.5);
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(double x = 0; x<10.0; x+=0.1)
  {
    double good;
    EXPECT_TRUE(hessian_check(pdf, x, 1e-4, good));
    EXPECT_LE(good, 0.5);
  }
}

TEST(TestLimitedGaussianPDF_DX, ParameterGradientCheck) {
  Eigen::VectorXd p(2);
  p << 3, 2; // mean, sigma
  Eigen::VectorXd dp(2);
  dp << 1e-6, 1e-6;
  LimitedGaussianPDF pdf(1.0,9.0,0.5);
  for(double x = 0; x<10.0; x+=0.1)
  {
    Eigen::VectorXd good(2);
    EXPECT_TRUE(gradient_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
  }
}

TEST(TestLimitedGaussianPDF_DX, ParameterHessianCheck) {
  Eigen::VectorXd p(2);
  p << 3, 2; // mean, sigma
  Eigen::VectorXd dp(2);
  dp << 1e-6, 1e-6;
  LimitedGaussianPDF pdf(1.0,9.0,0.5);
  for(double x = 0; x<10.0; x+=0.1)
  {
    Eigen::MatrixXd good(2,2);
    EXPECT_TRUE(hessian_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0,0), 0.5);
    EXPECT_LE(good(0,1), 0.5);
    EXPECT_LE(good(1,0), 0.5);
    EXPECT_LE(good(1,1), 0.5);
  }
}

TEST(TestLimitedExponentialPDF, GradientCheck) {
  Eigen::VectorXd p(1);
  p << 3; // scale
  LimitedExponentialPDF pdf(1.0,9.0);
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(double x = 0; x<10.0; x+=0.1)
  {
    double good;
    EXPECT_TRUE(gradient_check(pdf, x, 1e-6, good));
    EXPECT_LE(good, 0.5);
  }
}

TEST(TestLimitedExponentialPDF, HessianCheck) {
  Eigen::VectorXd p(1);
  p << 3; // scale
  LimitedExponentialPDF pdf(1.0,9.0);
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(double x = 0; x<10.0; x+=0.1)
  {
    double good;
    EXPECT_TRUE(hessian_check(pdf, x, 1e-6, good));
    EXPECT_LE(good, 0.5);
  }
}

TEST(TestLimitedExponentialPDF, ParameterGradientCheck) {
  Eigen::VectorXd p(1);
  p << 3; // scale
  Eigen::VectorXd dp(1);
  dp << 1e-6;
  LimitedExponentialPDF pdf(1.0,9.0);
  for(double x = 0; x<10.0; x+=0.1)
  {
    Eigen::VectorXd good(1);
    EXPECT_TRUE(gradient_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0), 0.5);
  }
}

TEST(TestLimitedExponentialPDF, ParameterHessianCheck) {
  Eigen::VectorXd p(1);
  p << 3; // scale
  Eigen::VectorXd dp(1);
  dp << 1e-6;
  LimitedExponentialPDF pdf(1.0,9.0);
  for(double x = 0; x<10.0; x+=0.1)
  {
    Eigen::MatrixXd good(1,1);
    EXPECT_TRUE(hessian_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0,0), 0.5);
  }
}

TEST(TestLimitedExponentialPDF, Binned_ParameterGradientCheck) {
  Eigen::VectorXd p(1);
  p << 3; // scale
  Eigen::VectorXd dp(1);
  dp << 1e-6;
  LimitedExponentialPDF pdf(1.0,9.0,0.1);
  for(double x = 0; x<10.0; x+=0.1)
  {
    Eigen::VectorXd good(1);
    EXPECT_TRUE(gradient_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0), 0.5);
  }
}

TEST(TestLimitedExponentialPDF, Binned_ParameterHessianCheck) {
  Eigen::VectorXd p(1);
  p << 3; // scale
  Eigen::VectorXd dp(1);
  dp << 1e-6;
  LimitedExponentialPDF pdf(1.0,9.0,0.1);
  for(double x = 0; x<10.0; x+=0.1)
  {
    Eigen::MatrixXd good(1,1);
    EXPECT_TRUE(hessian_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0,0), 0.5);
  }
}

TEST(TestTwoComponent1DPDF_ExpGauss, SetAndRecallParameters) {
  pdf_1d::LimitedExponentialPDF exp_pdf(1.0, 9.0);
  pdf_1d::LimitedGaussianPDF gauss_pdf(1.0, 9.0);
  pdf_1d::TwoComponent1DPDF pdf(&exp_pdf, "exp", &gauss_pdf, "gauss");

  EXPECT_EQ(pdf.num_parameters(), 4U);
  auto pvec = pdf.parameters();
  EXPECT_EQ(pvec[0].name, "exp_probability");
  EXPECT_EQ(pvec[1].name, "exp.scale");
  EXPECT_EQ(pvec[2].name, "gauss.mean");
  EXPECT_EQ(pvec[3].name, "gauss.rms");
  Eigen::VectorXd p(4);
  p << 0.5, 40.0, 100.0, 35.0;
  pdf.set_parameter_values(p);
  EXPECT_EQ(pdf.parameter_values(), p);
}

TEST(TestTwoComponent1DPDF_ExpGauss, GradientCheck) {
  pdf_1d::LimitedExponentialPDF exp_pdf(1.0, 9.0);
  pdf_1d::LimitedGaussianPDF gauss_pdf(1.0, 9.0);
  pdf_1d::TwoComponent1DPDF pdf(&exp_pdf, "exp", &gauss_pdf, "gauss");
  Eigen::VectorXd p(4);
  p << 0.5, 40.0, 100.0, 35.0;
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(double x = 0; x<10.0; x+=0.1)
  {
    double good;
    EXPECT_TRUE(gradient_check(pdf, x, 1e-6, good));
    EXPECT_LE(good, 0.5);
  }
}

TEST(TestTwoComponent1DPDF_ExpGauss, HessianCheck) {
  pdf_1d::LimitedExponentialPDF exp_pdf(1.0, 9.0);
  pdf_1d::LimitedGaussianPDF gauss_pdf(1.0, 9.0);
  pdf_1d::TwoComponent1DPDF pdf(&exp_pdf, "exp", &gauss_pdf, "gauss");
  Eigen::VectorXd p(4);
  p << 0.5, 40.0, 100.0, 35.0;
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(double x = 0; x<10.0; x+=0.1)
  {
    double good;
    EXPECT_TRUE(hessian_check(pdf, x, 1e-6, good));
    EXPECT_LE(good, 0.5);
  }
}

TEST(TestTwoComponent1DPDF_ExpGauss, ParameterGradientCheck) {
  pdf_1d::LimitedExponentialPDF exp_pdf(1.0, 9.0);
  pdf_1d::LimitedGaussianPDF gauss_pdf(1.0, 9.0);
  pdf_1d::TwoComponent1DPDF pdf(&exp_pdf, "exp", &gauss_pdf, "gauss");
  Eigen::VectorXd p(4);
  p << 0.5, 40.0, 100.0, 35.0;
  Eigen::VectorXd dp(4);
  dp << 1e-6, 1e-6, 1e-6, 1e-6;
  for(double x = 0; x<10.0; x+=0.1)
  {
    Eigen::VectorXd good(4);
    EXPECT_TRUE(gradient_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
    EXPECT_LE(good(2), 0.5);
    EXPECT_LE(good(3), 0.5);
  }
}

TEST(TestTwoComponent1DPDF_ExpGauss, ParameterHessianCheck) {
  pdf_1d::LimitedExponentialPDF exp_pdf(1.0, 9.0);
  pdf_1d::LimitedGaussianPDF gauss_pdf(1.0, 9.0);
  pdf_1d::TwoComponent1DPDF pdf(&exp_pdf, "exp", &gauss_pdf, "gauss");
  Eigen::VectorXd p(4);
  p << 0.5, 40.0, 100.0, 35.0;
  Eigen::VectorXd dp(4);
  dp << 1e-6, 1e-6, 1e-6, 1e-6;
  for(double x = 0; x<10.0; x+=0.1)
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

TEST(TestFreezeThaw, FreezeAndThaw) {
  pdf_1d::LimitedExponentialPDF exp_pdf(1.0, 9.0);
  pdf_1d::LimitedGaussianPDF gauss_pdf(1.0, 9.0);
  pdf_1d::TwoComponent1DPDF pdf_x(&exp_pdf, "exp", &gauss_pdf, "gauss");
  SingleToParameterizableMultiAxisFunctionAdapter maf(&pdf_x);
  function::PMAFReverser pdf(&maf);
  function::FreezeThawFunction freezer(&pdf);
  EXPECT_EQ(freezer.num_domain_axes(), 4U);
  EXPECT_EQ(freezer.num_parameters(), 0U);
  EXPECT_EQ(freezer.num_domain_axes(), freezer.domain_axes().size());
  EXPECT_EQ(freezer.num_parameters(), freezer.parameters().size());

  EXPECT_TRUE(freezer.freeze(1,50.0));
  EXPECT_EQ(freezer.num_domain_axes(), 3U);
  EXPECT_EQ(freezer.num_parameters(), 1U);
  EXPECT_EQ(freezer.num_domain_axes(), freezer.domain_axes().size());
  EXPECT_EQ(freezer.num_parameters(), freezer.parameters().size());
  EXPECT_EQ(freezer.free_axes(), std::vector<unsigned>({0,2,3}));
  EXPECT_EQ(freezer.frozen_axes(), std::vector<unsigned>({1}));
  Eigen::VectorXd p1(1);
  p1 << 50.0;
  EXPECT_EQ(freezer.parameter_values(), p1);

  EXPECT_TRUE(freezer.freeze(0,0.2));
  EXPECT_EQ(freezer.num_domain_axes(), 2U);
  EXPECT_EQ(freezer.num_parameters(), 2U);
  EXPECT_EQ(freezer.num_domain_axes(), freezer.domain_axes().size());
  EXPECT_EQ(freezer.num_parameters(), freezer.parameters().size());
  EXPECT_EQ(freezer.free_axes(), std::vector<unsigned>({2,3}));
  EXPECT_EQ(freezer.frozen_axes(), std::vector<unsigned>({0,1}));
  Eigen::VectorXd p01(2);
  p01 << 0.2, 50.0;
  EXPECT_EQ(freezer.parameter_values(), p01);

  EXPECT_TRUE(freezer.freeze(3,35.0));
  EXPECT_EQ(freezer.num_domain_axes(), 1U);
  EXPECT_EQ(freezer.num_parameters(), 3U);
  EXPECT_EQ(freezer.num_domain_axes(), freezer.domain_axes().size());
  EXPECT_EQ(freezer.num_parameters(), freezer.parameters().size());
  EXPECT_EQ(freezer.free_axes(), std::vector<unsigned>({2}));
  EXPECT_EQ(freezer.frozen_axes(), std::vector<unsigned>({0,1,3}));
  Eigen::VectorXd p013(3);
  p013 << 0.2, 50.0, 35.0;
  EXPECT_EQ(freezer.parameter_values(), p013);

  EXPECT_FALSE(freezer.freeze(3,35.0));
  EXPECT_EQ(freezer.num_domain_axes(), 1U);
  EXPECT_EQ(freezer.num_parameters(), 3U);
  EXPECT_EQ(freezer.num_domain_axes(), freezer.domain_axes().size());
  EXPECT_EQ(freezer.num_parameters(), freezer.parameters().size());
  EXPECT_EQ(freezer.free_axes(), std::vector<unsigned>({2}));
  EXPECT_EQ(freezer.frozen_axes(), std::vector<unsigned>({0,1,3}));
  EXPECT_EQ(freezer.parameter_values(), p013);

  EXPECT_TRUE(freezer.freeze(2,100.0));
  EXPECT_EQ(freezer.num_domain_axes(), 0U);
  EXPECT_EQ(freezer.num_parameters(), 4U);
  EXPECT_EQ(freezer.num_domain_axes(), freezer.domain_axes().size());
  EXPECT_EQ(freezer.num_parameters(), freezer.parameters().size());
  EXPECT_EQ(freezer.free_axes(), std::vector<unsigned>({}));
  EXPECT_EQ(freezer.frozen_axes(), std::vector<unsigned>({0,1,2,3}));
  Eigen::VectorXd p0123(4);
  p0123 << 0.2, 50.0, 100.0, 35.0;
  EXPECT_EQ(freezer.parameter_values(), p0123);

  EXPECT_TRUE(freezer.thaw(1));
  EXPECT_EQ(freezer.num_domain_axes(), 1U);
  EXPECT_EQ(freezer.num_parameters(), 3U);
  EXPECT_EQ(freezer.num_domain_axes(), freezer.domain_axes().size());
  EXPECT_EQ(freezer.num_parameters(), freezer.parameters().size());
  EXPECT_EQ(freezer.free_axes(), std::vector<unsigned>({1}));
  EXPECT_EQ(freezer.frozen_axes(), std::vector<unsigned>({0,2,3}));
  Eigen::VectorXd p023(3);
  p023 << 0.2, 100.0, 35.0;
  EXPECT_EQ(freezer.parameter_values(), p023);

  EXPECT_FALSE(freezer.thaw(1));
  EXPECT_EQ(freezer.num_domain_axes(), 1U);
  EXPECT_EQ(freezer.num_parameters(), 3U);
  EXPECT_EQ(freezer.num_domain_axes(), freezer.domain_axes().size());
  EXPECT_EQ(freezer.num_parameters(), freezer.parameters().size());
  EXPECT_EQ(freezer.free_axes(), std::vector<unsigned>({1}));
  EXPECT_EQ(freezer.frozen_axes(), std::vector<unsigned>({0,2,3}));
  EXPECT_EQ(freezer.parameter_values(), p023);

  EXPECT_TRUE(freezer.thaw(2));
  EXPECT_EQ(freezer.num_domain_axes(), 2U);
  EXPECT_EQ(freezer.num_parameters(), 2U);
  EXPECT_EQ(freezer.num_domain_axes(), freezer.domain_axes().size());
  EXPECT_EQ(freezer.num_parameters(), freezer.parameters().size());
  EXPECT_EQ(freezer.free_axes(), std::vector<unsigned>({1,2}));
  EXPECT_EQ(freezer.frozen_axes(), std::vector<unsigned>({0,3}));
  Eigen::VectorXd p03(2);
  p03 << 0.2, 35.0;
  EXPECT_EQ(freezer.parameter_values(), p03);

  EXPECT_TRUE(freezer.thaw(0));
  EXPECT_EQ(freezer.num_domain_axes(), 3U);
  EXPECT_EQ(freezer.num_parameters(), 1U);
  EXPECT_EQ(freezer.num_domain_axes(), freezer.domain_axes().size());
  EXPECT_EQ(freezer.num_parameters(), freezer.parameters().size());
  EXPECT_EQ(freezer.free_axes(), std::vector<unsigned>({0,1,2}));
  EXPECT_EQ(freezer.frozen_axes(), std::vector<unsigned>({3}));
  Eigen::VectorXd p3(1);
  p3 << 35.0;
  EXPECT_EQ(freezer.parameter_values(), p3);
}

TEST(TestFreezeThaw, GradientCheck) {
  pdf_1d::LimitedExponentialPDF exp_pdf(1.0, 9.0);
  pdf_1d::LimitedGaussianPDF gauss_pdf(1.0, 9.0);
  pdf_1d::TwoComponent1DPDF pdf_x(&exp_pdf, "exp", &gauss_pdf, "gauss");
  Eigen::VectorXd p(3);
  p << 0.5, 100.0, 35.0;
  Eigen::VectorXd dp(3);
  dp << 1e-6, 1e-6, 1e-6;
  SingleToParameterizableMultiAxisFunctionAdapter maf(&pdf_x);
  function::PMAFReverser pdf(&maf);
  function::FreezeThawFunction freezer(&pdf);
  freezer.freeze(1,40.0);
  Eigen::VectorXd x(1);
  for(x(0) = 0; x(0)<10.0; x(0)+=0.1)
  {
    pdf.set_parameter_values(x);
    EXPECT_EQ(x, pdf.parameter_values());
    Eigen::VectorXd good(3);
    EXPECT_TRUE(gradient_check(freezer, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
    EXPECT_LE(good(2), 0.5);
  }
}

TEST(TestFreezeThaw, HessianCheck) {
  pdf_1d::LimitedExponentialPDF exp_pdf(1.0, 9.0);
  pdf_1d::LimitedGaussianPDF gauss_pdf(1.0, 9.0);
  pdf_1d::TwoComponent1DPDF pdf_x(&exp_pdf, "exp", &gauss_pdf, "gauss");
  Eigen::VectorXd p(3);
  p << 0.5, 100.0, 35.0;
  Eigen::VectorXd dp(3);
  dp << 1e-6, 1e-6, 1e-6;
  SingleToParameterizableMultiAxisFunctionAdapter maf(&pdf_x);
  function::PMAFReverser pdf(&maf);
  function::FreezeThawFunction freezer(&pdf);
  freezer.freeze(1,40.0);

  Eigen::VectorXd x(1);
  for(x(0) = 0; x(0)<10.0; x(0)+=0.1)
  {
    pdf.set_parameter_values(x);
    EXPECT_EQ(x, pdf.parameter_values());
    Eigen::MatrixXd good(3,3);
    EXPECT_TRUE(hessian_check(freezer, p, dp, good));
    EXPECT_LE(good(0,0), 0.5);
    EXPECT_LE(good(0,1), 0.5);
    EXPECT_LE(good(0,2), 0.5);
    EXPECT_LE(good(1,0), 0.5);
    EXPECT_LE(good(1,1), 0.5);
    EXPECT_LE(good(1,2), 0.5);
    EXPECT_LE(good(2,0), 0.5);
    EXPECT_LE(good(2,1), 0.5);
    EXPECT_LE(good(2,2), 0.5);
  }
}

TEST(TestFreezeThaw, ParameterGradientCheck) {
  pdf_1d::LimitedExponentialPDF exp_pdf(1.0, 9.0);
  pdf_1d::LimitedGaussianPDF gauss_pdf(1.0, 9.0);
  pdf_1d::TwoComponent1DPDF pdf_x(&exp_pdf, "exp", &gauss_pdf, "gauss");
  Eigen::VectorXd p(2);
  p << 40, 100.0;
  Eigen::VectorXd dp(2);
  dp << 1e-6, 1e-6;
  SingleToParameterizableMultiAxisFunctionAdapter maf(&pdf_x);
  function::PMAFReverser pdf(&maf);
  function::FreezeThawFunction freezer(&pdf);
  freezer.freeze(1,40.0);
  freezer.freeze(2,100.0);
  function::PMAFReverser rev_freezer(&freezer);

  Eigen::VectorXd p2(2);
  p2 << 0.5, 35.0;
  rev_freezer.set_parameter_values(p2);

  Eigen::VectorXd x(1);
  for(x(0) = 0; x(0)<10.0; x(0)+=0.1)
  {
    pdf.set_parameter_values(x);
    EXPECT_EQ(x, pdf.parameter_values());
    Eigen::VectorXd good(2);
    EXPECT_TRUE(gradient_check(rev_freezer, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
  }
}

TEST(TestFreezeThaw, ParameterHessianCheck) {
  pdf_1d::LimitedExponentialPDF exp_pdf(1.0, 9.0);
  pdf_1d::LimitedGaussianPDF gauss_pdf(1.0, 9.0);
  pdf_1d::TwoComponent1DPDF pdf_x(&exp_pdf, "exp", &gauss_pdf, "gauss");
  Eigen::VectorXd p(2);
  p << 40, 100.0;
  Eigen::VectorXd dp(2);
  dp << 1e-6, 1e-6;
  SingleToParameterizableMultiAxisFunctionAdapter maf(&pdf_x);
  function::PMAFReverser pdf(&maf);
  function::FreezeThawFunction freezer(&pdf);
  freezer.freeze(1,40.0);
  freezer.freeze(2,100.0);
  function::PMAFReverser rev_freezer(&freezer);

  Eigen::VectorXd p2(2);
  p2 << 0.5, 35.0;
  rev_freezer.set_parameter_values(p2);

  Eigen::VectorXd x(1);
  for(x(0) = 0; x(0)<10.0; x(0)+=0.1)
  {
    pdf.set_parameter_values(x);
    EXPECT_EQ(x, pdf.parameter_values());
    Eigen::MatrixXd good(2,2);
    EXPECT_TRUE(hessian_check(rev_freezer, p, dp, good));
    EXPECT_LE(good(0,0), 0.5);
    EXPECT_LE(good(0,1), 0.5);
    EXPECT_LE(good(1,0), 0.5);
    EXPECT_LE(good(1,1), 0.5);
  }
}

TEST(TestLogQuadSpline, UnnormalisedValues) {
  Eigen::VectorXd xknot(3);
  xknot << -1.0, 0.0, 1.0;
  pdf_1d::LogQuadraticSpline1DPDF pdf(xknot, -1.5, 1.5, 0.0, ParamZeroType::SLOPE, ParamZeroLocation::LEFT, false);
  Eigen::VectorXd p(4);
  p << 0, 1.0, 2.0, 1.0;
  pdf.set_parameter_values(p);
  Eigen::VectorXd p_readback;
  p_readback = pdf.parameter_values();
  EXPECT_EQ(p, p_readback);
  for(unsigned i=0;i<xknot.size();i++)
    EXPECT_EQ(pdf.value_1d(xknot(i)), std::exp(p(i+1)));
}

TEST(TestLogQuadSpline, UnnormalisedGradientCheck) {
  Eigen::VectorXd xknot(3);
  xknot << -1.0, 0.0, 1.0;
  pdf_1d::LogQuadraticSpline1DPDF pdf(xknot, -1.5, 1.5, 0.0, ParamZeroType::SLOPE, ParamZeroLocation::LEFT, false);
  Eigen::VectorXd p(4);
  p << 0, 1.0, 2.0, 1.0;
  pdf.set_parameter_values(p);
  for(double x=-1.4; x<1.4; x+=0.1)
  {
    double good;
    EXPECT_TRUE(gradient_check(pdf, x, 0.00001, good));
    EXPECT_LE(good, 0.5);
  }
}

TEST(TestLogQuadSpline, UnnormalisedParameterGradientCheck) {
  Eigen::VectorXd xknot(3);
  xknot << -1.0, 0.0, 1.0;
  pdf_1d::LogQuadraticSpline1DPDF pdf(xknot, -1.5, 1.5, 0.0, ParamZeroType::SLOPE, ParamZeroLocation::LEFT, false);

  Eigen::VectorXd p(4);
  p << -0.5, 1.0, 2.0, 1.5;
  Eigen::VectorXd dp(4);
  dp << 0.001, 0.001, 0.001, 0.001;

  for(double x=-1.4; x<1.4; x+=0.1)
  {
    Eigen::VectorXd good(4);
    EXPECT_TRUE(gradient_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
    EXPECT_LE(good(2), 0.5);
    EXPECT_LE(good(3), 0.5);
  }
}

TEST(TestLogQuadSpline, NormalisedValues_ANeg) {
  Eigen::VectorXd xknot(3);
  xknot << -1.0, 0.0, 1.0;
  pdf_1d::LogQuadraticSpline1DPDF pdf(xknot, -1.5, 1.5);
  Eigen::VectorXd p(4);
  p << 2.0, 1.0, 2.0, 1.0;
  pdf.set_parameter_values(p);
  Eigen::VectorXd p_readback;
  p_readback = pdf.parameter_values();
  EXPECT_EQ(p, p_readback);
  ASSERT_EQ(pdf.a().size(), 2);
  ASSERT_EQ(pdf.b().size(), 2);
  ASSERT_LT(pdf.a()(0), 0);
  ASSERT_LT(pdf.a()(1), 0);
  double dx = 0.001;
  double sum = 0;
  for(double x=-1.5+dx/2; x<1.5;x+=dx)sum += pdf.value_1d(x);
  EXPECT_NEAR(sum*dx, 1.0, dx);
}

TEST(TestLogQuadSpline, NormalisedValues_APos) {
  Eigen::VectorXd xknot(3);
  xknot << -1.0, 0.0, 1.0;
  pdf_1d::LogQuadraticSpline1DPDF pdf(xknot, -1.5, 1.5);
  Eigen::VectorXd p(4);
  p << -1.5, 2.0, 1.0, 2.0;
  pdf.set_parameter_values(p);
  Eigen::VectorXd p_readback;
  p_readback = pdf.parameter_values();
  EXPECT_EQ(p, p_readback);
  ASSERT_EQ(pdf.a().size(), 2);
  ASSERT_EQ(pdf.b().size(), 2);
  ASSERT_GT(pdf.a()(0), 0);
  ASSERT_GT(pdf.a()(1), 0);
  double dx = 0.001;
  double sum = 0;
  for(double x=-1.5+dx/2; x<1.5;x+=dx)sum += pdf.value_1d(x);
  EXPECT_NEAR(sum*dx, 1.0, dx);
}

TEST(TestLogQuadSpline, NormalisedValues_AZero) {
  Eigen::VectorXd xknot(3);
  xknot << -1.0, 0.0, 1.0;
  pdf_1d::LogQuadraticSpline1DPDF pdf(xknot, -1.5, 1.5);
  Eigen::VectorXd p(4);
  p << 1.0, 1.0, 2.0, 3.0;
  pdf.set_parameter_values(p);
  Eigen::VectorXd p_readback;
  p_readback = pdf.parameter_values();
  EXPECT_EQ(p, p_readback);
  ASSERT_EQ(pdf.a().size(), 2);
  ASSERT_EQ(pdf.b().size(), 2);
  ASSERT_EQ(pdf.a()(0), 0);
  ASSERT_EQ(pdf.a()(1), 0);
  double dx = 0.001;
  double sum = 0;
  for(double x=-1.5+dx/2; x<1.5;x+=dx)sum += pdf.value_1d(x);
  EXPECT_NEAR(sum*dx, 1.0, dx);
}

TEST(TestLogQuadSpline, NormalisedValues_ABZero) {
  Eigen::VectorXd xknot(3);
  xknot << -1.0, 0.0, 1.0;
  pdf_1d::LogQuadraticSpline1DPDF pdf(xknot, -1.5, 1.5);
  Eigen::VectorXd p(4);
  p << 0.0, 1.0, 1.0, 1.0;
  pdf.set_parameter_values(p);
  Eigen::VectorXd p_readback;
  p_readback = pdf.parameter_values();
  EXPECT_EQ(p, p_readback);
  ASSERT_EQ(pdf.a().size(), 2);
  ASSERT_EQ(pdf.b().size(), 2);
  ASSERT_EQ(pdf.a()(0), 0);
  ASSERT_EQ(pdf.a()(1), 0);
  ASSERT_EQ(pdf.b()(0), 0);
  ASSERT_EQ(pdf.b()(1), 0);
  double dx = 0.001;
  double sum = 0;
  for(double x=-1.5+dx/2; x<1.5;x+=dx)sum += pdf.value_1d(x);
  EXPECT_NEAR(sum*dx, 1.0, dx);
}

TEST(TestLogQuadSpline, NormalisedParameterGradientCheck_ANeg) {
  Eigen::VectorXd xknot(3);
  xknot << -1.0, 0.0, 1.0;
  pdf_1d::LogQuadraticSpline1DPDF pdf(xknot, -1.5, 1.5);

  Eigen::VectorXd p(4);
  p << 2.0, 1.0, 2.0, 1.0;
  Eigen::VectorXd dp(4);
  dp << 0.001, 0.001, 0.001, 0.001;

  for(double x=-1.4; x<1.4; x+=0.1)
  {
    Eigen::VectorXd good(4);
    EXPECT_TRUE(gradient_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
    EXPECT_LE(good(2), 0.5);
    EXPECT_LE(good(3), 0.5);
  }
}

TEST(TestLogQuadSpline, NormalisedParameterGradientCheck_APos) {
  Eigen::VectorXd xknot(3);
  xknot << -1.0, 0.0, 1.0;
  pdf_1d::LogQuadraticSpline1DPDF pdf(xknot, -1.5, 1.5);

  Eigen::VectorXd p(4);
  p << -1.5, 2.0, 1.0, 2.0;
  Eigen::VectorXd dp(4);
  dp << 0.001, 0.001, 0.001, 0.001;

  for(double x=-1.4; x<1.4; x+=0.1)
  {
    Eigen::VectorXd good(4);
    EXPECT_TRUE(gradient_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
    EXPECT_LE(good(2), 0.5);
    EXPECT_LE(good(3), 0.5);
  }
}

TEST(TestLogQuadSpline, NormalisedParameterGradientCheck_AZero) {
  Eigen::VectorXd xknot(3);
  xknot << -1.0, 0.0, 1.0;
  pdf_1d::LogQuadraticSpline1DPDF pdf(xknot, -1.5, 1.5);

  Eigen::VectorXd p(4);
  p << 1.0, 1.0, 2.0, 3.0;
  Eigen::VectorXd dp(4);
  dp << 0.001, 0.001, 0.001, 0.001;

  for(double x=-1.4; x<1.4; x+=0.1)
  {
    Eigen::VectorXd good(4);
    EXPECT_TRUE(gradient_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
    EXPECT_LE(good(2), 0.5);
    EXPECT_LE(good(3), 0.5);
  }
}

TEST(TestLogQuadSpline, NormalisedParameterGradientCheck_ABZero) {
  Eigen::VectorXd xknot(3);
  xknot << -1.0, 0.0, 1.0;
  pdf_1d::LogQuadraticSpline1DPDF pdf(xknot, -1.5, 1.5);

  Eigen::VectorXd p(4);
  p << 0.0, 1.0, 1.0, 1.0;
  Eigen::VectorXd dp(4);
  dp << 0.001, 0.001, 0.001, 0.001;

  for(double x=-1.4; x<1.4; x+=0.1)
  {
    Eigen::VectorXd good(4);
    EXPECT_TRUE(gradient_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
    EXPECT_LE(good(2), 0.5);
    EXPECT_LE(good(3), 0.5);
  }
}

TEST(TestLogQuadSpline, BinnedValues_ANeg) {
  double dx = 0.001;
  Eigen::VectorXd xknot(3);
  xknot << -1.0, 0.0, 1.0;
  pdf_1d::LogQuadraticSpline1DPDF pdf(xknot, -1.5, 1.5, dx);
  Eigen::VectorXd p(4);
  p << 2.0, 1.0, 2.0, 1.0;
  pdf.set_parameter_values(p);
  Eigen::VectorXd p_readback;
  p_readback = pdf.parameter_values();
  EXPECT_EQ(p, p_readback);
  ASSERT_EQ(pdf.a().size(), 2);
  ASSERT_EQ(pdf.b().size(), 2);
  ASSERT_LT(pdf.a()(0), 0);
  ASSERT_LT(pdf.a()(1), 0);
  double sum = 0;
  for(double x=-1.5+dx/2; x<1.5;x+=dx)sum += pdf.value_1d(x);
  EXPECT_NEAR(sum*dx, 1.0, dx);
}

TEST(TestLogQuadSpline, BinnedValues_APos) {
  double dx = 0.001;
  Eigen::VectorXd xknot(3);
  xknot << -1.0, 0.0, 1.0;
  pdf_1d::LogQuadraticSpline1DPDF pdf(xknot, -1.5, 1.5, dx);
  Eigen::VectorXd p(4);
  p << -1.5, 2.0, 1.0, 2.0;
  pdf.set_parameter_values(p);
  Eigen::VectorXd p_readback;
  p_readback = pdf.parameter_values();
  EXPECT_EQ(p, p_readback);
  ASSERT_EQ(pdf.a().size(), 2);
  ASSERT_EQ(pdf.b().size(), 2);
  ASSERT_GT(pdf.a()(0), 0);
  ASSERT_GT(pdf.a()(1), 0);
  double sum = 0;
  for(double x=-1.5+dx/2; x<1.5;x+=dx)sum += pdf.value_1d(x);
  EXPECT_NEAR(sum*dx, 1.0, dx);
}

TEST(TestLogQuadSpline, BinnedValues_AZero) {
  double dx = 0.001;
  Eigen::VectorXd xknot(3);
  xknot << -1.0, 0.0, 1.0;
  pdf_1d::LogQuadraticSpline1DPDF pdf(xknot, -1.5, 1.5, dx);
  Eigen::VectorXd p(4);
  p << 1.0, 1.0, 2.0, 3.0;
  pdf.set_parameter_values(p);
  Eigen::VectorXd p_readback;
  p_readback = pdf.parameter_values();
  EXPECT_EQ(p, p_readback);
  ASSERT_EQ(pdf.a().size(), 2);
  ASSERT_EQ(pdf.b().size(), 2);
  ASSERT_EQ(pdf.a()(0), 0);
  ASSERT_EQ(pdf.a()(1), 0);
  double sum = 0;
  for(double x=-1.5+dx/2; x<1.5;x+=dx)sum += pdf.value_1d(x);
  EXPECT_NEAR(sum*dx, 1.0, dx);
}

TEST(TestLogQuadSpline, BinnedValues_ABZero) {
  double dx = 0.001;
  Eigen::VectorXd xknot(3);
  xknot << -1.0, 0.0, 1.0;
  pdf_1d::LogQuadraticSpline1DPDF pdf(xknot, -1.5, 1.5, dx);
  Eigen::VectorXd p(4);
  p << 0.0, 1.0, 1.0, 1.0;
  pdf.set_parameter_values(p);
  Eigen::VectorXd p_readback;
  p_readback = pdf.parameter_values();
  EXPECT_EQ(p, p_readback);
  ASSERT_EQ(pdf.a().size(), 2);
  ASSERT_EQ(pdf.b().size(), 2);
  ASSERT_EQ(pdf.a()(0), 0);
  ASSERT_EQ(pdf.a()(1), 0);
  ASSERT_EQ(pdf.b()(0), 0);
  ASSERT_EQ(pdf.b()(1), 0);
  double sum = 0;
  for(double x=-1.5+dx/2; x<1.5;x+=dx)sum += pdf.value_1d(x);
  EXPECT_NEAR(sum*dx, 1.0, dx);
}

TEST(TestLogQuadSpline, BinnedParameterGradientCheck_ANeg) {
  Eigen::VectorXd xknot(3);
  xknot << -1.0, 0.0, 1.0;
  pdf_1d::LogQuadraticSpline1DPDF pdf(xknot, -1.5, 1.5, 0.1);

  Eigen::VectorXd p(4);
  p << 2.0, 1.0, 2.0, 1.0;
  Eigen::VectorXd dp(4);
  dp << 0.001, 0.001, 0.001, 0.001;

  for(double x=-1.5; x<1.5; x+=0.1)
  {
    Eigen::VectorXd good(4);
    EXPECT_TRUE(gradient_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
    EXPECT_LE(good(2), 0.5);
    EXPECT_LE(good(3), 0.5);
  }
}

TEST(TestLogQuadSpline, BinnedParameterGradientCheck_APos) {
  Eigen::VectorXd xknot(3);
  xknot << -1.0, 0.0, 1.0;
  pdf_1d::LogQuadraticSpline1DPDF pdf(xknot, -1.5, 1.5, 0.1);

  Eigen::VectorXd p(4);
  p << -1.5, 2.0, 1.0, 2.0;
  Eigen::VectorXd dp(4);
  dp << 0.001, 0.001, 0.001, 0.001;

  for(double x=-1.5; x<1.5; x+=0.1)
  {
    Eigen::VectorXd good(4);
    EXPECT_TRUE(gradient_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
    EXPECT_LE(good(2), 0.5);
    EXPECT_LE(good(3), 0.5);
  }
}

TEST(TestLogQuadSpline, BinnedParameterGradientCheck_AZero) {
  Eigen::VectorXd xknot(3);
  xknot << -1.0, 0.0, 1.0;
  pdf_1d::LogQuadraticSpline1DPDF pdf(xknot, -1.5, 1.5, 0.1);

  Eigen::VectorXd p(4);
  p << 1.0, 1.0, 2.0, 3.0;
  Eigen::VectorXd dp(4);
  dp << 0.001, 0.001, 0.001, 0.001;

  for(double x=-1.5; x<1.5; x+=0.1)
  {
    Eigen::VectorXd good(4);
    EXPECT_TRUE(gradient_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
    EXPECT_LE(good(2), 0.5);
    EXPECT_LE(good(3), 0.5);
  }
}

TEST(TestLogQuadSpline, BinnedParameterGradientCheck_ABZero) {
  Eigen::VectorXd xknot(3);
  xknot << -1.0, 0.0, 1.0;
  pdf_1d::LogQuadraticSpline1DPDF pdf(xknot, -1.5, 1.5, 0.1);

  Eigen::VectorXd p(4);
  p << 0.0, 1.0, 1.0, 1.0;
  Eigen::VectorXd dp(4);
  dp << 0.001, 0.001, 0.001, 0.001;

  for(double x=-1.5; x<1.5; x+=0.1)
  {
    Eigen::VectorXd good(4);
    EXPECT_TRUE(gradient_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
    EXPECT_LE(good(2), 0.5);
    EXPECT_LE(good(3), 0.5);
  }
}


TEST(TestTwoComponent1DConstraintPDF, SetAndRecallParameters) {
  pdf_1d::LimitedGaussianPDF gauss1_pdf(1.0, 9.0);
  pdf_1d::LimitedGaussianPDF gauss2_pdf(1.0, 9.0);
  pdf_1d::TwoComponent1DPDF TwoComponent(&gauss1_pdf, "gauss_low", &gauss2_pdf, "gauss_high");
  pdf_1d::TwoComponent1DConstraintPDF pdf(&TwoComponent);

  EXPECT_EQ(pdf.num_parameters(), 5U);
  auto pvec = pdf.parameters();
  EXPECT_EQ(pvec[0].name, "probability");
  EXPECT_EQ(pvec[1].name, "mean_low_gauss");
  EXPECT_EQ(pvec[2].name, "resolution");
  EXPECT_EQ(pvec[3].name, "mean_high_gauss");
  EXPECT_EQ(pvec[4].name, "Width_Ratio");
  Eigen::VectorXd p(5);
  p << 0.4, 0, 0.4, 5.0, 0.5;
  pdf.set_parameter_values(p);
  EXPECT_EQ(pdf.parameter_values(), p);
}

TEST(TestTwoComponent1DConstraintPDF, GradientCheck) {
  pdf_1d::LimitedGaussianPDF gauss1_pdf(1.0, 9.0);
  pdf_1d::LimitedGaussianPDF gauss2_pdf(1.0, 9.0);
  pdf_1d::TwoComponent1DPDF TwoComponent(&gauss1_pdf, "gauss_low", &gauss2_pdf, "gauss_high");
  pdf_1d::TwoComponent1DConstraintPDF pdf(&TwoComponent);
  Eigen::VectorXd p(5);
  p << 0.4, 0, 0.4, 5.0, 0.5;
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(double x = 0; x<10.0; x+=0.1)
  {
    double good;
    EXPECT_TRUE(gradient_check(pdf, x, 1e-6, good));
    EXPECT_LE(good, 0.5);
  }
}

TEST(TestTwoComponent1DConstraintPDF, HessianCheck) {
  pdf_1d::LimitedGaussianPDF gauss1_pdf(1.0, 9.0);
  pdf_1d::LimitedGaussianPDF gauss2_pdf(1.0, 9.0);
  pdf_1d::TwoComponent1DPDF TwoComponent(&gauss1_pdf, "gauss_low", &gauss2_pdf, "gauss_high");
  pdf_1d::TwoComponent1DConstraintPDF pdf(&TwoComponent);
  Eigen::VectorXd p(5);
  p << 0.4, 0, 0.4, 5.0, 0.5;
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(double x = 0; x<10.0; x+=0.1)
  {
    double good;
    EXPECT_TRUE(hessian_check(pdf, x, 1e-6, good));
    EXPECT_LE(good, 0.5);
  }
}

TEST(TestTwoComponent1DConstraintPDF, ParameterGradientCheck) {
  pdf_1d::LimitedGaussianPDF gauss1_pdf(1.0, 9.0);
  pdf_1d::LimitedGaussianPDF gauss2_pdf(1.0, 9.0);
  pdf_1d::TwoComponent1DPDF TwoComponent(&gauss1_pdf, "gauss_low", &gauss2_pdf, "gauss_high");
  pdf_1d::TwoComponent1DConstraintPDF pdf(&TwoComponent);
  Eigen::VectorXd p(5);
  p << 0.4, 0, 0.4, 5.0, 0.5;
  Eigen::VectorXd dp(5);
  dp << 1e-6, 1e-6, 1e-6, 1e-6, 1e-6;
  for(double x = 0; x<10.0; x+=0.1)
  {
    Eigen::VectorXd good(5);
    EXPECT_TRUE(gradient_check_par(pdf, x, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
    EXPECT_LE(good(2), 0.5);
    EXPECT_LE(good(3), 0.5);
    EXPECT_LE(good(4), 0.5);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
