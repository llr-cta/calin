#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>
#include <vector>

#include "Eigen/Dense"
#include "math/function.hpp"

using namespace calin::math;

TEST(TestGaussianPDF, GradientCheck) {
  Eigen::VectorXd p(2);
  p << 3, 2; // mean, sigma
  Eigen::VectorXd x(1);
  Eigen::VectorXd dx(1);
  dx << 1e-6;
  GaussianPDF pdf;
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(x(0) = 0; x(0)<10.0; x(0)+=0.1)
  {
    Eigen::VectorXd good(1);
    EXPECT_TRUE(gradient_check(pdf, x, dx, good));
    EXPECT_LE(good(0), 0.5);
  }
}

TEST(TestGaussianPDF, HessianCheck) {
  Eigen::VectorXd p(2);
  p << 3, 2; // mean, sigma
  Eigen::VectorXd x(1);
  Eigen::VectorXd dx(1);
  dx << 1e-6;
  GaussianPDF pdf;
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(x(0) = 0; x(0)<10.0; x(0)+=0.1)
  {
    Eigen::MatrixXd good(1,1);
    EXPECT_TRUE(hessian_check(pdf, x, dx, good));
    EXPECT_LE(good(0,0), 0.5);
  }
}

TEST(TestGaussianPDF, ParameterGradientCheck) {
  Eigen::VectorXd p(2);
  p << 3, 2; // mean, sigma
  Eigen::VectorXd dp(2);
  dp << 1e-6, 1e-6;
  GaussianPDF pdf_x;
  function::PMAFReverser pdf(&pdf_x);
  Eigen::VectorXd x(1);
  for(x(0) = 0; x(0)<10.0; x(0)+=0.1)
  {
    pdf.set_parameter_values(x);
    EXPECT_EQ(x, pdf.parameter_values());
    Eigen::VectorXd good(2);
    EXPECT_TRUE(gradient_check(pdf, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
  }
}

TEST(TestGaussianPDF, ParameterHessianCheck) {
  Eigen::VectorXd p(2);
  p << 3, 2; // mean, sigma
  Eigen::VectorXd dp(2);
  dp << 1e-6, 1e-6;
  GaussianPDF pdf_x;
  function::PMAFReverser pdf(&pdf_x);
  Eigen::VectorXd x(1);
  for(x(0) = 0; x(0)<10.0; x(0)+=0.1)
  {
    pdf.set_parameter_values(x);
    EXPECT_EQ(x, pdf.parameter_values());
    Eigen::MatrixXd good(2,2);
    EXPECT_TRUE(hessian_check(pdf, p, dp, good));
    EXPECT_LE(good(0,0), 0.5);
    EXPECT_LE(good(0,1), 0.5);
    EXPECT_LE(good(1,0), 0.5);
    EXPECT_LE(good(1,1), 0.5);
  }
}

TEST(TestLimitedGaussianPDF, GradientCheck) {
  Eigen::VectorXd p(2);
  p << 3, 2; // mean, sigma
  Eigen::VectorXd x(1);
  Eigen::VectorXd dx(1);
  dx << 1e-6;
  LimitedGaussianPDF pdf(1.0,9.0);
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(x(0) = 0; x(0)<10.0; x(0)+=0.1)
  {
    Eigen::VectorXd good(1);
    EXPECT_TRUE(gradient_check(pdf, x, dx, good));
    EXPECT_LE(good(0), 0.5);
  }
}

TEST(TestLimitedGaussianPDF, HessianCheck) {
  Eigen::VectorXd p(2);
  p << 3, 2; // mean, sigma
  Eigen::VectorXd x(1);
  Eigen::VectorXd dx(1);
  dx << 1e-6;
  LimitedGaussianPDF pdf(1.0,9.0);
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(x(0) = 0; x(0)<10.0; x(0)+=0.1)
  {
    Eigen::MatrixXd good(1,1);
    EXPECT_TRUE(hessian_check(pdf, x, dx, good));
    EXPECT_LE(good(0,0), 0.5);
  }
}

TEST(TestLimitedGaussianPDF, ParameterGradientCheck) {
  Eigen::VectorXd p(2);
  p << 3, 2; // mean, sigma
  Eigen::VectorXd dp(2);
  dp << 1e-6, 1e-6;
  LimitedGaussianPDF pdf_x(1.0,9.0);
  function::PMAFReverser pdf(&pdf_x);
  Eigen::VectorXd x(1);
  for(x(0) = 0; x(0)<10.0; x(0)+=0.1)
  {
    pdf.set_parameter_values(x);
    EXPECT_EQ(x, pdf.parameter_values());
    Eigen::VectorXd good(2);
    EXPECT_TRUE(gradient_check(pdf, p, dp, good));
    EXPECT_LE(good(0), 0.5);
    EXPECT_LE(good(1), 0.5);
  }
}

TEST(TestLimitedGaussianPDF, ParameterHessianCheck) {
  Eigen::VectorXd p(2);
  p << 3, 2; // mean, sigma
  Eigen::VectorXd dp(2);
  dp << 1e-6, 1e-6;
  LimitedGaussianPDF pdf_x(1.0,9.0);
  function::PMAFReverser pdf(&pdf_x);
  Eigen::VectorXd x(1);
  for(x(0) = 0; x(0)<10.0; x(0)+=0.1)
  {
    pdf.set_parameter_values(x);
    EXPECT_EQ(x, pdf.parameter_values());
    Eigen::MatrixXd good(2,2);
    EXPECT_TRUE(hessian_check(pdf, p, dp, good));
    EXPECT_LE(good(0,0), 0.5);
    EXPECT_LE(good(0,1), 0.5);
    EXPECT_LE(good(1,0), 0.5);
    EXPECT_LE(good(1,1), 0.5);
  }
}

TEST(TestLimitedExponentialPDF, GradientCheck) {
  Eigen::VectorXd p(1);
  p << 3; // scale
  Eigen::VectorXd x(1);
  Eigen::VectorXd dx(1);
  dx << 1e-6;
  LimitedExponentialPDF pdf(1.0,9.0);
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(x(0) = 0; x(0)<10.0; x(0)+=0.1)
  {
    Eigen::VectorXd good(1);
    EXPECT_TRUE(gradient_check(pdf, x, dx, good));
    EXPECT_LE(good(0), 0.5);
  }
}

TEST(TestLimitedExponentialPDF, HessianCheck) {
  Eigen::VectorXd p(1);
  p << 3; // scale
  Eigen::VectorXd x(1);
  Eigen::VectorXd dx(1);
  dx << 1e-6;
  LimitedExponentialPDF pdf(1.0,9.0);
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(x(0) = 0; x(0)<10.0; x(0)+=0.1)
  {
    Eigen::MatrixXd good(1,1);
    EXPECT_TRUE(hessian_check(pdf, x, dx, good));
    EXPECT_LE(good(0,0), 0.5);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
