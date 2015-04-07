#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>
#include <vector>

#include "Eigen/Dense"
#include "math/pdf_1d.hpp"

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

TEST(TestLimitedExponentialPDF, ParameterGradientCheck) {
  Eigen::VectorXd p(1);
  p << 3; // scale
  Eigen::VectorXd dp(1);
  dp << 1e-6;
  LimitedExponentialPDF pdf_x(1.0,9.0);
  function::PMAFReverser pdf(&pdf_x);
  Eigen::VectorXd x(1);
  for(x(0) = 0; x(0)<10.0; x(0)+=0.1)
  {
    pdf.set_parameter_values(x);
    EXPECT_EQ(x, pdf.parameter_values());
    Eigen::VectorXd good(1);
    EXPECT_TRUE(gradient_check(pdf, p, dp, good));
    EXPECT_LE(good(0), 0.5);
  }
}

TEST(TestLimitedExponentialPDF, ParameterHessianCheck) {
  Eigen::VectorXd p(1);
  p << 3; // scale
  Eigen::VectorXd dp(1);
  dp << 1e-6;
  LimitedExponentialPDF pdf_x(1.0,9.0);
  function::PMAFReverser pdf(&pdf_x);
  Eigen::VectorXd x(1);
  for(x(0) = 0; x(0)<10.0; x(0)+=0.1)
  {
    pdf.set_parameter_values(x);
    EXPECT_EQ(x, pdf.parameter_values());
    Eigen::MatrixXd good(1,1);
    EXPECT_TRUE(hessian_check(pdf, p, dp, good));
    EXPECT_LE(good(0,0), 0.5);
  }
}

TEST(TestTwoComponentPDF_ExpGauss, SetAndRecallParameters) {
  pdf_1d::LimitedExponentialPDF exp_pdf(1.0, 9.0);
  pdf_1d::LimitedGaussianPDF gauss_pdf(1.0, 9.0);
  pdf_1d::TwoComponentPDF pdf(&exp_pdf, "exp", &gauss_pdf, "gauss");

  EXPECT_EQ(pdf.num_parameters(), 4);
  auto pvec = pdf.parameters();
  EXPECT_EQ(pvec[0].name, "exp probability");
  EXPECT_EQ(pvec[1].name, "exp.scale");
  EXPECT_EQ(pvec[2].name, "gauss.mean");
  EXPECT_EQ(pvec[3].name, "gauss.rms");
  Eigen::VectorXd p(4);
  p << 0.5, 40.0, 100.0, 35.0;
  pdf.set_parameter_values(p);
  EXPECT_EQ(pdf.parameter_values(), p);
}

TEST(TestTwoComponentPDF_ExpGauss, GradientCheck) {
  pdf_1d::LimitedExponentialPDF exp_pdf(1.0, 9.0);
  pdf_1d::LimitedGaussianPDF gauss_pdf(1.0, 9.0);
  pdf_1d::TwoComponentPDF pdf(&exp_pdf, "exp", &gauss_pdf, "gauss");
  Eigen::VectorXd p(4);
  p << 0.5, 40.0, 100.0, 35.0;
  Eigen::VectorXd x(1);
  Eigen::VectorXd dx(1);
  dx << 1e-6;
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(x(0) = 0; x(0)<10.0; x(0)+=0.1)
  {
    Eigen::VectorXd good(1);
    EXPECT_TRUE(gradient_check(pdf, x, dx, good));
    EXPECT_LE(good(0), 0.5);
  }
}

TEST(TestTwoComponentPDF_ExpGauss, HessianCheck) {
  pdf_1d::LimitedExponentialPDF exp_pdf(1.0, 9.0);
  pdf_1d::LimitedGaussianPDF gauss_pdf(1.0, 9.0);
  pdf_1d::TwoComponentPDF pdf(&exp_pdf, "exp", &gauss_pdf, "gauss");
  Eigen::VectorXd p(4);
  p << 0.5, 40.0, 100.0, 35.0;
  Eigen::VectorXd x(1);
  Eigen::VectorXd dx(1);
  dx << 1e-6;
  pdf.set_parameter_values(p);
  EXPECT_EQ(p, pdf.parameter_values());
  for(x(0) = 0; x(0)<10.0; x(0)+=0.1)
  {
    Eigen::MatrixXd good(1,1);
    EXPECT_TRUE(hessian_check(pdf, x, dx, good));
    EXPECT_LE(good(0,0), 0.5);
  }
}

TEST(TestTwoComponentPDF_ExpGauss, ParameterGradientCheck) {
  pdf_1d::LimitedExponentialPDF exp_pdf(1.0, 9.0);
  pdf_1d::LimitedGaussianPDF gauss_pdf(1.0, 9.0);
  pdf_1d::TwoComponentPDF pdf_x(&exp_pdf, "exp", &gauss_pdf, "gauss");
  Eigen::VectorXd p(4);
  p << 0.5, 40.0, 100.0, 35.0;
  Eigen::VectorXd dp(4);
  dp << 1e-6, 1e-6, 1e-6, 1e-6;
  function::PMAFReverser pdf(&pdf_x);
  Eigen::VectorXd x(1);
  for(x(0) = 0; x(0)<10.0; x(0)+=0.1)
  {
    pdf.set_parameter_values(x);
    EXPECT_EQ(x, pdf.parameter_values());
    Eigen::VectorXd good(4);
    EXPECT_TRUE(gradient_check(pdf, p, dp, good));
    EXPECT_LE(good(0), 0.5);
  }
}

TEST(TestTwoComponentPDF_ExpGauss, ParameterHessianCheck) {
  pdf_1d::LimitedExponentialPDF exp_pdf(1.0, 9.0);
  pdf_1d::LimitedGaussianPDF gauss_pdf(1.0, 9.0);
  pdf_1d::TwoComponentPDF pdf_x(&exp_pdf, "exp", &gauss_pdf, "gauss");
  Eigen::VectorXd p(4);
  p << 0.5, 40.0, 100.0, 35.0;
  Eigen::VectorXd dp(4);
  dp << 1e-6, 1e-6, 1e-6, 1e-6;
  function::PMAFReverser pdf(&pdf_x);
  Eigen::VectorXd x(1);
  for(x(0) = 0; x(0)<10.0; x(0)+=0.1)
  {
    pdf.set_parameter_values(x);
    EXPECT_EQ(x, pdf.parameter_values());
    Eigen::MatrixXd good(4,4);
    EXPECT_TRUE(hessian_check(pdf, p, dp, good));
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
