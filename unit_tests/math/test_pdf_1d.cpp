#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>
#include <vector>

#include "Eigen/Dense"
#include "math/pdf_1d.hpp"

using namespace calin::math;
using namespace calin::math::function;
using namespace calin::math::pdf_1d;

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

TEST(TestLimitedExponentialPDF, Binned_ParameterGradientCheck) {
  Eigen::VectorXd p(1);
  p << 3; // scale
  Eigen::VectorXd dp(1);
  dp << 1e-6;
  LimitedExponentialPDF pdf_x(1.0,9.0,0.1);
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

TEST(TestLimitedExponentialPDF, Binned_ParameterHessianCheck) {
  Eigen::VectorXd p(1);
  p << 3; // scale
  Eigen::VectorXd dp(1);
  dp << 1e-6;
  LimitedExponentialPDF pdf_x(1.0,9.0,0.1);
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
    EXPECT_LE(good(1), 0.5);
    EXPECT_LE(good(2), 0.5);
    EXPECT_LE(good(3), 0.5);
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

TEST(TestFreezeThaw, FreezeAndThaw) {
  pdf_1d::LimitedExponentialPDF exp_pdf(1.0, 9.0);
  pdf_1d::LimitedGaussianPDF gauss_pdf(1.0, 9.0);
  pdf_1d::TwoComponentPDF pdf_x(&exp_pdf, "exp", &gauss_pdf, "gauss");
  function::PMAFReverser pdf(&pdf_x);
  function::FreezeThawFunction freezer(&pdf);
  EXPECT_EQ(freezer.num_domain_axes(), 4);
  EXPECT_EQ(freezer.num_parameters(), 0);
  EXPECT_EQ(freezer.num_domain_axes(), freezer.domain_axes().size());
  EXPECT_EQ(freezer.num_parameters(), freezer.parameters().size());
  
  EXPECT_TRUE(freezer.freeze(1,50.0));
  EXPECT_EQ(freezer.num_domain_axes(), 3);
  EXPECT_EQ(freezer.num_parameters(), 1);
  EXPECT_EQ(freezer.num_domain_axes(), freezer.domain_axes().size());
  EXPECT_EQ(freezer.num_parameters(), freezer.parameters().size());
  EXPECT_EQ(freezer.free_axes(), std::vector<unsigned>({0,2,3}));
  EXPECT_EQ(freezer.frozen_axes(), std::vector<unsigned>({1}));
  Eigen::VectorXd p1(1);
  p1 << 50.0;
  EXPECT_EQ(freezer.parameter_values(), p1);

  EXPECT_TRUE(freezer.freeze(0,0.2));
  EXPECT_EQ(freezer.num_domain_axes(), 2);
  EXPECT_EQ(freezer.num_parameters(), 2);
  EXPECT_EQ(freezer.num_domain_axes(), freezer.domain_axes().size());
  EXPECT_EQ(freezer.num_parameters(), freezer.parameters().size());
  EXPECT_EQ(freezer.free_axes(), std::vector<unsigned>({2,3}));
  EXPECT_EQ(freezer.frozen_axes(), std::vector<unsigned>({0,1}));
  Eigen::VectorXd p01(2);
  p01 << 0.2, 50.0;
  EXPECT_EQ(freezer.parameter_values(), p01);

  EXPECT_TRUE(freezer.freeze(3,35.0));
  EXPECT_EQ(freezer.num_domain_axes(), 1);
  EXPECT_EQ(freezer.num_parameters(), 3);
  EXPECT_EQ(freezer.num_domain_axes(), freezer.domain_axes().size());
  EXPECT_EQ(freezer.num_parameters(), freezer.parameters().size());
  EXPECT_EQ(freezer.free_axes(), std::vector<unsigned>({2}));
  EXPECT_EQ(freezer.frozen_axes(), std::vector<unsigned>({0,1,3}));
  Eigen::VectorXd p013(3);
  p013 << 0.2, 50.0, 35.0;
  EXPECT_EQ(freezer.parameter_values(), p013);
  
  EXPECT_FALSE(freezer.freeze(3,35.0));
  EXPECT_EQ(freezer.num_domain_axes(), 1);
  EXPECT_EQ(freezer.num_parameters(), 3);
  EXPECT_EQ(freezer.num_domain_axes(), freezer.domain_axes().size());
  EXPECT_EQ(freezer.num_parameters(), freezer.parameters().size());
  EXPECT_EQ(freezer.free_axes(), std::vector<unsigned>({2}));
  EXPECT_EQ(freezer.frozen_axes(), std::vector<unsigned>({0,1,3}));
  EXPECT_EQ(freezer.parameter_values(), p013);

  EXPECT_TRUE(freezer.freeze(2,100.0));
  EXPECT_EQ(freezer.num_domain_axes(), 0);
  EXPECT_EQ(freezer.num_parameters(), 4);
  EXPECT_EQ(freezer.num_domain_axes(), freezer.domain_axes().size());
  EXPECT_EQ(freezer.num_parameters(), freezer.parameters().size());
  EXPECT_EQ(freezer.free_axes(), std::vector<unsigned>({}));
  EXPECT_EQ(freezer.frozen_axes(), std::vector<unsigned>({0,1,2,3}));
  Eigen::VectorXd p0123(4);
  p0123 << 0.2, 50.0, 100.0, 35.0;
  EXPECT_EQ(freezer.parameter_values(), p0123);

  EXPECT_TRUE(freezer.thaw(1));
  EXPECT_EQ(freezer.num_domain_axes(), 1);
  EXPECT_EQ(freezer.num_parameters(), 3);
  EXPECT_EQ(freezer.num_domain_axes(), freezer.domain_axes().size());
  EXPECT_EQ(freezer.num_parameters(), freezer.parameters().size());
  EXPECT_EQ(freezer.free_axes(), std::vector<unsigned>({1}));
  EXPECT_EQ(freezer.frozen_axes(), std::vector<unsigned>({0,2,3}));
  Eigen::VectorXd p023(3);
  p023 << 0.2, 100.0, 35.0;
  EXPECT_EQ(freezer.parameter_values(), p023);

  EXPECT_FALSE(freezer.thaw(1));
  EXPECT_EQ(freezer.num_domain_axes(), 1);
  EXPECT_EQ(freezer.num_parameters(), 3);
  EXPECT_EQ(freezer.num_domain_axes(), freezer.domain_axes().size());
  EXPECT_EQ(freezer.num_parameters(), freezer.parameters().size());
  EXPECT_EQ(freezer.free_axes(), std::vector<unsigned>({1}));
  EXPECT_EQ(freezer.frozen_axes(), std::vector<unsigned>({0,2,3}));
  EXPECT_EQ(freezer.parameter_values(), p023);
  
  EXPECT_TRUE(freezer.thaw(2));
  EXPECT_EQ(freezer.num_domain_axes(), 2);
  EXPECT_EQ(freezer.num_parameters(), 2);
  EXPECT_EQ(freezer.num_domain_axes(), freezer.domain_axes().size());
  EXPECT_EQ(freezer.num_parameters(), freezer.parameters().size());
  EXPECT_EQ(freezer.free_axes(), std::vector<unsigned>({1,2}));
  EXPECT_EQ(freezer.frozen_axes(), std::vector<unsigned>({0,3}));
  Eigen::VectorXd p03(2);
  p03 << 0.2, 35.0;
  EXPECT_EQ(freezer.parameter_values(), p03);

  EXPECT_TRUE(freezer.thaw(0));
  EXPECT_EQ(freezer.num_domain_axes(), 3);
  EXPECT_EQ(freezer.num_parameters(), 1);
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
  pdf_1d::TwoComponentPDF pdf_x(&exp_pdf, "exp", &gauss_pdf, "gauss");
  Eigen::VectorXd p(3);
  p << 0.5, 100.0, 35.0;
  Eigen::VectorXd dp(3);
  dp << 1e-6, 1e-6, 1e-6;
  function::PMAFReverser pdf(&pdf_x);
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
  pdf_1d::TwoComponentPDF pdf_x(&exp_pdf, "exp", &gauss_pdf, "gauss");
  Eigen::VectorXd p(3);
  p << 0.5, 100.0, 35.0;
  Eigen::VectorXd dp(3);
  dp << 1e-6, 1e-6, 1e-6;
  function::PMAFReverser pdf(&pdf_x);
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
  pdf_1d::TwoComponentPDF pdf_x(&exp_pdf, "exp", &gauss_pdf, "gauss");
  Eigen::VectorXd p(2);
  p << 40, 100.0;
  Eigen::VectorXd dp(2);
  dp << 1e-6, 1e-6;
  function::PMAFReverser pdf(&pdf_x);
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
  pdf_1d::TwoComponentPDF pdf_x(&exp_pdf, "exp", &gauss_pdf, "gauss");
  Eigen::VectorXd p(2);
  p << 40, 100.0;
  Eigen::VectorXd dp(2);
  dp << 1e-6, 1e-6;
  function::PMAFReverser pdf(&pdf_x);
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

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
