/*
   calin/calib/pmt_ses_models.cpp -- Stephen Fegan -- 2017-04-24

   PMT single-electron spectrum models

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

#include <limits>

#include "calib/pmt_ses_models.hpp"

using namespace calin::calib::pmt_ses_models;

static constexpr double inf = std::numeric_limits<double>::infinity();

TwoGaussianSES::TwoGaussianSES(double dx):
  pattern::delegation::Delegator<ParameterizableSingleAxisFunction>(
    new calin::math::pdf_1d::TwoComponent1DPDF(
      new calin::math::pdf_1d::LimitedGaussianPDF(0, inf, dx), "gauss_lo",
      new calin::math::pdf_1d::LimitedGaussianPDF(0, inf, dx), "gauss_hi",
      true, true), true),
  calin::math::function::BasicReducedSpaceParameterizable<calin::math::function::ParameterizableSingleAxisFunction>(delegate_)
{
  remove_parameter_from_subspace(1, 0.0); // Mean of gauss-lo
}

TwoGaussianSES::~TwoGaussianSES()
{
  // nothing to see here
}

TwoGaussianSESConstrained::TwoGaussianSESConstrained(double dx):
  pattern::delegation::Delegator<ParameterizableSingleAxisFunction>(
    new calin::math::pdf_1d::TwoComponent1DConstraintPDF(
      new calin::math::pdf_1d::LimitedGaussianPDF(0, inf, dx), "gauss_lo",
      new calin::math::pdf_1d::LimitedGaussianPDF(0, inf, dx), "gauss_hi",
      true, true), true),
  calin::math::function::BasicReducedSpaceParameterizable<calin::math::function::ParameterizableSingleAxisFunction>(delegate_)
{
  remove_parameter_from_subspace(1, 0.0); // Mean of gauss-lo
}

TwoGaussianSESConstrained::~TwoGaussianSESConstrained()
{
  // nothing to see here
}

TwoGaussianSESConstrained_Fast::TwoGaussianSESConstrained_Fast(double dx):
  pattern::delegation::Delegator<ParameterizableSingleAxisFunction>(
    new calin::math::pdf_1d::TwoComponent1DConstraintPDF(
      new calin::math::pdf_1d::LimitedGaussianPDF(0, inf, dx), "gauss_lo",
      new calin::math::pdf_1d::LimitedGaussianPDF(0, inf, dx), "gauss_hi",
      true, true, true), true),
  calin::math::function::BasicReducedSpaceParameterizable<calin::math::function::ParameterizableSingleAxisFunction>(delegate_)
{
  remove_parameter_from_subspace(1, 0.0); // Mean of gauss-lo
}

TwoGaussianSESConstrained_Fast::~TwoGaussianSESConstrained_Fast()
{
  // nothing to see here
}
