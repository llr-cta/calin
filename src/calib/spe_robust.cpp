/*
   calin/calib/spe_fit_robust.cpp -- Stephen Fegan -- 2017-04-19

   Functions to do fit to multi-electron spectrum in the "single PE"
   domain. Robust estimation.

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
#include <cmath>
#include <limits>
#include <cassert>
#include <memory>

#include <math/special.hpp>
#include <calib/spe_fit.hpp>

using namespace calin::calib::spe_fit;
using calin::math::special::SQR;
using calin::math::function::assign_parameters;
using calin::calib::spe_fit;

MESSigAdapter::MESSigAdapter(MultiElectronSpectrum* mes):
  Parameterizable1DPDF(), mes_(mes)
{
  // nothing to see here
}

MESSigAdapter::~MESSigAdapter()
{
  // nothing to see here
}

unsigned MESSigAdapter::num_parameters()
{
  return mes_->num_parameters()
}

std::vector<function::ParameterAxis> MESSigAdapter::parameters()
{
  return mes_->parameters();
}

Eigen::VectorXd MESSigAdapter::parameter_values()
{
  return mes_->parameter_values();
}

void MESSigAdapter::set_parameter_values(ConstVecRef values)
{
  // do nothing
}

function::DomainAxis MESSigAdapter::domain_axis()
{
  return { "charge", "", false, 0, false, 0 };
}

bool MESSigAdapter::can_calculate_gradient()
{
  return false;
}

bool MESSigAdapter::can_calculate_hessian()
{
  return false;
}

bool can_calculate_parameter_gradient() override;
bool can_calculate_parameter_hessian() override;

  double value_1d(double x) override;
  double value_and_gradient_1d(double x,  double& dfdx) override;
  double value_gradient_and_hessian_1d(double x, double& dfdx,
                            double& d2fdx2) override;
  double value_and_parameter_gradient_1d(double x,
                                         VecRef gradient) override;
  double value_parameter_gradient_and_hessian_1d(double x, VecRef gradient,
                                              MatRef hessian) override;

  double error_up() override;

  bool can_calculate_mean_and_variance() override;
  void mean_and_variance(double& mean, double& var) override;
protected:
  MultiElectronSpectrum* mes_ = nullptr;
};
