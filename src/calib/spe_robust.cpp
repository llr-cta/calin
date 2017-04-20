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

using namespace calin;
using namespace calin::calib::spe_fit;
using calin::math::special::SQR;
using calin::math::function::assign_parameters;
using namespace calin::calib::spe_fit;
using namespace calin::math::m_estimate;
using namespace calin::math::data_modeling;

namespace {

class MESSigAdapter : public calin::math::pdf_1d::Parameterizable1DPDF
{
public:
  MESSigAdapter(MultiElectronSpectrum* mes);
  virtual ~MESSigAdapter();

  // Reiterate functions from ParameterizableSingleAxisFunction

  unsigned num_parameters() override;
  std::vector<calin::math::function::ParameterAxis> parameters() override;
  Eigen::VectorXd parameter_values() override;
  void set_parameter_values(ConstVecRef values) override;

  calin::math::function::DomainAxis domain_axis() override;

  bool can_calculate_gradient() override;
  bool can_calculate_hessian() override;
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

class MESPedAdapter : public MESSigAdapter
{
public:
  using MESSigAdapter::MESSigAdapter;
  virtual ~MESPedAdapter();

  double value_1d(double x) override;
  double value_and_parameter_gradient_1d(double x,
                                         VecRef gradient) override;
  double value_parameter_gradient_and_hessian_1d(double x, VecRef gradient,
                                              MatRef hessian) override;
};

} // anonoymous namespace

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
  return mes_->num_parameters();
}

std::vector<calin::math::function::ParameterAxis> MESSigAdapter::parameters()
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

calin::math::function::DomainAxis MESSigAdapter::domain_axis()
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

bool MESSigAdapter::can_calculate_parameter_gradient()
{
  return mes_->can_calculate_parameter_gradient();
}

bool MESSigAdapter::can_calculate_parameter_hessian()
{
  return mes_->can_calculate_parameter_hessian();
}

double MESSigAdapter::value_and_gradient_1d(double x,  double& dfdx)
{
  assert(0);
  dfdx = 0;
  return 0;
}

double MESSigAdapter::
value_gradient_and_hessian_1d(double x, double& dfdx, double& d2fdx2)
{
  assert(0);
  d2fdx2 = 0;
  dfdx = 0;
  return 0;
}

double MESSigAdapter::value_1d(double x)
{
  return mes_->pdf_mes(x);
}

double MESSigAdapter::
value_and_parameter_gradient_1d(double x, VecRef gradient)
{
  return mes_->pdf_gradient_mes(x, gradient);
}

double MESSigAdapter::
value_parameter_gradient_and_hessian_1d(double x, VecRef gradient, MatRef hessian)
{
  return mes_->pdf_gradient_hessian_mes(x, gradient, hessian);
}

double MESSigAdapter::error_up()
{
  return 0;
}

bool MESSigAdapter::can_calculate_mean_and_variance()
{
  return false;
}

void MESSigAdapter::mean_and_variance(double& mean, double& var)
{
  assert(0);
  mean = 0;
  var = 0;
  return;
}

MESPedAdapter::~MESPedAdapter()
{
  // nothing to see here
}

double MESPedAdapter::value_1d(double x)
{
  return mes_->pdf_ped(x);
}

double MESPedAdapter::
value_and_parameter_gradient_1d(double x, VecRef gradient)
{
  return mes_->pdf_gradient_ped(x, gradient);
}

double MESPedAdapter::
value_parameter_gradient_and_hessian_1d(double x, VecRef gradient, MatRef hessian)
{
  return mes_->pdf_gradient_hessian_ped(x, gradient, hessian);
}

// -----------------------------------------------------------------------------
// *****************************************************************************
// -----------------------------------------------------------------------------
//
// SPERobust
//
// -----------------------------------------------------------------------------
// *****************************************************************************
// -----------------------------------------------------------------------------

SPERobust::SPERobust(MultiElectronSpectrum& mes_model,
    const calin::math::histogram::SimpleHist& mes_data,
    LikelihoodRhoFunction* rho, bool adopt_rho):
  calin::math::function::MultiAxisFunction(),
  mes_model_(&mes_model), npar_(mes_model_->num_parameters()),
  rho_(rho==nullptr ? new NullLikelihoodRhoFunction() : rho),
  adopt_rho_(rho==nullptr ? true : adopt_rho),
  mes_cost_(new IID1DDataMEstimateLikelihoodFunction(new MESPedAdapter(mes_model_)),
    rho_, mes_data, true, false)
{
  // nothing to see here
}

SPERobust::SPERobust(MultiElectronSpectrum& mes_model,
    const calin::math::histogram::SimpleHist& mes_data,
    const calin::math::histogram::SimpleHist& ped_data,
    LikelihoodRhoFunction* rho, bool adopt_rho):
  calin::math::function::MultiAxisFunction(),
  mes_model_(&mes_model), npar_(mes_model_->num_parameters()),
  rho_(rho==nullptr ? new NullLikelihoodRhoFunction() : rho),
  adopt_rho_(rho==nullptr ? true : adopt_rho),
  mes_cost_(new IID1DDataMEstimateLikelihoodFunction(new MESPedAdapter(mes_model_)),
    rho_, mes_data, true, false),
  ped_cost_(new IID1DDataMEstimateLikelihoodFunction(new MESPedAdapter(mes_model_)),
    rho_, ped_data, true, false)
{
  // nothing to see here
}

SPERobust::~SPERobust()
{
  delete mes_cost_;
  delete ped_cost_;
  if(adopt_rho_)delete rho_;
}

unsigned SPERobust::num_domain_axes()
{

}

std::vector<calin::math::function::DomainAxis> SPERobust::domain_axes()
{

}

double SPERobust::value(ConstVecRef x)
{

}

bool SPERobust::can_calculate_gradient()
{

}

double SPERobust::value_and_gradient(ConstVecRef x, VecRef gradient)
{

}

bool SPERobust::can_calculate_hessian()
{

}
double SPERobust::
value_gradient_and_hessian(ConstVecRef x, VecRef gradient, MatRef hessian)
{

}
