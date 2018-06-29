/*

   calin/calib/pmt_ses_models.hpp -- Stephen Fegan -- 2017-04-24

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

#pragma once

#include "math/function.hpp"
#include "math/pdf_1d.hpp"

namespace calin { namespace calib { namespace pmt_ses_models {

class TwoGaussianSES:
  public calin::math::function::ReducedSpaceParameterizableSingleAxisFunction
{
public:
  TwoGaussianSES(double dx = 0);
  virtual ~TwoGaussianSES();
};

class TwoGaussianSESConstrained:
  public calin::math::function::ReducedSpaceParameterizableSingleAxisFunction
{
public:
  TwoGaussianSESConstrained(double dx = 0);
  virtual ~TwoGaussianSESConstrained();
};

class TwoGaussianSESConstrained_Fast:
  public calin::math::function::ReducedSpaceParameterizableSingleAxisFunction
{
public:
  TwoGaussianSESConstrained_Fast(double dx = 0);
  virtual ~TwoGaussianSESConstrained_Fast();
};

} } } // namespace calin::calib::pmt_ses_models
