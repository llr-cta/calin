/*

   calin/math/fftw_util.hpp -- Stephen Fegan -- 2016-03-29

   Utility functions for fftw HC types

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#pragma once

namespace calin { namespace math { namespace fftw_util {

void hcvec_scale_and_multiply(double* ovec, const double* ivec1,
  const double* ivec2, unsigned nsample, double scale = 1.0);

void hcvec_scale_and_multiply_conj(double* ovec, const double* ivec1,
  const double* ivec2_conj, unsigned nsample, double scale = 1.0);

void hcvec_scale_and_add(double* ovec, const double* ivec, unsigned nsample,
  double scale = 1.0);


} } } // namespace calin::math::fftw_util
