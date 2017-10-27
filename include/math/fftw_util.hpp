/*

   calin/math/fftw_util.hpp -- Stephen Fegan -- 2016-03-29

   Utility functions for fftw HC types

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include<math/fftw_util.pb.h>

namespace calin { namespace math { namespace fftw_util {

#ifndef SWIG

void hcvec_scale_and_multiply(double* ovec, const double* ivec1,
  const double* ivec2, unsigned nsample, double scale = 1.0);

void hcvec_scale_and_multiply_conj(double* ovec, const double* ivec1,
  const double* ivec2_conj, unsigned nsample, double scale = 1.0);

void hcvec_scale_and_add(double* ovec, const double* ivec, unsigned nsample,
  double scale = 1.0);

#endif

int proto_planning_enum_to_fftw_flag(calin::ix::math::fftw_util::FFTWPlanningRigor x);

bool load_wisdom_from_file(std::string filename = "~/.calin_fft_wisdom");
bool load_wisdom_from_proto(const calin::ix::math::fftw_util::FFTWWisdom& proto);

bool save_wisdom_to_file(std::string filename = "~/.calin_fft_wisdom");
bool save_wisdom_to_proto(calin::ix::math::fftw_util::FFTWWisdom& proto);

} } } // namespace calin::math::fftw_util
