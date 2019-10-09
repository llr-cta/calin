/*

   calin/math/fftw_util.cpp -- Stephen Fegan -- 2016-03-29

   Utility functions for fftw HC types

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <cstdlib>

#include <fftw3.h>

#include <util/vcl.hpp>
#include <util/file.hpp>
#include <math/fftw_util.hpp>

#if INSTRSET >= 7
void calin::math::fftw_util::hcvec_multiply_and_add_real(double* ovec, const double* ivec1,
  const double* ivec2, double real_addand, unsigned nsample)
{
  hcvec_multiply_and_add_real_vcl<calin::util::vcl::VCLDoubleReal<calin::util::vcl::VCL256Architecture> >(ovec,
    ivec1, ivec2, real_addand, nsample);
  // hcvec_multiply_and_add_real<double>(ovec, ivec1, ivec2, real_addand, nsample);
}

void calin::math::fftw_util::hcvec_polynomial(double* ovec, const double* ivec,
  const std::vector<double>& p, unsigned nsample)
{
  hcvec_polynomial_vcl<calin::util::vcl::VCLDoubleReal<calin::util::vcl::VCL256Architecture> >(ovec, ivec, p, nsample);
  //hcvec_polynomial<double>(ovec, ivec, p, nsample);
}

#endif

bool calin::math::fftw_util::load_wisdom_from_file(std::string filename)
{
  calin::util::file::expand_filename_in_place(filename);
  return fftw_import_wisdom_from_filename(filename.c_str()) != 0;
}

bool calin::math::fftw_util::
load_wisdom_from_proto(const calin::ix::math::fftw_util::FFTWWisdom& proto)
{
  return fftw_import_wisdom_from_string(proto.wisdom().c_str()) != 0;
}

bool calin::math::fftw_util::save_wisdom_to_file(std::string filename)
{
  calin::util::file::expand_filename_in_place(filename);
  return fftw_export_wisdom_to_filename(filename.c_str()) != 0;
}

bool calin::math::fftw_util::
save_wisdom_to_proto(calin::ix::math::fftw_util::FFTWWisdom& proto)
{
  char* data = fftw_export_wisdom_to_string();
  if(!data)return false;
  proto.set_wisdom(data);
  std::free(data);
  return true;
}

int calin::math::fftw_util::
proto_planning_enum_to_fftw_flag(calin::ix::math::fftw_util::FFTWPlanningRigor x)
{
  switch(x) {
    case calin::ix::math::fftw_util::ESTIMATE:    return FFTW_ESTIMATE;
    case calin::ix::math::fftw_util::MEASURE:     return FFTW_MEASURE;
    case calin::ix::math::fftw_util::PATIENT:     return FFTW_PATIENT;
    case calin::ix::math::fftw_util::EXHAUSTIVE:  return FFTW_EXHAUSTIVE;
    default: throw std::invalid_argument("Unknown planning rigor type.");
  }
}
