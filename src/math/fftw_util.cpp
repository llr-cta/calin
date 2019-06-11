/*

   calin/math/fftw_util.cpp -- Stephen Fegan -- 2016-03-29

   Utility functions for fftw HC types

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

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

#include <util/file.hpp>
#include <math/fftw_util.hpp>

#if 0
void calin::math::fftw_util::
hcvec_scale_and_multiply(double* ovec, const double* ivec1,
                         const double* ivec2, unsigned nsample, double scale)
{
  double *ro = ovec;
  double *co = ovec + nsample-1;
  const double *ri1 = ivec1;
  const double *ci1 = ivec1 + nsample-1;
  const double *ri2 = ivec2;
  const double *ci2 = ivec2 + nsample-1;

  (*ro++) = (*ri1++) * (*ri2++) * scale;
  if(ro==ri1 or ro==ri2)
  {
    while(ro < co)
    {
      double vri1 = *ri1++;
      double vci1 = *ci1--;
      double vri2 = *ri2++;
      double vci2 = *ci2--;
      (*ro++) = (vri1*vri2 - vci1*vci2)*scale;
      (*co--) = (vri1*vci2 + vci1*vri2)*scale;
    }
   }
  else
  {
    while(ro < co)
    {
      (*ro++) = ((*ri1)*(*ri2) - (*ci1)*(*ci2)) * scale;
      (*co--) = ((*ri1++)*(*ci2--) + (*ci1--)*(*ri2++)) * scale;
    }
  }
  if(ro==co)(*ro) = (*ri1) * (*ri2) * scale;
}

void calin::math::fftw_util::
hcvec_scale_and_multiply_conj(double* ovec, const double* ivec1,
                       const double* ivec2_conj, unsigned nsample, double scale)
{
  double *ro = ovec;
  double *co = ovec + nsample-1;
  const double *ri1 = ivec1;
  const double *ci1 = ivec1 + nsample-1;
  const double *ri2 = ivec2_conj;
  const double *ci2 = ivec2_conj + nsample-1;

  (*ro++) = (*ri1++) * (*ri2++) * scale;
  if(ro==ri1 or ro==ri2)
  {
    while(ro < co)
    {
      double vri1 = *ri1++;
      double vci1 = *ci1--;
      double vri2 = *ri2++;
      double vci2 = *ci2--;
      (*ro++) = (vri1*vri2 + vci1*vci2)*scale;
      (*co--) = (vci1*vri2 - vri1*vci2)*scale;
    }
   }
  else
  {
    while(ro < co)
    {
      (*ro++) = ((*ri1)*(*ri2) + (*ci1)*(*ci2)) * scale;
      (*co--) = ((*ci1--)*(*ri2++) - (*ri1++)*(*ci2--)) * scale;
    }
  }
  if(ro==co)(*ro) = (*ri1) * (*ri2) * scale;
}

void calin::math::fftw_util::
hcvec_scale_and_add(double* ovec, const double* ivec, unsigned nsample,
  double scale)
{
  double *ro = ovec;
  double *re = ovec + nsample;
  const double *ri = ivec;
  while(ro<re)
    *(ro++) += *(ri++)*scale;
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
