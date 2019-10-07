/*

   calin/math/fftw_util.hpp -- Stephen Fegan -- 2016-03-29

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

#pragma once

#include<util/log.hpp>
#include<math/fftw_util.pb.h>

namespace calin { namespace math { namespace fftw_util {

#ifndef SWIG

template<typename T>
void hcvec_scale_and_multiply(T* ovec, const T* ivec1,
  const T* ivec2, unsigned nsample, T scale = 1.0)
{
  T *ro = ovec;
  T *co = ovec + nsample-1;
  const T *ri1 = ivec1;
  const T *ci1 = ivec1 + nsample-1;
  const T *ri2 = ivec2;
  const T *ci2 = ivec2 + nsample-1;

  (*ro++) = (*ri1++) * (*ri2++) * scale;
  if(ro==ri1 or ro==ri2)
  {
    while(ro < co)
    {
      T vri1 = *ri1++;
      T vci1 = *ci1--;
      T vri2 = *ri2++;
      T vci2 = *ci2--;
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

template<typename T>
void hcvec_scale_and_multiply_conj(T* ovec, const T* ivec1,
  const T* ivec2_conj, unsigned nsample, T scale = 1.0)
{
  T *ro = ovec;
  T *co = ovec + nsample-1;
  const T *ri1 = ivec1;
  const T *ci1 = ivec1 + nsample-1;
  const T *ri2 = ivec2_conj;
  const T *ci2 = ivec2_conj + nsample-1;

  (*ro++) = (*ri1++) * (*ri2++) * scale;
  if(ro==ri1 or ro==ri2)
  {
    while(ro < co)
    {
      T vri1 = *ri1++;
      T vci1 = *ci1--;
      T vri2 = *ri2++;
      T vci2 = *ci2--;
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

template<typename T>
void hcvec_scale(T* ovec, unsigned nsample,
  T scale = 1.0)
{
  T *ro = ovec;
  T *re = ovec + nsample;
  while(ro<re)
    *(ro++) *= scale;
}

template<typename T>
void hcvec_scale_and_add(T* ovec, const T* ivec, unsigned nsample,
  T scale = 1.0)
{
  T *ro = ovec;
  T *re = ovec + nsample;
  const T *ri = ivec;
  while(ro<re)
    *(ro++) += *(ri++)*scale;
}

template<typename T>
void hcvec_set_real(T* ovec, T real_value, unsigned nsample)
{
  T *ro = ovec;
  T *co = ovec + nsample-1;
  (*ro++) = real_value;
  while(ro < co)
  {
    (*ro++) = real_value;
    (*co--) = T(0);
  }
  if(ro==co)(*ro) = real_value;
}

template<typename T>
void hcvec_multiply_and_add_real(T* ovec, const T* ivec1,
  const T* ivec2, T real_addand, unsigned nsample, T scale = 1.0)
{
  T *ro = ovec;
  T *co = ovec + nsample-1;
  const T *ri1 = ivec1;
  const T *ci1 = ivec1 + nsample-1;
  const T *ri2 = ivec2;
  const T *ci2 = ivec2 + nsample-1;

  (*ro++) = (*ri1++) * (*ri2++) * scale + real_addand;
  if(ro==ri1 or ro==ri2)
  {
    while(ro < co)
    {
      T vri1 = *ri1++;
      T vci1 = *ci1--;
      T vri2 = *ri2++;
      T vci2 = *ci2--;
      (*ro++) = (vri1*vri2 - vci1*vci2) * scale + real_addand;
      (*co--) = (vri1*vci2 + vci1*vri2) * scale;
    }
   }
  else
  {
    while(ro < co)
    {
      (*ro++) = ((*ri1)*(*ri2) - (*ci1)*(*ci2)) * scale + real_addand;
      (*co--) = ((*ri1++)*(*ci2--) + (*ci1--)*(*ri2++)) * scale;
    }
  }
  if(ro==co)(*ro) = (*ri1) * (*ri2) * scale + real_addand;
}

template<typename T>
void hcvec_polynomial(T* ovec, const T* ivec, const std::vector<T>& p, unsigned nsample, T scale = 1.0)
{
  if(p.empty()) {
    hcvec_set_real(ovec, T(0), nsample);
    return;
  }

  auto pi = p.end();

  --pi;
  hcvec_set_real(ovec, *pi, nsample);

  while(pi != p.begin()) {
    --pi;
    hcvec_multiply_and_add_real(ovec, ovec, ivec, *pi, nsample, scale);
  }
}

#endif

int proto_planning_enum_to_fftw_flag(calin::ix::math::fftw_util::FFTWPlanningRigor x);

bool load_wisdom_from_file(std::string filename = "~/.calin_fft_wisdom");
bool load_wisdom_from_proto(const calin::ix::math::fftw_util::FFTWWisdom& proto);

bool save_wisdom_to_file(std::string filename = "~/.calin_fft_wisdom");
bool save_wisdom_to_proto(calin::ix::math::fftw_util::FFTWWisdom& proto);

} } } // namespace calin::math::fftw_util
