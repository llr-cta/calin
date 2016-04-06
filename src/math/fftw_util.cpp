/*

   calin/math/fftw_util.cpp -- Stephen Fegan -- 2016-03-29

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

#include <math/fftw_util.hpp>

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
      (*ro++) = ((*ri1)*(*ri2) - (*ci1)*(*ci2)) * scale;
      (*co--) = ((*ri1++)*(*ci2--) + (*ci1--)*(*ri2++)) * scale;
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
