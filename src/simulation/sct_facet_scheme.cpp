/*

   calin/simulation/sct_facet_scheme.cpp -- Stephen Fegan -- 2021-04-18

   Class for mirror facet finding

   Copyright 2021, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <cmath>
#include <algorithm>

#include <simulation/sct_facet_scheme.hpp>

using namespace calin::simulation::sct_optics;

SCTFacetScheme::~SCTFacetScheme()
{
  // nothing to see here
}

SCTPrimaryFacetScheme::~SCTPrimaryFacetScheme()
{
  // nothing to see here
}

int SCTPrimaryFacetScheme::find_facet(double x, double y)
{
  unsigned sector = 0;

  double x1 = x;
  double y1 = y;
  sector |= x1>=0;
  x1 = std::abs(x1);

  double x2 = y1;
  double y2 = x1;
  sector = (sector<<1) | (x2>=0);
  x2 = std::abs(x2);

  double x4 = COS_PI_4*x2 - SIN_PI_4*y2;
  double y4 = COS_PI_4*y2 + SIN_PI_4*x2;
  sector = (sector<<1) | (x4>=0);
  x4 = std::abs(x4);

  double x8 = COS_PI_8*x4 - SIN_PI_8*y4;
  double y8 = COS_PI_8*y4 + SIN_PI_8*x4;
  sector = (sector<<1) | (x8>=0);
  x8 = std::abs(x8);

  double x16 = COS_PI_16*x8 - SIN_PI_16*y8;
  double y16 = COS_PI_16*y8 + SIN_PI_16*x8;
  sector = (sector<<1) | (x16>=0);
  x16 = std::abs(x16);

  // double x32 = COS_PI_32*x16 - SIN_PI_32*y16;
  double y32 = COS_PI_32*y16 + SIN_PI_32*x16;
  // x32 = std::abs(x32);

  if(y16>=r1i_ and y32<r1o_) {
    if(std::min(std::min(x1, x2), std::min(x4, x8))>=gap_2_) {
      return sector>>1;
    }
  } else if(y32>=r2i_ and y32<r2o_) {
    if(std::min(std::min(std::min(x1, x2), std::min(x4, x8)), x16)>=gap_2_) {
      return sector + 16;
    }
  }

  return -1;
}

SCTSecondaryFacetScheme::~SCTSecondaryFacetScheme()
{
  // nothing to see here
}

int SCTSecondaryFacetScheme::find_facet(double x, double y)
{
  unsigned sector = 0;

  double x1 = x;
  double y1 = y;
  sector |= x1>=0;
  x1 = std::abs(x1);

  double x2 = y1;
  double y2 = x1;
  sector = (sector<<1) | (x2>=0);
  x2 = std::abs(x2);

  double x4 = COS_PI_4*x2 - SIN_PI_4*y2;
  double y4 = COS_PI_4*y2 + SIN_PI_4*x2;
  sector = (sector<<1) | (x4>=0);
  x4 = std::abs(x4);

  double x8 = COS_PI_8*x4 - SIN_PI_8*y4;
  double y8 = COS_PI_8*y4 + SIN_PI_8*x4;
  sector = (sector<<1) | (x8>=0);
  x8 = std::abs(x8);

  // double x16 = COS_PI_16*x8 - SIN_PI_16*y8;
  double y16 = COS_PI_16*y8 + SIN_PI_16*x8;
  // x16 = std::abs(x16);

  if(y8>=r1i_ and y16<r1o_) {
    if(std::min(std::min(x1, x2), x4)>=gap_2_) {
      return sector>>1;
    }
  } else if(y16>=r2i_ and y16<r2o_) {
    if(std::min(std::min(x1, x2), std::min(x4, x8))>=gap_2_) {
      return sector + 8;
    }
  }

  return -1;
}
