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
#include <stdexcept>

#include <simulation/sct_facet_scheme.hpp>

using namespace calin::simulation::sct_optics;

SCTFacetScheme::~SCTFacetScheme()
{
  // nothing to see here
}

Eigen::VectorXi SCTFacetScheme::
bulk_find_facet(const Eigen::VectorXd& x, const Eigen::VectorXd& z)
{
  if(x.size() != z.size()) {
    throw std::runtime_error("bulk_find_facet: X and Z must have same size");
  }
  Eigen::VectorXi id(z.size());
  for(unsigned i=0;i<id.size(); ++i) {
    id(i) = this->find_facet(x(i), z(i));
  }
  return id;
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

  double x2 = -y1;
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

unsigned SCTPrimaryFacetScheme::num_facets()
{
  return 48;
}

double SCTPrimaryFacetScheme::facet_area(int ifacet)
{
  if(ifacet<0) {
    return 0;
  } else if(ifacet<16) {
    const double r0 = gap_2_/SIN_PI_16;
    const double ri = std::max(r1i_ - r0, 0.0);
    const double ro = std::max(r1o_ - r0*COS_PI_32, 0.0);
    return ro*ro*SIN_PI_16/COS_PI_32/COS_PI_32
     - ri*ri*SIN_PI_16/COS_PI_16;
  }  else if(ifacet<48) {
    const double r0 = gap_2_/SIN_PI_32;
    const double ri = std::max(r2i_ - r0, 0.0);
    const double ro = std::max(r2o_ - r0, 0.0);
    return (ro*ro - ri*ri)*SIN_PI_32/COS_PI_32;
  }
  return 0;
}

namespace {
  void rotate_in_place(double& x, double&y, double C, double S) {
    double x_new = C*x + S*y;
    y = C*y - S*x;
    x = x_new;
  }
}

bool SCTPrimaryFacetScheme::facet_centroid(int ifacet, double& x, double& y)
{
  if(ifacet<0) {
    return false;
  } else if(ifacet < 16) {
    const double r0 = gap_2_/SIN_PI_16;
    const double ri = std::max(r1i_ - r0, 0.0);
    const double ro = std::max(r1o_ - r0*COS_PI_32, 0.0);
    x = 0;
    y = 2.0*(ro*ro*ro/COS_PI_32 - ri*ri*ri/COS_PI_16)/
      (3.0*(ro*ro/COS_PI_32/COS_PI_32 - ri*ri/COS_PI_16)) + r0;
    rotate_in_place(x, y, COS_PI_16, SIN_PI_16);
    if((ifacet & 1) == 0)x = -x;
    rotate_in_place(x, y, COS_PI_8, SIN_PI_8);
    if((ifacet & 2) == 0)x = -x;
    rotate_in_place(x, y, COS_PI_4, SIN_PI_4);
    if((ifacet & 4) == 0)x = -x;
    rotate_in_place(x, y, COS_PI_2, SIN_PI_2);
    if((ifacet & 8) == 0)x = -x;
    return true;
  } else if(ifacet < 48) {
    const double r0 = gap_2_/SIN_PI_32;
    const double ri = std::max(r2i_ - r0, 0.0);
    const double ro = std::max(r2o_ - r0, 0.0);
    x = 0;
    // y = 0.5*(r2i_ + r2o_);
    y = 2.0*(ro*ro*ro-ri*ri*ri)/(3.0*(ro*ro-ri*ri)) + r0;
    ifacet -= 16;
    rotate_in_place(x, y, COS_PI_32, SIN_PI_32);
    if((ifacet & 1) == 0)x = -x;
    rotate_in_place(x, y, COS_PI_16, SIN_PI_16);
    if((ifacet & 2) == 0)x = -x;
    rotate_in_place(x, y, COS_PI_8, SIN_PI_8);
    if((ifacet & 4) == 0)x = -x;
    rotate_in_place(x, y, COS_PI_4, SIN_PI_4);
    if((ifacet & 8) == 0)x = -x;
    rotate_in_place(x, y, COS_PI_2, SIN_PI_2);
    if((ifacet & 16) == 0)x = -x;
    return true;
  }
  return false;
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

  double x2 = -y1;
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

unsigned SCTSecondaryFacetScheme::num_facets()
{
  return 24;
}

double SCTSecondaryFacetScheme::facet_area(int ifacet)
{
  if(ifacet<0) {
    return 0;
  } else if(ifacet<8) {
    const double r0 = gap_2_/SIN_PI_8;
    const double ri = std::max(r1i_ - r0, 0.0);
    const double ro = std::max(r1o_ - r0*COS_PI_16, 0.0);
    return ro*ro*SIN_PI_8/COS_PI_16/COS_PI_16
     - ri*ri*SIN_PI_8/COS_PI_8;
  }  else if(ifacet<24) {
    const double r0 = gap_2_/SIN_PI_16;
    const double ri = std::max(r2i_ - r0, 0.0);
    const double ro = std::max(r2o_ - r0, 0.0);
    return (ro*ro - ri*ri)*SIN_PI_16/COS_PI_16;
  }
  return 0;
}

bool SCTSecondaryFacetScheme::facet_centroid(int ifacet, double& x, double& y)
{
  if(ifacet<0) {
    return false;
  } else if(ifacet < 8) {
    const double r0 = gap_2_/SIN_PI_8;
    const double ri = std::max(r1i_ - r0, 0.0);
    const double ro = std::max(r1o_ - r0*COS_PI_16, 0.0);
    x = 0;
    y = 2.0*(ro*ro*ro/COS_PI_16 - ri*ri*ri/COS_PI_8)/
      (3.0*(ro*ro/COS_PI_16/COS_PI_16 - ri*ri/COS_PI_8)) + r0;
    rotate_in_place(x, y, COS_PI_8, SIN_PI_8);
    if((ifacet & 1) == 0)x = -x;
    rotate_in_place(x, y, COS_PI_4, SIN_PI_4);
    if((ifacet & 2) == 0)x = -x;
    rotate_in_place(x, y, COS_PI_2, SIN_PI_2);
    if((ifacet & 4) == 0)x = -x;
    return true;
  } else if(ifacet < 24) {
    const double r0 = gap_2_/SIN_PI_16;
    const double ri = std::max(r2i_ - r0, 0.0);
    const double ro = std::max(r2o_ - r0, 0.0);
    x = 0;
    y = 2.0*(ro*ro*ro-ri*ri*ri)/(3.0*(ro*ro-ri*ri)) + r0;
    ifacet -= 8;
    rotate_in_place(x, y, COS_PI_16, SIN_PI_16);
    if((ifacet & 1) == 0)x = -x;
    rotate_in_place(x, y, COS_PI_8, SIN_PI_8);
    if((ifacet & 2) == 0)x = -x;
    rotate_in_place(x, y, COS_PI_4, SIN_PI_4);
    if((ifacet & 4) == 0)x = -x;
    rotate_in_place(x, y, COS_PI_2, SIN_PI_2);
    if((ifacet & 8) == 0)x = -x;
    return true;
  }
  return false;
}
