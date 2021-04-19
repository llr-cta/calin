/*

   calin/simulation/sct_facet_scheme.hpp -- Stephen Fegan -- 2021-04-16

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

#pragma once

#include <cmath>

namespace calin { namespace simulation { namespace sct_optics {

constexpr double COS_PI_4 = M_SQRT1_2;
constexpr double SIN_PI_4 = M_SQRT1_2;
constexpr double COS_PI_8 = 9.2387953251128674e-01;
constexpr double SIN_PI_8 = 3.8268343236508978e-01;
constexpr double COS_PI_16 = 9.8078528040323043e-01;
constexpr double SIN_PI_16 = 1.9509032201612825e-01;
constexpr double COS_PI_32 = 9.9518472667219693e-01;
constexpr double SIN_PI_32 = 9.8017140329560604e-02;

class SCTFacetScheme
{
public:
  virtual ~SCTFacetScheme();
  virtual int find_facet(double x, double z) = 0;
};

class SCTPrimaryFacetScheme: public SCTFacetScheme
{
public:
  SCTPrimaryFacetScheme(): SCTFacetScheme() { };
  SCTPrimaryFacetScheme(double r1i, double r1o, double r2i, double r2o, double gap_2):
    SCTFacetScheme(), r1i_(r1i), r1o_(r1o), r2i_(r2i), r2o_(r2o), gap_2_(gap_2) { }
  virtual ~SCTPrimaryFacetScheme();
  virtual int find_facet(double x, double z);
private:
  double r1i_ = 219.350*COS_PI_16;        // inner radius of inner panel (flat edge distance from origin)
  double r1o_ = 340.000*COS_PI_32 - 0.7; // outer radius of inner panel
  double r2i_ = 340.000*COS_PI_32 + 0.7; // inner radius of outer panel
  double r2o_ = 483.1875*COS_PI_32;      // outer radius of outer panel
  double gap_2_ = 0.7;                   // half the gap between panels along edge
};

class SCTSecondaryFacetScheme: public SCTFacetScheme
{
public:
  SCTSecondaryFacetScheme(): SCTFacetScheme() { };
  SCTSecondaryFacetScheme(double r1i, double r1o, double r2i, double r2o, double gap_2):
    SCTFacetScheme(), r1i_(r1i), r1o_(r1o), r2i_(r2i), r2o_(r2o), gap_2_(gap_2) { }
  virtual ~SCTSecondaryFacetScheme();
  virtual int find_facet(double x, double z);
private:
  double r1i_ = 39.45*COS_PI_8;        // inner radius of inner panel (flat edge distance from origin)
  double r1o_ = 159.65*COS_PI_16 - 0.7; // outer radius of inner panel
  double r2i_ = 159.65*COS_PI_16 + 0.7; // inner radius of outer panel
  double r2o_ = 270.83*COS_PI_16;       // outer radius of outer panel
  double gap_2_ = 0.7;                  // half the gap between panels along edge
};

} } } // namespace calin::simulations::sct_optics
