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

constexpr double COS_PI = -1;
constexpr double SIN_PI = 0;
constexpr double COS_PI_2 = 0;
constexpr double SIN_PI_2 = 1;
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
  virtual double facet_area(int ifacet) = 0;
  virtual bool facet_centroid(int ifacet, double& x_out, double& z_out) = 0;
};

class SCTPrimaryFacetScheme: public SCTFacetScheme
{
public:
  SCTPrimaryFacetScheme(): SCTFacetScheme() { };
  SCTPrimaryFacetScheme(double r1i, double r1o, double r2i, double r2o, double gap_2):
    SCTFacetScheme(), r1i_(r1i), r1o_(r1o), r2i_(r2i), r2o_(r2o), gap_2_(gap_2) { }
  virtual ~SCTPrimaryFacetScheme();
  virtual int find_facet(double x, double z) override;
  virtual double facet_area(int ifacet) override;
  virtual bool facet_centroid(int ifacet, double& x_out, double& z_out) override;

  static double default_r1i() { return 219.350*COS_PI_16; }
  static double default_r1o() { return 340.000*COS_PI_32 - 0.7; }
  static double default_r2i() { return 340.000*COS_PI_32 + 0.7; }
  static double default_r2o() { return 483.1875*COS_PI_32; }
  static double default_gap_2() { return 0.7; }

private:
  double r1i_ = default_r1i(); // inner radius of inner panel (flat edge distance from origin)
  double r1o_ = default_r1o(); // outer radius of inner panel
  double r2i_ = default_r2i(); // inner radius of outer panel
  double r2o_ = default_r2o(); // outer radius of outer panel
  double gap_2_ = default_gap_2(); // half the gap between panels along edge
};

class SCTSecondaryFacetScheme: public SCTFacetScheme
{
public:
  SCTSecondaryFacetScheme(): SCTFacetScheme() { };
  SCTSecondaryFacetScheme(double r1i, double r1o, double r2i, double r2o, double gap_2):
    SCTFacetScheme(), r1i_(r1i), r1o_(r1o), r2i_(r2i), r2o_(r2o), gap_2_(gap_2) { }
  virtual ~SCTSecondaryFacetScheme();
  virtual int find_facet(double x, double z) override;
  virtual double facet_area(int ifacet) override;
  virtual bool facet_centroid(int ifacet, double& x_out, double& z_out) override;

  static double default_r1i() { return 39.45*COS_PI_8; }
  static double default_r1o() { return 159.65*COS_PI_16 - 0.7; }
  static double default_r2i() { return 159.65*COS_PI_16 + 0.7; }
  static double default_r2o() { return 270.83*COS_PI_16; }
  static double default_gap_2() { return 0.7; }

private:
  double r1i_ = default_r1i(); // inner radius of inner panel (flat edge distance from origin)
  double r1o_ = default_r1o(); // outer radius of inner panel
  double r2i_ = default_r2i(); // inner radius of outer panel
  double r2o_ = default_r2o(); // outer radius of outer panel
  double gap_2_ = default_gap_2(); // half the gap between panels along edge
};

} } } // namespace calin::simulations::sct_optics
