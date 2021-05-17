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
#include <Eigen/Dense>
#include <simulation/sct_optics.pb.h>

namespace calin { namespace simulation { namespace sct_optics {

// from mpmath import mp
// mp.dps = 50 # far more digits than needed
// c = -1
// s = 0
// for i in range(0,6):
//     print("constexpr double COS_PI_%d = %s;"%(2**i,c))
//     print("constexpr double SIN_PI_%d = %s;"%(2**i,s))
//     s = mp.sqrt((1-c)/2)
//     c = mp.sqrt((1+c)/2)

constexpr double COS_PI_1 = -1.0;
constexpr double SIN_PI_1 = 0.0;
constexpr double COS_PI_2 = 0.0;
constexpr double SIN_PI_2 = 1.0;
constexpr double COS_PI_4 = 0.70710678118654752440084436210484903928483593768847;
constexpr double SIN_PI_4 = 0.70710678118654752440084436210484903928483593768847;
constexpr double COS_PI_8 = 0.92387953251128675612818318939678828682241662586364;
constexpr double SIN_PI_8 = 0.38268343236508977172845998403039886676134456248563;
constexpr double COS_PI_16 = 0.98078528040323044912618223613423903697393373089334;
constexpr double SIN_PI_16 = 0.19509032201612826784828486847702224092769161775195;
constexpr double COS_PI_32 = 0.99518472667219688624483695310947992157547486872986;
constexpr double SIN_PI_32 = 0.0980171403295606019941955638886418458611366731675;

class SCTFacetScheme
{
public:
  virtual ~SCTFacetScheme();
  virtual int find_facet(double x, double z) = 0;
  virtual unsigned num_facets() = 0;
  virtual double facet_area(int ifacet) = 0;
  virtual bool facet_centroid(int ifacet, double& x_out, double& z_out) = 0;
  virtual bool facet_vertices(int ifacet, Eigen::VectorXd& x_out, Eigen::VectorXd& z_out) = 0;
  virtual double inner_radius() = 0;
  virtual double outer_radius() = 0;
  Eigen::VectorXi bulk_find_facet(const Eigen::VectorXd& x, const Eigen::VectorXd& z);
#ifndef SWIG
  Eigen::Vector3d facet_centroid_3d(int ifacet, const double* p, unsigned pn);
#endif
  Eigen::Vector3d facet_centroid_3d(int ifacet, const Eigen::VectorXd p) {
    return facet_centroid_3d(ifacet, p.data(), p.size());
  }
};

class SCTPrimaryFacetScheme: public SCTFacetScheme
{
public:
  SCTPrimaryFacetScheme(): SCTFacetScheme() { };
  SCTPrimaryFacetScheme(double r1i, double r1o, double r2i, double r2o, double gap_2):
    SCTFacetScheme(), r1i_(r1i), r1o_(r1o), r2i_(r2i), r2o_(r2o), gap_2_(gap_2) { }
  SCTPrimaryFacetScheme(const calin::ix::simulation::sct_optics::SCTFacetScheme& scheme):
    SCTFacetScheme(),
    r1i_(scheme.inner_ring_inner_radius()), r1o_(scheme.inner_ring_outer_radius()),
    r2i_(scheme.outer_ring_inner_radius()), r2o_(scheme.outer_ring_outer_radius()),
    gap_2_(scheme.long_edge_half_gap()) { }
  virtual ~SCTPrimaryFacetScheme();
  virtual int find_facet(double x, double z) override;
  virtual unsigned num_facets() override;
  virtual double facet_area(int ifacet) override;
  virtual bool facet_centroid(int ifacet, double& x_out, double& z_out) override;
  virtual bool facet_vertices(int ifacet, Eigen::VectorXd& x_out, Eigen::VectorXd& z_out) override;
  virtual double inner_radius() override;
  virtual double outer_radius() override;

  static double default_r1i() { return 219.350*COS_PI_16; }
  static double default_r1o() { return 340.000*COS_PI_32 - 0.7; }
  static double default_r2i() { return 340.000*COS_PI_32 + 0.7; }
  static double default_r2o() { return 483.1875*COS_PI_32; }
  static double default_gap_2() { return 0.7; }

  static Eigen::VectorXi P1();
  static Eigen::VectorXi P2();

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
  SCTSecondaryFacetScheme(const calin::ix::simulation::sct_optics::SCTFacetScheme& scheme):
    SCTFacetScheme(),
    r1i_(scheme.inner_ring_inner_radius()), r1o_(scheme.inner_ring_outer_radius()),
    r2i_(scheme.outer_ring_inner_radius()), r2o_(scheme.outer_ring_outer_radius()),
    gap_2_(scheme.long_edge_half_gap()) { }
  virtual ~SCTSecondaryFacetScheme();
  virtual int find_facet(double x, double z) override;
  virtual unsigned num_facets() override;
  virtual double facet_area(int ifacet) override;
  virtual bool facet_centroid(int ifacet, double& x_out, double& z_out) override;
  virtual bool facet_vertices(int ifacet, Eigen::VectorXd& x_out, Eigen::VectorXd& z_out) override;
  virtual double inner_radius() override;
  virtual double outer_radius() override;

  static double default_r1i() { return 39.45*COS_PI_8; }
  static double default_r1o() { return 159.65*COS_PI_16 - 0.7; }
  static double default_r2i() { return 159.65*COS_PI_16 + 0.7; }
  static double default_r2o() { return 270.83*COS_PI_16; }
  static double default_gap_2() { return 0.7; }

  static Eigen::VectorXi S1();
  static Eigen::VectorXi S2();

private:
  double r1i_ = default_r1i(); // inner radius of inner panel (flat edge distance from origin)
  double r1o_ = default_r1o(); // outer radius of inner panel
  double r2i_ = default_r2i(); // inner radius of outer panel
  double r2o_ = default_r2o(); // outer radius of outer panel
  double gap_2_ = default_gap_2(); // half the gap between panels along edge
};

} } } // namespace calin::simulations::sct_optics
