/*

   calin/simulation/atmosphere.hpp -- Stephen Fegan -- 2015-06-11

   Classes to model density and refractive index of atmosphere

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

// Originally: Atmposphere.hpp - Classes to handle atmosphere
// Stephen Fegan - sfegan@llr.in2p3.fr - July 2012
// Id: Atmosphere.hpp 5416 2013-06-24 13:46:00Z sfegan

#pragma once

#include<vector>
#include<cassert>

#include"calin_global_definitions.hpp"

#include <math/special.hpp>
#include <math/ray.hpp>
#include <math/spline_interpolation.hpp>
#include <simulation/atmosphere.pb.h>

// Units:
// Height, z:    cm
// Density, rho: g/cm^3
// Thickness, t: g/cm^2

namespace calin { namespace simulation { namespace atmosphere {

#if 0
class AtmComposition
{
 public:
  double c_N2   = 0.7807914486; // 0.78084;
  double c_O2   = 0.20947;
  double c_Ar   = 0.00934;
  double c_CO2  = 0.00037;
  double c_H20  = 5E-6;
  double c_Ne   = 0.00001818;
  double c_He   = 0.00000524;
  double c_O3   = 1.314E-7;

  double frac_N2() const { return c_N2/norm(); }
  double frac_O2() const { return c_O2/norm(); }
  double frac_Ar() const { return c_Ar/norm(); }
  double frac_CO2() const { return c_CO2/norm(); }
  double frac_H20() const { return c_H20/norm(); }
  double frac_Ne() const { return c_Ne/norm(); }
  double frac_He() const { return c_He/norm(); }
  double frac_O3() const { return c_O3/norm(); }

  double norm() const { return c_N2+c_O2+c_Ar+c_CO2+c_H20+c_Ne+c_He+c_O3; }
};
#endif

class AtmSlice
{
public:
  AtmSlice(double _v=0): zb(_v), zt(_v), tb(_v), tt(_v), rho(_v),
			 n_minus_one(_v)
#if 0
                       , temperature(_v), pressure(_v), composition()
#endif
  { }

  double zb;
  double zt;
  double tb;
  double tt;
  double rho;
  double n_minus_one;
#if 0
  double temperature;
  double pressure;
  AtmComposition composition;
#endif

  bool operator< (const AtmSlice& o) const { return zt < o.zt; }

#ifndef SWIG
  class CmpZAsc
  { public: bool operator() (const AtmSlice& a, const AtmSlice& b)
    { return a.zt < b.zt; } };
  /*
  class CmpTDec
  { public: bool operator() (const AtmSlice& a, const AtmSlice& b)
  { return b.tb < a.tb; } };*/
#endif
};

inline double sin2_thetac_for_gamma_sq(double g2, double n)
{
  const double b2 = 1.0 - 1.0/g2;                    // beta^2
  return 1.0 - 1.0/(b2*calin::math::special::SQR(n));
}

inline double sin2_thetac_for_energy(double energy, double mass, double n)
{
  double g2 = calin::math::special::SQR(energy/mass); // gamma^2
  return sin2_thetac_for_gamma_sq(g2, n);
}

// Base class for all atmospheric models
class Atmosphere
{
 public:
  virtual ~Atmosphere();
  virtual double rho(double z) = 0;
  virtual double thickness(double z) = 0;
  virtual double n_minus_one(double z) = 0;
  virtual double dn_dz(double z, double& n_minus_one) = 0;
  virtual void cherenkov_parameters(double z,
    double& n_minus_one, double& propagation_ct_correction);
#if 0
  virtual double pressure(double z) = 0;
  virtual double temperature(double z) = 0;
  virtual AtmComposition composition(double z) = 0;
#endif
  virtual double propagation_ct_correction(double z) = 0;
  virtual double z_for_thickness(double t) = 0;
  virtual double top_of_atmosphere() = 0;
  std::vector<AtmSlice> make_atm_slices(unsigned nslice,
                                        double zmax, double zmin);
  std::vector<AtmSlice> make_atm_slices(unsigned nslice);

  double sin2_thetac_for_gamma_sq(double g2, double z) {
    return calin::simulation::atmosphere::
      sin2_thetac_for_gamma_sq(g2, 1.0+this->n_minus_one(z)); }
  double sin2_thetac_for_energy(double energy, double mass, double z) {
    return calin::simulation::atmosphere::
      sin2_thetac_for_energy(energy, mass, 1.0+this->n_minus_one(z)); }
};

// Simple (i.e. unrealistic) isothermal atmosphere
class IsothermalAtmosphere: public Atmosphere
{
 public:
  IsothermalAtmosphere(double rho0=1.2e-3, double zs = 8.5e5,
		       double zmax = 1.2e7, double nmo0=2.75e-4,
                       double temperature = 273.);
  virtual ~IsothermalAtmosphere();
  double rho(double z) override;
  double thickness(double z) override;
  double n_minus_one(double z) override;
  double dn_dz(double z, double& n_minus_one) override;
#if 0
  double pressure(double z) override;
  double temperature(double z) override;
  AtmComposition composition(double z) override;
#endif
  double propagation_ct_correction(double z) override;
  double z_for_thickness(double t) override;
  double top_of_atmosphere() override;
private:
  double m_ttoa;
  double m_rho0;
  double m_zs;
  double m_zmax;
  double m_nmo0;
  //double temperature_;
};

struct LayeredAtmosphereLevel
{
  double z;
  double rho;
  double t;
  double nmo;

#ifndef SWIG
  class CmpZAsc { public: bool operator() (const LayeredAtmosphereLevel& a,
                                           const LayeredAtmosphereLevel& b) {
    return a.z < b.z; } };

  class CmpTDec { public: bool operator() (const LayeredAtmosphereLevel& a,
                                           const LayeredAtmosphereLevel& b) {
    return b.t < a.t; } };
#endif // ifndef SWIG
};

#ifndef SWIG
struct LayeredAtmosphereLayer
{
  LayeredAtmosphereLayer(double _v = 0):
      zb(_v), zt(_v), rho0(_v), rhozs(_v), t0(_v), tzs(_v), tb(_v), tt(_v),
      nmo0(_v), nmozs(_v), ptc0(_v) { }
  double zb;
  double zt;
  double rho0;
  double rhozs;
  double t0;
  double tzs;
  double tb;
  double tt;
  double nmo0;
  double nmozs;
  double ptc0;
  bool operator< (const LayeredAtmosphereLayer& o) const { return zt<o.zt; }

  class CmpTDec { public: bool operator() (const LayeredAtmosphereLayer& a,
                                           const LayeredAtmosphereLayer& b) {
    return b.tt<a.tt; } };
};
#endif // ifndef SWIG

std::vector<LayeredAtmosphereLevel> us76_levels();
std::vector<LayeredAtmosphereLevel> load_levels(const std::string& filename);

// Parameterized, layered model
class LayeredAtmosphere: public Atmosphere
{
public:
  CALIN_TYPEALIAS(Level, LayeredAtmosphereLevel);

  LayeredAtmosphere(const std::string& filename);
  LayeredAtmosphere(const std::vector<Level> levels);

  virtual ~LayeredAtmosphere();
  double rho(double z) override;
  double thickness(double z) override;
  double n_minus_one(double z) override;
  double dn_dz(double z, double& n_minus_one) override;
#if 0
  double pressure(double z) override;
  double temperature(double z) override;
  AtmComposition composition(double z) override;
#endif
  double propagation_ct_correction(double z) override;
  void cherenkov_parameters(double z,
    double& n_minus_one, double& propagation_ct_correction) override;
  double z_for_thickness(double t) override;
  double top_of_atmosphere() override;

  const std::vector<Level>& getLevels() const { return m_levels; }

  static LayeredAtmosphere* us76();

private:
  CALIN_TYPEALIAS(Layer, LayeredAtmosphereLayer);

  void initialize();
  inline std::vector<Layer>::const_iterator findZ(double z) const;

  double m_ztoa;
  double m_ttoa;
  std::vector<Level> m_levels;
  std::vector<Layer> m_layers;
  mutable std::vector<Layer>::const_iterator m_ilayer;
};

#if 1
// Parameterized, layered model that is capable of calculating refraction
// corrections in impact time and point
class LayeredRefractiveAtmosphere: public Atmosphere
{
public:
  CALIN_TYPEALIAS(Level, LayeredAtmosphereLevel);

  LayeredRefractiveAtmosphere(const std::string& filename,
    const std::vector<double>& obs_levels = {},
    const calin::ix::simulation::atmosphere::LayeredRefractiveAtmosphereConfig& config = default_config());
  LayeredRefractiveAtmosphere(const std::string& filename,
    double obs_level,
    const calin::ix::simulation::atmosphere::LayeredRefractiveAtmosphereConfig& config = default_config());
  LayeredRefractiveAtmosphere(const std::vector<Level>& levels,
    const std::vector<double>& obs_levels = {},
    const calin::ix::simulation::atmosphere::LayeredRefractiveAtmosphereConfig& config = default_config());

  virtual ~LayeredRefractiveAtmosphere();

  double num_obs_levels() const { return zobs_.size(); }
  double zobs(unsigned iground) const { return zobs_[iground]; }

  double rho(double z) override;
  double thickness(double z) override;
  double n_minus_one(double z) override;
  double dn_dz(double z, double& n_minus_one) override;

  double propagation_ct_correction(double z) override;
  void cherenkov_parameters(double z,
    double& n_minus_one, double& propagation_ct_correction) override;

  double z_for_thickness(double t) override;
  double top_of_atmosphere() override;

  bool propagate_ray_with_refraction(calin::math::ray::Ray& ray, unsigned iobs=0,
    bool time_reversal_ok=true);

  const std::vector<Level>& get_levels() const { return levels_; }

  static LayeredRefractiveAtmosphere* us76(const std::vector<double>& obs_levels = {});

  const calin::math::spline_interpolation::CubicMultiSpline* spline() const { return s_; }

  const Eigen::VectorXd& test_ray_emi_zn() { return test_ray_emi_zn_; }
  const Eigen::VectorXd& test_ray_emi_z() { return test_ray_emi_z_; }
  const Eigen::MatrixXd& test_ray_emi_x() { return test_ray_emi_x_; }
  const Eigen::MatrixXd& test_ray_emi_ct() { return test_ray_emi_ct_; }

  const Eigen::MatrixXd& test_ray_obs_x(unsigned iobs) { return test_ray_obs_x_[iobs]; }
  const Eigen::MatrixXd& test_ray_obs_ct(unsigned iobs) { return test_ray_obs_ct_[iobs]; }

  const Eigen::MatrixXd& test_ray_boa_x() { return test_ray_boa_x_; }
  const Eigen::MatrixXd& test_ray_boa_ct() { return test_ray_boa_ct_; }

  const Eigen::MatrixXd& model_ray_obs_x(unsigned iobs) { return model_ray_obs_x_[iobs]; }
  const Eigen::MatrixXd& model_ray_obs_ct(unsigned iobs) { return model_ray_obs_ct_[iobs]; }

  static calin::ix::simulation::atmosphere::LayeredRefractiveAtmosphereConfig default_config();

private:
  void initialize();

  std::vector<Level> levels_;
  calin::math::spline_interpolation::CubicMultiSpline* s_ = nullptr;

  std::vector<double> zobs_;

  struct ObsLevelData {
    double z = 0.0;
    double n_inv = 1.0;
    double x_ba_ratio = 0.0;
    double ct_ba_ratio = 0.0;
    double ct_x_ratio = 0.0;
  };

  std::vector<ObsLevelData> obs_level_data_;
  bool high_accuracy_mode_ = false;

  Eigen::VectorXd test_ray_emi_zn_;
  Eigen::VectorXd test_ray_emi_z_;
  Eigen::MatrixXd test_ray_emi_x_;
  Eigen::MatrixXd test_ray_emi_ct_;
  std::vector<Eigen::MatrixXd> test_ray_obs_x_;
  std::vector<Eigen::MatrixXd> test_ray_obs_ct_;
  Eigen::MatrixXd test_ray_boa_x_;
  Eigen::MatrixXd test_ray_boa_ct_;
  std::vector<Eigen::MatrixXd> model_ray_obs_x_;
  std::vector<Eigen::MatrixXd> model_ray_obs_ct_;
};
#endif
} } } // namespace calin::simulation::atmosphere
