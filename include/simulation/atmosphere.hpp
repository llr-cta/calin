/*

   calin/simulation/atmosphere.hpp -- Stephen Fegan -- 2015-06-11

   Classes to model density and refractive index of atmosphere

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

// Originally: Atmposphere.hpp - Classes to handle atmosphere
// Stephen Fegan - sfegan@llr.in2p3.fr - July 2012
// Id: Atmosphere.hpp 5416 2013-06-24 13:46:00Z sfegan

#pragma once

#include<vector>
#include<cassert>

#include"calin_global_definitions.hpp"

#include <math/special.hpp>
#include <math/ray.hpp>
#include <util/vcl.hpp>
#include <math/ray_vcl.hpp>
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
  static double solve_for_thickness_at_toa(const std::vector<Level>& levels);

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

// Parameterized, layered model that is capable of calculating refraction
// corrections in impact time and point
// Also supports SIMD interface
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
  LayeredRefractiveAtmosphere(const LayeredRefractiveAtmosphere& o) : Atmosphere(),
    levels_(o.levels_), s_(new calin::math::spline_interpolation::CubicMultiSpline(*o.s_)),
    zobs_(o.zobs_), thickness_to_toa_(o.thickness_to_toa_), obs_level_data_(o.obs_level_data_),
    high_accuracy_mode_(o.high_accuracy_mode_), test_ray_emi_zn_(o.test_ray_emi_zn_),
    test_ray_emi_z_(o.test_ray_emi_z_), test_ray_emi_x_(o.test_ray_emi_x_),
    test_ray_emi_ct_(o.test_ray_emi_ct_), test_ray_obs_x_(o.test_ray_obs_x_),
    test_ray_obs_ct_(o.test_ray_obs_ct_), test_ray_boa_x_(o.test_ray_boa_x_),
    test_ray_boa_ct_(o.test_ray_boa_ct_), model_ray_obs_x_(o.model_ray_obs_x_),
    model_ray_obs_ct_(o.model_ray_obs_ct_) { }

  virtual ~LayeredRefractiveAtmosphere();

  LayeredRefractiveAtmosphere* clone() const { return new LayeredRefractiveAtmosphere(*this); }

  unsigned num_obs_levels() const { return zobs_.size(); }
  double zobs(unsigned iground) const { return zobs_[iground]; }

  double rho(double z) override;
  double thickness(double z) override;
  double n_minus_one(double z) override;
  double dn_dz(double z, double& n_minus_one) override;

  double propagation_ct_correction(double z) override;
  double propagation_ct_correction_to_iobs(double z, unsigned iobs);

  void cherenkov_parameters(double z,
    double& n_minus_one, double& propagation_ct_correction) override;

  double z_for_thickness(double t) override;
  double top_of_atmosphere() override;

  bool propagate_ray_with_refraction(calin::math::ray::Ray& ray, unsigned iobs=0,
    bool time_reversal_ok=true);

  double refraction_offset(double z, double theta_rad, unsigned iobs=0);
  double refraction_safety_radius(double zenith_rad, unsigned iobs=0);

  template<typename VCLArchitecture> inline typename VCLArchitecture::double_vt
  vcl_n_minus_one(typename VCLArchitecture::double_vt z) const
  {
    return vcl::exp(s_->vcl_value<VCLArchitecture>(z, 2));
  }

  template<typename VCLArchitecture> inline typename VCLArchitecture::double_vt
  vcl_dn_dz(typename VCLArchitecture::double_vt z,
    typename VCLArchitecture::double_vt& n_minus_one) const
  {
    typename VCLArchitecture::double_vt dlogn_dz =
      s_->vcl_derivative_and_value<VCLArchitecture>(z, 2, n_minus_one);
    n_minus_one = vcl::exp(n_minus_one);
    return n_minus_one * dlogn_dz;
  }

  template<typename VCLArchitecture> inline typename VCLArchitecture::double_vt
  vcl_dlognmo_dz(typename VCLArchitecture::double_vt z,
    typename VCLArchitecture::double_vt& n_minus_one) const
  {
    typename VCLArchitecture::double_vt dlognmo_dz =
      s_->vcl_derivative_and_value<VCLArchitecture>(z, 2, n_minus_one);
    n_minus_one = vcl::exp(n_minus_one);
    return dlognmo_dz;
  }

  template<typename VCLArchitecture> inline typename VCLArchitecture::double_vt
  vcl_propagation_ct_correction_to_iobs(typename VCLArchitecture::double_vt z, unsigned iobs=0) const
  {
    return s_->vcl_value<VCLArchitecture>(z, 4+iobs*5);
  }

  template<typename VCLArchitecture> typename VCLArchitecture::double_bvt
  vcl_propagate_ray_with_refraction_and_mask(
    calin::math::ray::VCLRay<typename VCLArchitecture::double_real>& ray,
    typename VCLArchitecture::double_bvt ray_mask,
    unsigned iobs=0, bool time_reversal_ok=true);

  const std::vector<Level>& get_levels() const { return levels_; }

  static LayeredRefractiveAtmosphere* us76(const std::vector<double>& obs_levels = {});

  const calin::math::spline_interpolation::CubicMultiSpline* spline() const { return s_; }
  void regularize_internal_spline(double dx = 0);

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

  bool test_vcl_propagate_ray_with_refraction_and_mask(
    calin::math::ray::Ray& ray, bool ray_mask, unsigned iobs=0, bool time_reversal_ok=true);

  static calin::ix::simulation::atmosphere::LayeredRefractiveAtmosphereConfig default_config();

private:
  std::vector<Level> levels_;
  calin::math::spline_interpolation::CubicMultiSpline* s_ = nullptr;

  std::vector<double> zobs_;
  double thickness_to_toa_ = 0;

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

template<typename R>
inline void calculate_refraction_angular_terms(
  R sin_i, R sec_i, R n_i_over_n_r,
  R& sin_r, R& cos_r, R& t1_x, R& t2_x, R& t1_ct, R& t2_ct)
{
  using std::sqrt;
  using vcl::sqrt;
  sin_r = sin_i * n_i_over_n_r;
  const R cos2_r = 1.0 - calin::math::special::SQR(sin_r);
  cos_r = sqrt(cos2_r);
  const R sec2_i = calin::math::special::SQR(sec_i);

  t1_x = sin_i * sec2_i / cos_r;
  t2_x = t1_x * calin::math::special::SQR(sin_i) * sec2_i;

  t1_ct = t1_x * sin_i;
  t2_ct = t2_x * sin_i;
}

template<typename VCLArchitecture> typename VCLArchitecture::double_bvt
LayeredRefractiveAtmosphere::vcl_propagate_ray_with_refraction_and_mask(
  calin::math::ray::VCLRay<typename VCLArchitecture::double_real>& ray,
  typename VCLArchitecture::double_bvt ray_mask,
  unsigned iobs, bool time_reversal_ok)
{
  typename VCLArchitecture::double_vt dz = ray.z() - obs_level_data_[iobs].z;

  const typename VCLArchitecture::double_vt cos_i = -ray.uz();

  ray_mask &= cos_i>0;
  if(not time_reversal_ok) {
    ray_mask &= dz>0;
  }

  if(not vcl::horizontal_or(ray_mask)) {
    // We outie ...
    return ray_mask;
  }

  typename VCLArchitecture::double_vt n_i;
  typename VCLArchitecture::double_vt v_ct;
  typename VCLArchitecture::double_vt a_x;
  typename VCLArchitecture::double_vt b_x;
  typename VCLArchitecture::double_vt a_ct;
  typename VCLArchitecture::double_vt b_ct;

  unsigned ispline = 4+iobs*5;
  if(high_accuracy_mode_) {
    // 6 spline interpolations for all required data
    s_->vcl_value<VCLArchitecture>(ray.z(),
      /* 1 */ 2, n_i,
      /* 2 */ ispline, v_ct,
      /* 3 */ ispline+1, a_x,
      /* 4 */ ispline+2, b_x,
      /* 5 */ ispline+3, a_ct,
      /* 6 */ ispline+4, b_ct);
  } else {
    // 3 spline interpolations for n, vertical ct and horizontal x corrections
    s_->vcl_value<VCLArchitecture>(ray.z(),
      /* 1 */ 2, n_i,
      /* 2 */ ispline, v_ct,
      /* 3 */ ispline+1, a_x);
    b_x = obs_level_data_[iobs].x_ba_ratio * a_x;
    a_ct = obs_level_data_[iobs].ct_x_ratio * a_x;
    b_ct = obs_level_data_[iobs].ct_ba_ratio * a_ct;
  }
  n_i = 1.0 + vcl::exp(n_i);

  const typename VCLArchitecture::double_vt n_r_inv = obs_level_data_[iobs].n_inv;
  const typename VCLArchitecture::double_vt n_i_over_n_r = n_i * n_r_inv;

  const typename VCLArchitecture::double_vt sec_i = 1.0/cos_i;
  const typename VCLArchitecture::double_vt sin2_i = calin::math::special::SQR(ray.ux()) + calin::math::special::SQR(ray.uy());
  const typename VCLArchitecture::double_vt sin_i = vcl::sqrt(sin2_i);

  typename VCLArchitecture::double_vt sin_r;
  typename VCLArchitecture::double_vt cos_r;
  typename VCLArchitecture::double_vt t1_x;
  typename VCLArchitecture::double_vt t2_x;
  typename VCLArchitecture::double_vt t1_ct;
  typename VCLArchitecture::double_vt t2_ct;

  calculate_refraction_angular_terms(sin_i, sec_i, n_i_over_n_r,
    sin_r, cos_r, t1_x, t2_x, t1_ct, t2_ct);

  ray.propagate_dist_with_mask(ray_mask, dz*sec_i);

  // Bend direction for Snell's law
  ray.mutable_uz() = vcl::select(ray_mask, -cos_r, ray.uz());
  ray.mutable_ux() *= vcl::select(ray_mask, n_i_over_n_r, 1.0);
  ray.mutable_uy() *= vcl::select(ray_mask, n_i_over_n_r, 1.0);

  typename VCLArchitecture::double_vt dx = a_x * t1_x + b_x * t2_x;
  typename VCLArchitecture::double_vt dx_csc_r = vcl::select(sin_r == 0, 0, dx / vcl::abs(sin_r));

  ray.mutable_x() -= vcl::select(ray_mask, dx_csc_r * ray.ux(), 0);
  ray.mutable_y() -= vcl::select(ray_mask, dx_csc_r * ray.uy(), 0);

  ray.mutable_ct() += vcl::select(ray_mask, v_ct * sec_i - (a_ct * t1_ct + b_ct * t2_ct), 0);

#if 0
  std::cout << v_ct << ' ' << a_x << ' ' << b_x << ' ' << a_ct << ' ' << b_ct << '\n';
  std::cout << t1_x << ' ' << t2_x << ' ' << t1_ct << ' ' << t2_ct << '\n';
  std::cout << sin_i << ' ' << cos_i << ' ' << calin::math::special::SQR(sin_i)+calin::math::special::SQR(cos_i) << ' '
    << sin_r << ' ' << cos_r << ' ' << calin::math::special::SQR(sin_r)+calin::math::special::SQR(cos_r) << '\n';
#endif

  return ray_mask;
}

// Constant density "Atmosphere" - for doing some simple simulations in fixed
// substances (such as Iron etc)

// Units:
// Height, z:    cm
// Density, rho: g/cm^3
// Thickness, t: g/cm^2

class IsotropicAtmosphere: public Atmosphere
{
public:
  IsotropicAtmosphere(): Atmosphere() { /* nothing to see here */ }
  IsotropicAtmosphere(double rho, double ztop, double zground, double nmo):
    Atmosphere(),
    rho_(rho), zground_(zground), ztop_(ztop), n_minus_one_(nmo) { }
  virtual ~IsotropicAtmosphere();
  virtual double rho(double z) override;
  virtual double thickness(double z) override;
  virtual double n_minus_one(double z) override;
  virtual double dn_dz(double z, double& n_minus_one) override;

  virtual double propagation_ct_correction(double z) override;
  virtual double z_for_thickness(double t) override;
  virtual double top_of_atmosphere() override;

private:
  double rho_         = 0.12219E-02; // g/cm^3
  double zground_     = 0;           // cm;
  double ztop_        = 100e5;       // cm;
  double n_minus_one_ = 0.28232E-03;
};

} } } // namespace calin::simulation::atmosphere
