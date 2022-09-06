/*

   calin/simulation/detector_efficiency.hpp -- Stephen Fegan -- 2016-10-14

   Classes to do 1D linear or exponential interpolation over detector efficiency
   and Cherenkov bandwidth tables.

   Originally from EGS5 ACT code, see header below.

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

// DetectorEfficiency.hpp - Classes to handle efficiency of detector
// - inculding absoption of Cherenkov light along the path
// Stephen Fegan - sfegan@llr.in2p3.fr - November 2012
// $Id: DetectorEfficiency.hpp 5422 2013-06-26 14:01:03Z sfegan $

#pragma once

// #define ACT_LIGHT_YIELD_TAYLOR_SERIES_IN_LOG

#include <cmath>
#include <tuple>
#include <ostream>

#include <util/vcl.hpp>
#include <math/rng.hpp>
#include <math/rng_vcl.hpp>
#include <math/interpolation_1d.hpp>
#include <math/spline_interpolation.hpp>

namespace calin { namespace simulation { namespace detector_efficiency {

struct CherenkovBandwidthTaylorCoefficients
{
  CherenkovBandwidthTaylorCoefficients() { /* nothing to see here */ }
  CherenkovBandwidthTaylorCoefficients(double n_, double dn_dw_, double d2n_dw2_):
    n(n_), dn_dw(dn_dw_), d2n_dw2(d2n_dw2_) { /* nothing to see here */ }

  double n         = 0;
  double dn_dw     = 0;
  double d2n_dw2   = 0;

  CherenkovBandwidthTaylorCoefficients&
  operator+= (const CherenkovBandwidthTaylorCoefficients o) {
    n         += o.n;
    dn_dw     += o.dn_dw;
    d2n_dw2   += o.d2n_dw2;
    return *this;
  }

  CherenkovBandwidthTaylorCoefficients&
  operator-= (const CherenkovBandwidthTaylorCoefficients o) {
    n         -= o.n;
    dn_dw     -= o.dn_dw;
    d2n_dw2   -= o.d2n_dw2;
    return *this;
  }

  CherenkovBandwidthTaylorCoefficients&
  operator*= (const CherenkovBandwidthTaylorCoefficients o) {
    n         *= o.n;
    dn_dw     *= o.dn_dw;
    d2n_dw2   *= o.d2n_dw2;
    return *this;
  }

  CherenkovBandwidthTaylorCoefficients&
  operator/= (const CherenkovBandwidthTaylorCoefficients o) {
    n         /= o.n;
    dn_dw     /= o.dn_dw;
    d2n_dw2   /= o.d2n_dw2;
    return *this;
  }

  CherenkovBandwidthTaylorCoefficients& operator*= (double c) {
    n         *= c;
    dn_dw     *= c;
    d2n_dw2   *= c;
    return *this;
  }

  CherenkovBandwidthTaylorCoefficients
  operator+ (const CherenkovBandwidthTaylorCoefficients o) const {
    CherenkovBandwidthTaylorCoefficients r(*this);
    r += o;
    return r;
  }

  CherenkovBandwidthTaylorCoefficients
  operator- (const CherenkovBandwidthTaylorCoefficients o) const {
    CherenkovBandwidthTaylorCoefficients r(*this);
    r -= o;
    return r;
  }

  CherenkovBandwidthTaylorCoefficients operator* (double c) const {
    CherenkovBandwidthTaylorCoefficients r(*this);
    r *= c;
    return r;
  }

  bool operator< (const CherenkovBandwidthTaylorCoefficients o) const {
    return false; // All coefficients equal in sort
  }

  bool operator== (const CherenkovBandwidthTaylorCoefficients o) const {
    return n==o.n and dn_dw==o.dn_dw and d2n_dw2==o.d2n_dw2;
  }
};

inline std::ostream& operator<<(std::ostream& stream, const CherenkovBandwidthTaylorCoefficients& c)
{
  stream << c.n << ' ' << c.dn_dw << ' ' << c.d2n_dw2;
  return stream;
}

inline std::istream& operator>>(std::istream& stream, CherenkovBandwidthTaylorCoefficients& c)
{
  stream >> c.n >> c.dn_dw >> c.d2n_dw2;
  return stream;
}

typedef CherenkovBandwidthTaylorCoefficients bandwidth_t;

inline bandwidth_t operator/ (bandwidth_t a, const bandwidth_t& b)
{
  a /= b;
  return a;
}

inline bandwidth_t operator* (bandwidth_t a, const bandwidth_t& b)
{
  a *= b;
  return a;
}

class DetectionEfficiency: public calin::math::interpolation_1d::InterpLinear1D
{
public:
  DetectionEfficiency(double const_eff = 1.0);
  DetectionEfficiency(const std::string& filename);
  void scaleEff(const calin::math::interpolation_1d::InterpLinear1D& eff);
  void scaleEffByConst(double c);
  void scaleEffFromFile(const std::string& filename);
  void scaleEffFromOldStyleFile(const std::string& filename,
		double lambda0_nm=180.0, double dlambda_nm=5.0);
  bool operator==(const DetectionEfficiency& o) const {
    return static_cast<const calin::math::interpolation_1d::InterpLinear1D&>(*this) ==
      static_cast<const calin::math::interpolation_1d::InterpLinear1D&>(o);
  }
};

class AngularEfficiency: public calin::math::interpolation_1d::InterpLinear1D
{
public:
  AngularEfficiency(double const_eff = 1.0);
  AngularEfficiency(const std::string& filename);
  void scaleEff(const calin::math::interpolation_1d::InterpLinear1D& eff);
  void scaleEffByConst(double c);
  void scaleEffFromFile(const std::string& filename);
  bool operator==(const AngularEfficiency& o) const {
    return static_cast<const calin::math::interpolation_1d::InterpLinear1D&>(*this) ==
      static_cast<const calin::math::interpolation_1d::InterpLinear1D&>(o);
  }
};

#if 0
class ACTILYInterpolator
  : private ExpInterpolator<double>, private LinearInterpolator<double>
{
private:
  typedef ExpInterpolator<double> EXP;
  typedef LinearInterpolator<double> LIN;
public:
  inline bandwidth_t interpolate(double x,
			     double x0, const bandwidth_t& y0,
			     double x1, const bandwidth_t& y1) const
  {
    bandwidth_t y;
    y.first  = EXP::interpolate(x, x0, y0.first,  x1, y1.first);
    y.second = LIN::interpolate(x, x0, y0.second, x1, y1.second);
    y.third  = LIN::interpolate(x, x0, y0.third,  x1, y1.third);
    return y;
  }

  inline bandwidth_t integrate(double x0, const bandwidth_t& y0,
			   double x1, const bandwidth_t& y1) const
  {
    bandwidth_t y;
    y.first  = EXP::integrate(x0, y0.first,  x1, y1.first);
    y.second = LIN::integrate(x0, y0.second, x1, y1.second);
    y.third  = LIN::integrate(x0, y0.third,  x1, y1.third);
    return y;
  }
};

#ifdef SWIG
%template(Interp1DACT) Interpolation1D<bandwidth_t, ACTILYInterpolator>;
#endif
#endif

#ifdef SWIG
} } }
%template(LinearInterpolatorBandwidthT)
  calin::math::interpolation_1d::LinearInterpolator<
    calin::simulation::detector_efficiency::bandwidth_t>;
%template(InterpLinear1DBandwidthT)
  calin::math::interpolation_1d::Interpolation1D<
    calin::simulation::detector_efficiency::bandwidth_t,
    calin::math::interpolation_1d::LinearInterpolator<
      calin::simulation::detector_efficiency::bandwidth_t> >;
namespace calin { namespace simulation { namespace detector_efficiency {
#endif

class ACTEffectiveBandwidth:
  public calin::math::interpolation_1d::Interpolation1D<bandwidth_t, calin::math::interpolation_1d::LinearInterpolator<bandwidth_t> >
  //public Interpolation1D<bandwidth_t, ExpInterpolator<bandwidth_t> >
  //public Interpolation1D<bandwidth_t, ACTILYInterpolator>
{
public:
  ACTEffectiveBandwidth(double w0=0);
  double bandwidth(double h, double w) const;
  double w0() const { return m_w0; }
private:
  double m_w0;
};

struct OldStyleAtmObsFlag { };

class AtmosphericAbsorption
{
public:
  AtmosphericAbsorption();
  AtmosphericAbsorption(const std::string& filename, OldStyleAtmObsFlag flag,
    double ground_level_km = 0.0, double spacing_km=1.0);
  AtmosphericAbsorption(const std::string& filename,
    std::vector<double> levels_cm = {});
  calin::math::interpolation_1d::InterpLinear1D optical_depth_for_altitude(double h) const;
  double optical_depth_for_altitude_and_energy(double h, double e) const;
  ACTEffectiveBandwidth integrateBandwidth(double h0, double w0,
    const DetectionEfficiency& eff) const;
  ACTEffectiveBandwidth integrateBandwidth(double h0, double w0,
    const DetectionEfficiency& eff, double emin, double emax) const;
  calin::math::spline_interpolation::TwoDimensionalCubicSpline* integrate_bandwidth_to_spline(
    double h0, const DetectionEfficiency& eff, std::vector<double> h = { },
    std::vector<double> w = { 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0 }, double emin=0, double emax=0) const;
  const std::vector<double>& energy_ev() const { return e_ev_; }
  std::vector<double> levels_cm() const;
  void set_zref(double zref);
  const calin::math::interpolation_1d::InterpLinear1D& absorption(unsigned ie) {
    return absorption_.at(ie);
  }
private:
  std::vector<double>                                        e_ev_;
  std::vector<calin::math::interpolation_1d::InterpLinear1D> absorption_;
};

class PEAmplitudeGenerator
{
public:
  virtual ~PEAmplitudeGenerator();
  virtual double generate_amplitude() = 0;
  Eigen::VectorXd bulk_generate_amplitude(unsigned n) {
    Eigen::VectorXd rvs(n);
    for(unsigned i=0;i<n;i++) {
      rvs(i) = this->generate_amplitude();
    }
    return rvs;
  }
};

enum SplineMode {
  SM_LINEAR,
  SM_LOG,
  SM_SQRT_LOG
};

class SplinePEAmplitudeGenerator: public PEAmplitudeGenerator
{
public:
  SplinePEAmplitudeGenerator(const Eigen::VectorXd& q, const Eigen::VectorXd& dp_dq,
    SplineMode spline_mode, calin::math::rng::RNG* rng = nullptr, bool adopt_rng = false);
  SplinePEAmplitudeGenerator(const calin::math::spline_interpolation::CubicSpline& spline,
    SplineMode spline_mode, calin::math::rng::RNG* rng = nullptr, bool adopt_rng = false);
  virtual ~SplinePEAmplitudeGenerator();
  double generate_amplitude() final;
  template<typename VCLArchitecture> typename VCLArchitecture::double_vt vcl_generate_amplitude(
    calin::math::rng::VCLRNG<VCLArchitecture>& rng)
  {
    typename VCLArchitecture::double_vt x = rng->uniform_double();
    switch(spline_mode_) {
    case SM_LINEAR:
      break;
    case SM_LOG:
      x = -vcl::log(x);
      break;
    case SM_SQRT_LOG:
      x = vcl::sqrt(-vcl::log(x));
      break;
    }
    return spline_->vcl_value(x);
  }
  const calin::math::spline_interpolation::CubicSpline& spline() const { return *spline_; }
  SplineMode spline_mode() const { return spline_mode_; }
  static calin::math::spline_interpolation::CubicSpline* make_spline(
    const Eigen::VectorXd& q, const Eigen::VectorXd& dp_dq, SplineMode spline_mode,
    bool regularize_spline = true, bool extend_linear_rhs = true,
    unsigned regularize_ninterval = 0,
    double norm = 1.0-std::numeric_limits<double>::epsilon());
protected:
  calin::math::spline_interpolation::CubicSpline* spline_ = nullptr;
  SplineMode spline_mode_ = SM_LINEAR;
  calin::math::rng::RNG* rng_ = nullptr;
  bool adopt_rng_ = false;
};

} } } // namespace calin::simulation::detector_efficiency
