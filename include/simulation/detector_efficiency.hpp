/*

   calin/simulation/detector_efficiency.hpp -- Stephen Fegan -- 2016-10-14

   Classes to do 1D linear or exponential interpolation over detector efficiency
   and Cherenkov bandwidth tables.

   Originally from EGS5 ACT code, see header below.

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

// DetectorEfficiency.hpp - Classes to handle efficiency of detector
// - inculding absoption of Cherenkov light along the path
// Stephen Fegan - sfegan@llr.in2p3.fr - November 2012
// $Id: DetectorEfficiency.hpp 5422 2013-06-26 14:01:03Z sfegan $

#pragma once

// #define ACT_LIGHT_YIELD_TAYLOR_SERIES_IN_LOG

#include <cmath>
#include <tuple>

#include <math/interpolation_1d.hpp>

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
};

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
  DetectionEfficiency();
  void scaleEff(const calin::math::interpolation_1d::InterpLinear1D& eff);
  void scaleEffFromFile(const std::string& filename);
  void scaleEffFromOldStyleFile(const std::string& filename,
			double lambda0_nm=180.0, double dlambda_nm=5.0);
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
  AtmosphericAbsorption(const std::string& filename, OldStyleAtmObsFlag flag,
    double ground_level_km = 0.0, double spacing_km=1.0);
  AtmosphericAbsorption(const std::string& filename,
    std::vector<double> levels_cm = {});
  calin::math::interpolation_1d::InterpLinear1D absorptionForAltitude(double h) const;
  ACTEffectiveBandwidth integrateBandwidth(double h0, double w0,
    const DetectionEfficiency& eff) const;
  const std::vector<double>& energy_ev() const { return e_ev_; }
  std::vector<double> levels_cm() const;
private:
  std::vector<double>                                        e_ev_;
  std::vector<calin::math::interpolation_1d::InterpLinear1D> absorption_;
};

} } } // namespace calin::simulation::detector_efficiency
