/*

   calin/simulation/detector_efficiency.hpp -- Stephen Fegan -- 2016-10-14

   Classes to do 1D linear or exponential interpolation over detector efficiency
   and Cherenkov yield tables.

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

struct CherenkovYieldTaylorCoefficients
{
  CherenkovYieldTaylorCoefficients() { /* nothing to see here */ }
  CherenkovYieldTaylorCoefficients(double n_, double dn_dw_, double d2n_dw2_):
    n(n_), dn_dw(dn_dw_), d2n_dw2(d2n_dw2_) { /* nothing to see here */ }

  double n         = 0;
  double dn_dw     = 0;
  double d2n_dw2   = 0;

  CherenkovYieldTaylorCoefficients&
  operator+= (const CherenkovYieldTaylorCoefficients o) {
    n         += o.n;
    dn_dw     += o.dn_dw;
    d2n_dw2   += o.d2n_dw2;
    return *this;
  }

  CherenkovYieldTaylorCoefficients&
  operator-= (const CherenkovYieldTaylorCoefficients o) {
    n         -= o.n;
    dn_dw     -= o.dn_dw;
    d2n_dw2   -= o.d2n_dw2;
    return *this;
  }

  CherenkovYieldTaylorCoefficients& operator*= (double c) {
    n         *= c;
    dn_dw     *= c;
    d2n_dw2   *= c;
    return *this;
  }

  CherenkovYieldTaylorCoefficients
  operator+ (const CherenkovYieldTaylorCoefficients o) const {
    CherenkovYieldTaylorCoefficients r(*this);
    r += o;
    return r;
  }

  CherenkovYieldTaylorCoefficients
  operator- (const CherenkovYieldTaylorCoefficients o) const {
    CherenkovYieldTaylorCoefficients r(*this);
    r -= o;
    return r;
  }

  CherenkovYieldTaylorCoefficients operator* (double c) const {
    CherenkovYieldTaylorCoefficients r(*this);
    r *= c;
    return r;
  }

  bool operator< (const CherenkovYieldTaylorCoefficients o) const {
    return false; // All coefficients equal in sort
  }
};

typedef CherenkovYieldTaylorCoefficients yield_t;

#if 0
inline yield_t operator/ (const yield_t& a, const yield_t& b)
{
  return yield_t(a.first/b.first, a.second/b.second, a.third/b.third);
}

inline yield_t operator* (const yield_t& a, const yield_t& b)
{
  return yield_t(a.first*b.first, a.second*b.second, a.third*b.third);
}
#endif

#ifdef SWIG
%template(InterpLinear1DDouble) Interpolation1D<double, LinearInterpolator<double> >;
%template(InterpLinear1DYieldT) Interpolation1D<yield_t, LinearInterpolator<yield_t> >;
#endif

class TelescopeEfficiency: public math::interpolation_1d::InterpLinear1D
{
public:
  TelescopeEfficiency();
  void scaleEff(const math::interpolation_1d::InterpLinear1D& eff);
  void scaleEffFromFile(const std::string& filename,
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
  inline yield_t interpolate(double x,
			     double x0, const yield_t& y0,
			     double x1, const yield_t& y1) const
  {
    yield_t y;
    y.first  = EXP::interpolate(x, x0, y0.first,  x1, y1.first);
    y.second = LIN::interpolate(x, x0, y0.second, x1, y1.second);
    y.third  = LIN::interpolate(x, x0, y0.third,  x1, y1.third);
    return y;
  }

  inline yield_t integrate(double x0, const yield_t& y0,
			   double x1, const yield_t& y1) const
  {
    yield_t y;
    y.first  = EXP::integrate(x0, y0.first,  x1, y1.first);
    y.second = LIN::integrate(x0, y0.second, x1, y1.second);
    y.third  = LIN::integrate(x0, y0.third,  x1, y1.third);
    return y;
  }
};

#ifdef SWIG
%template(Interp1DACT) Interpolation1D<yield_t, ACTILYInterpolator>;
#endif
#endif

class ACTIntegratedLightYield:
  public math::interpolation_1d::Interpolation1D<yield_t, math::interpolation_1d::LinearInterpolator<yield_t> >
  //public Interpolation1D<yield_t, ExpInterpolator<yield_t> >
  //public Interpolation1D<yield_t, ACTILYInterpolator>
{
public:
  ACTIntegratedLightYield(double w0=0);
  double yield(double h, double w) const;
  double w0() const { return m_w0; }
private:
  double m_w0;
};

class AtmosphericAbsorption
{
public:
  AtmosphericAbsorption(const std::string& filename, double spacing_km=1.0);
  math::interpolation_1d::InterpLinear1D absorptionForAltitude(double h) const;
  ACTIntegratedLightYield integrateYield(double h0, double w0,
					 const TelescopeEfficiency& eff);
private:
  std::vector<double>         m_e_ev;
  std::vector<math::interpolation_1d::InterpLinear1D> m_absorption;
};

} } } // namespace calin::simulation::detector_efficiency
