/*

   calin/simulation/detector_efficiency.cpp -- Stephen Fegan -- 2016-10-15

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
// $Id: DetectorEfficiency.cpp 5422 2013-06-26 14:01:03Z sfegan $

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include <io/log.hpp>
#include <math/special.hpp>
#include <simulation/detector_efficiency.hpp>

using calin::math::special::SQR;
using namespace calin::io::log;
using namespace calin::simulation::detector_efficiency;
using namespace calin::math::interpolation_1d;

#define EV_NM 1239.84193009239 // gunits: c/(ev/h) -> nm

// ----------------------------------------------------------------------------
// AtmosphericAbsorption
// ----------------------------------------------------------------------------

AtmosphericAbsorption::
AtmosphericAbsorption(const std::string& filename, double spacing_km)
{
  std::ifstream stream(filename.c_str());
  std::string line;
  std::getline(stream,line); // skip first line
  double e = 0;
  double z = 0;
  InterpLinear1D abs;
  std::getline(stream,line);
  while(stream)
  {
    if(line[1] != ' ')
    {
      std::istringstream lstream(line);
      if(e)
      {
        m_e_ev.push_back(e);
        m_absorption.push_back(abs);
      }
      double lambda;
      lstream >> lambda;
      e = EV_NM/lambda;
      z = 0;
      abs.clear();
    }
    else
    {
      std::istringstream lstream(line);
      double a;
      lstream >> a;
      while(lstream)
      {
        abs.insert(z*100000.0, a);
        z += spacing_km;
        lstream >> a;
      }
    }
    std::getline(stream,line);
  }
  if(e)
  {
    m_e_ev.push_back(e);
    m_absorption.push_back(abs);
  }
}

InterpLinear1D AtmosphericAbsorption::absorptionForAltitude(double h) const
{
  InterpLinear1D abs;
  for(unsigned ie=0;ie<m_e_ev.size();ie++)
    abs.insert(m_e_ev[ie], m_absorption[ie](h));
  return abs;
}

ACTIntegratedLightYield AtmosphericAbsorption::
integrateYield(double h0, double w0, const TelescopeEfficiency& eff)
{
  InterpLinear1D abs0 = absorptionForAltitude(h0);
  ACTIntegratedLightYield yield(w0);

#if 1
  bool obslevel = false;
  for(unsigned ih=0;ih<m_absorption.front().nXY();ih++)
  {
    double h = m_absorption.front().xi(ih);
    if(h==h0 && !obslevel)
    {
      obslevel = true;
    }
    else if(h>h0 && !obslevel)
    {
      ih--;
      h = h0;
      obslevel = true;
    }
#else
  for(double h=m_absorption.front().xmin();
      h<m_absorption.front().xmax(); h+=10000)
  {
#endif
    InterpLinear1D Y0;
    InterpLinear1D Y1;
    InterpLinear1D Y2;
    for(unsigned ie=0;ie<m_e_ev.size();ie++)
  	{
  	  double e = m_e_ev[ie];
  	  double mfp = std::fabs(m_absorption[ie](h)-abs0(e));
  	  double abs = std::exp(-mfp/w0);
  	  Y0.insert(e, abs);
  	  Y1.insert(e, mfp/(w0*w0)*abs);
  	  Y2.insert(e, mfp*(0.5*mfp/w0-1.0)/(w0*w0*w0)*abs);
  	}
    Y0 *= eff;
    Y1 *= eff;
    Y2 *= eff;

    yield_t y(Y0.integrate(), Y1.integrate(), Y2.integrate());
#ifdef ACT_LIGHT_YIELD_TAYLOR_SERIES_IN_LOG
    y.dn_dw   = y.dn_dw/y.n;
    y.d2n_dw2 = y.d2n_dw2/y.n - 0.5*SQR(y.dn_dw);
    y.n       = std::log(y.n);
#endif
    yield.insert(h, y);
  }
  return yield;
}

// ----------------------------------------------------------------------------
// TelescopeEfficiency
// ----------------------------------------------------------------------------

TelescopeEfficiency::TelescopeEfficiency(): InterpLinear1D(1.0)
{
  // nothing to see here
}

void TelescopeEfficiency::scaleEff(const InterpLinear1D& eff)
{
  *static_cast<InterpLinear1D*>(this) *= eff;
}

void TelescopeEfficiency::
scaleEffFromFile(const std::string& filename)
{
  std::ifstream stream(filename.c_str());
  std::string line;
  std::getline(stream,line);
  while((line.empty() or line[0] == '#') and stream)std::getline(stream,line);
  if(!stream)return;
  InterpLinear1D eff_fn;
  do
  {
    double lambda;
    double eff;
    std::istringstream lstream(line);
    lstream >> lambda >> eff;
    LOG(INFO) << lambda << ' ' << eff;
    if(lstream) {
      double e = EV_NM / lambda;
      eff_fn.insert(e, eff);
    }
    std::getline(stream,line);
  }while(stream);
  scaleEff(eff_fn);
}

void TelescopeEfficiency::
scaleEffFromOldStyleFile(const std::string& filename,
		 double lambda0_nm, double dlambda_nm)
{
  std::ifstream stream(filename.c_str());
  InterpLinear1D eff_fn;
  std::string line;
  std::getline(stream,line); // skip first line
  double lambda = lambda0_nm;
  double eff;
  stream >> eff;
  while(stream)
  {
    double e = EV_NM / lambda;
    eff_fn.insert(e, eff);
    lambda += dlambda_nm;
    stream >> eff;
  }
  scaleEff(eff_fn);
}

// ----------------------------------------------------------------------------
// ACTIntegratedLightYield
// ----------------------------------------------------------------------------

ACTIntegratedLightYield::ACTIntegratedLightYield(double w0)
  : Interpolation1D(), m_w0(w0)
{
  // nothing to see here
}

double ACTIntegratedLightYield::yield(double h, double w) const
{
  yield_t _y = y(h);
  double dw = w-m_w0;
#ifdef ACT_LIGHT_YIELD_TAYLOR_SERIES_IN_LOG
  return std::exp(_y.n + dw*(_y.dn_dw + _y.d2n_dw2*dw));
#else
  return _y.n + dw*(_y.dn_dw + _y.d2n_dw2*dw);
#endif
}
