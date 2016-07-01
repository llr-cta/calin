/*

   calin/simulation/atmosphere.cpp -- Stephen Fegan -- 2015-06-11

   Classes to model density and refractive index of atmosphere

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

// Originally: Atmposphere.cpp - Classes to handle atmosphere
// Stephen Fegan - sfegan@llr.in2p3.fr - July 2012
// Id: Atmosphere.cpp 4710 2012-11-06 15:00:55Z sfegan

#include<cmath>
#include<cctype>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<algorithm>

#include <simulation/atmosphere.hpp>

using namespace calin::simulation::atmosphere;

Atmosphere::~Atmosphere()
{
  // nothing to see here
}

std::vector<AtmSlice> Atmosphere::make_atm_slices(unsigned nslice)
{
  return make_atm_slices(nslice,this->top_of_atmosphere(),0);
}

std::vector<AtmSlice> Atmosphere::make_atm_slices(unsigned nslice,
                                                  double zmax, double zmin)
{
  std::vector<AtmSlice> atm(nslice);

  double tmax = this->thickness(zmin);
  double tmin = this->thickness(zmax);
  double dt = tmax-tmin;

  double z = zmin;
  for(unsigned islice=0;islice<nslice;islice++)
    {
      AtmSlice& s(atm[islice]);
      double t = dt * double(nslice-islice)/double(nslice);
      s.zb = z;
      s.tb = t;
      t = dt * double(nslice-islice-1)/double(nslice);
      z = this->z_for_thickness(t+tmin);
      s.zt = z;
      s.tt = t;
      s.rho = (s.tb-s.tt)/(s.zt-s.zb);
    }

  return atm;
}

// ****************************************************************************
// ISOTHERMAL ATMOSPHERE
// ****************************************************************************

IsothermalAtmosphere::
IsothermalAtmosphere(double rho0, double zs, double zmax, double nmo0,
                     double temperature):
  Atmosphere(), m_ttoa(rho0*zs*std::exp(-zmax/zs)), m_rho0(rho0),
  m_zs(zs), m_zmax(zmax), m_nmo0(nmo0) /*, temperature_(temperature)*/
{
  // nothing to see here
}

IsothermalAtmosphere::~IsothermalAtmosphere()
{
  // nothing to see here
}

double IsothermalAtmosphere::rho(double z)
{
  return m_rho0*std::exp(-z/m_zs);
}

double IsothermalAtmosphere::thickness(double z)
{
  return m_rho0*m_zs*std::exp(-z/m_zs)-m_ttoa;
}

double IsothermalAtmosphere::n_minus_one(double z)
{
  return m_nmo0*std::exp(-z/m_zs);
}

double IsothermalAtmosphere::propagation_time_correction(double z)
{
  return m_zs*m_nmo0*(1.0-exp(-z/m_zs));
}

double IsothermalAtmosphere::z_for_thickness(double t)
{
  return -std::log((t+m_ttoa)/(m_rho0*m_zs))*m_zs;
}

double IsothermalAtmosphere::top_of_atmosphere()
{
  return m_zmax;
}

// ****************************************************************************
// LAYERED ATMOSPHERE
// ****************************************************************************

LayeredAtmosphere::LayeredAtmosphere(const std::string& filename):
  Atmosphere(), m_ztoa(), m_ttoa(), m_levels(), m_layers(), m_ilayer()
{
  std::ifstream stream(filename.c_str());
  if(!stream)
    throw std::string("LayeredAtmosphere: could not open: ")+filename;

  std::string line;
  while(std::getline(stream, line))
    {
      unsigned ichar=0;
      while(isspace(line[ichar]))ichar++;
      if(line[ichar] == '#')continue;
      std::istringstream lstream(line);
      Level l;
      lstream >> l.z >> l.rho >> l.t >> l.nmo;
      l.z *= 1e5;
      m_levels.push_back(l);
    }
  initialize();
}

LayeredAtmosphere::LayeredAtmosphere(const std::vector<Level> levels):
  Atmosphere(), m_ztoa(), m_ttoa(), m_levels(levels), m_layers(), m_ilayer()
{
  initialize();
}

void LayeredAtmosphere::initialize()
{
  if(m_levels.size()<2)
    throw std::string("LayeredAtmosphere: A minimum of 2 levels required.");
  std::sort(m_levels.begin(), m_levels.end(), Level::CmpZAsc());
  m_layers.resize(m_levels.size()-1);

  m_ztoa = m_levels.back().z;
  m_ttoa = m_levels.back().t;

  if(m_ttoa <= 0.0)
    {
      // If the thickness at the top of the atmosphere is zero
      // (i.e. the thickness from infinity to that point has been
      // subtracted) then solve for the thickness there by assuming
      // that the scale height between the final three levels is
      // constant

      if(m_levels.size()<3)
	throw std::string("LayeredAtmosphere: A minimum of 3 levels required "
			  "to solve for thickness.");

      const double y2 = m_levels[m_levels.size() - 2].t;
      const double y3 = m_levels[m_levels.size() - 3].t;
      const double x1 = m_levels[m_levels.size() - 1].z;
      const double x2 = m_levels[m_levels.size() - 2].z;
      const double x3 = m_levels[m_levels.size() - 3].z;

      double H = (x2-x3)/(std::log(y3)-std::log(y2)); // initial guess
      double num = std::exp(-x3/H)-std::exp(-x1/H);
      double den = std::exp(-x2/H)-std::exp(-x1/H);
      double df = y3/y2 - num/den;
      unsigned niter = 10;
      while(std::fabs(df)/(y3/y2) > 1e-8)
	{
	  // Newton-Ralphson to find value of H giving agreement with data
	  if(niter-- == 0)
	    throw std::string("LayeredAtmosphere: max number of iterations exceeded");
	  double dfdH =
	    (den*(x3*std::exp(-x3/H)-x1*std::exp(-x1/H))
	     - num*(x2*std::exp(-x2/H)-x1*std::exp(-x1/H)))/(den*den*H*H);
	  H += df/dfdH;
	  num = std::exp(-x3/H)-std::exp(-x1/H);
	  den = std::exp(-x2/H)-std::exp(-x1/H);
	  df = y3/y2 - num/den;
	}

      double t0 = y3/(std::exp(-x3/H)-std::exp(-x1/H));
      m_ttoa = t0*std::exp(-x1/H);

      for(unsigned ilevel=0;ilevel<m_levels.size();ilevel++)
	m_levels[ilevel].t += m_ttoa;
    }

  double ptc = 0;
  for(unsigned ilayer=0;ilayer<m_layers.size();ilayer++)
    {
      Layer& layer(m_layers[ilayer]);
      const Level& t(m_levels[ilayer+1]);
      const Level& b(m_levels[ilayer]);
      layer.zb    = b.z;
      layer.zt    = t.z;
      layer.rhozs = -(t.z-b.z)/(std::log(t.rho)-std::log(b.rho));
      layer.rho0  = b.rho/std::exp(-b.z/layer.rhozs);
      layer.tzs   = -(t.z-b.z)/(std::log(t.t)-std::log(b.t));
      layer.t0    = b.t/std::exp(-b.z/layer.tzs);
      layer.tb    = b.t;
      layer.tt    = t.t;
      layer.nmozs = -(t.z-b.z)/(std::log(t.nmo)-std::log(b.nmo));
      layer.nmo0  = b.nmo/std::exp(-b.z/layer.nmozs);
      ptc += layer.nmozs*layer.nmo0*std::exp(-b.z/layer.nmozs);
      layer.ptc0   = ptc;
      ptc -= layer.nmozs*layer.nmo0*std::exp(-t.z/layer.nmozs);
    }
  m_ilayer = m_layers.end()-1;
}

LayeredAtmosphere::~LayeredAtmosphere()
{
  // nothing to see here
}

inline std::vector<calin::simulation::atmosphere::LayeredAtmosphereLayer>::const_iterator
LayeredAtmosphere::findZ(double z) const
{
  std::vector<Layer>::const_iterator ilayer = m_ilayer;

  //std::cout << "A: " << z << ' ' << ilayer-m_layers.begin() << ' ' << ilayer->zb << ' ' << ilayer->zt << '\n';

  if(z<=ilayer->zt)
  {
    if(z<=ilayer->zb && ilayer!=m_layers.begin())
    {
      ilayer--;
      if(z<=ilayer->zb && ilayer!=m_layers.begin())
      {
        ilayer = std::lower_bound(m_layers.begin(), ilayer, z);
        if(z<=ilayer->zb)ilayer = m_layers.begin();
      }
    }
  }
  else
  {
    if(++ilayer == m_layers.end())ilayer--;
    else if(ilayer->zt < z)
    {
      ++ilayer;
      ilayer = std::lower_bound(ilayer,m_layers.end(), z);
      if(ilayer == m_layers.end())ilayer--;
    }
  }

  //std::cout << "B: " << z << ' ' << ilayer-m_layers.begin() << ' ' << ilayer->zb << ' ' << ilayer->zt << '\n';

  assert(ilayer<m_layers.end());
  assert(z<=ilayer->zt || ilayer==m_layers.end()-1);
  assert(ilayer->zb<z || ilayer==m_layers.begin());

  m_ilayer = ilayer;
  return ilayer;
}

double LayeredAtmosphere::rho(double z)
{
  auto ilayer = findZ(z);
  return ilayer->rho0 * std::exp(-z/ilayer->rhozs);
}

double LayeredAtmosphere::thickness(double z)
{
  auto ilayer = findZ(z);
  return ilayer->t0 * std::exp(-z/ilayer->tzs) - m_ttoa;
}

double LayeredAtmosphere::n_minus_one(double z)
{
  auto ilayer = findZ(z);
  return ilayer->nmo0 * std::exp(-z/ilayer->nmozs);
}

double LayeredAtmosphere::propagation_time_correction(double z)
{
  auto ilayer = findZ(z);
  return ilayer->ptc0 -
      ilayer->nmozs * ilayer->nmo0 * std::exp(-z/ilayer->nmozs);
}

double LayeredAtmosphere::z_for_thickness(double t)
{
  t += m_ttoa;
  std::vector<Layer>::const_iterator ilayer
    = std::lower_bound(m_layers.begin(),m_layers.end(), t, Layer::CmpTDec());
  if(ilayer == m_layers.end())ilayer--;
  return -std::log(t/ilayer->t0)*ilayer->tzs;
}

double LayeredAtmosphere::top_of_atmosphere()
{
  return m_ztoa;
}

 LayeredAtmosphere* LayeredAtmosphere::us76()
{
  std::vector<Level> levels;
  levels.push_back({ 0.0, 0.12219E-02, 0.10350E+04, 0.28232E-03 });
  levels.push_back({ 100000.0, 0.11099E-02, 0.91853E+03, 0.25634E-03 });
  levels.push_back({ 200000.0, 0.10054E-02, 0.81286E+03, 0.23214E-03 });
  levels.push_back({ 300000.0, 0.90839E-03, 0.71725E+03, 0.20975E-03 });
  levels.push_back({ 400000.0, 0.81888E-03, 0.63097E+03, 0.18904E-03 });
  levels.push_back({ 500000.0, 0.73643E-03, 0.55328E+03, 0.16994E-03 });
  levels.push_back({ 600000.0, 0.66012E-03, 0.48352E+03, 0.15235E-03 });
  levels.push_back({ 700000.0, 0.59048E-03, 0.42105E+03, 0.13620E-03 });
  levels.push_back({ 800000.0, 0.52609E-03, 0.36529E+03, 0.12136E-03 });
  levels.push_back({ 900000.0, 0.46741E-03, 0.31567E+03, 0.10782E-03 });
  levels.push_back({ 1000000.0, 0.41370E-03, 0.27167E+03, 0.95426E-04 });
  levels.push_back({ 1100000.0, 0.36499E-03, 0.23278E+03, 0.84194E-04 });
  levels.push_back({ 1200000.0, 0.31209E-03, 0.19900E+03, 0.71987E-04 });
  levels.push_back({ 1300000.0, 0.26674E-03, 0.17012E+03, 0.61523E-04 });
  levels.push_back({ 1400000.0, 0.22792E-03, 0.14543E+03, 0.52581E-04 });
  levels.push_back({ 1500000.0, 0.19479E-03, 0.12434E+03, 0.44937E-04 });
  levels.push_back({ 1600000.0, 0.16651E-03, 0.10631E+03, 0.38406E-04 });
  levels.push_back({ 1700000.0, 0.14236E-03, 0.90902E+02, 0.32840E-04 });
  levels.push_back({ 1800000.0, 0.12168E-03, 0.77727E+02, 0.28071E-04 });
  levels.push_back({ 1900000.0, 0.10403E-03, 0.66465E+02, 0.23997E-04 });
  levels.push_back({ 2000000.0, 0.88928E-04, 0.56837E+02, 0.20516E-04 });
  levels.push_back({ 2100000.0, 0.75750E-04, 0.48620E+02, 0.17475E-04 });
  levels.push_back({ 2200000.0, 0.64544E-04, 0.41621E+02, 0.14887E-04 });
  levels.push_back({ 2300000.0, 0.55021E-04, 0.35655E+02, 0.12695E-04 });
  levels.push_back({ 2400000.0, 0.46965E-04, 0.30566E+02, 0.10833E-04 });
  levels.push_back({ 2500000.0, 0.40097E-04, 0.26222E+02, 0.92494E-05 });
  levels.push_back({ 2750000.0, 0.27126E-04, 0.17925E+02, 0.62570E-05 });
  levels.push_back({ 3000000.0, 0.18420E-04, 0.12302E+02, 0.42495E-05 });
  levels.push_back({ 3250000.0, 0.12139E-04, 0.85361E+01, 0.28004E-05 });
  levels.push_back({ 3500000.0, 0.84696E-05, 0.59874E+01, 0.19537E-05 });
  levels.push_back({ 3750000.0, 0.59542E-05, 0.42029E+01, 0.13738E-05 });
  levels.push_back({ 4000000.0, 0.39967E-05, 0.29752E+01, 0.92196E-06 });
  levels.push_back({ 4250000.0, 0.27910E-05, 0.21358E+01, 0.64379E-06 });
  levels.push_back({ 4500000.0, 0.19671E-05, 0.15470E+01, 0.45379E-06 });
  levels.push_back({ 4750000.0, 0.14044E-05, 0.11295E+01, 0.32390E-06 });
  levels.push_back({ 5000000.0, 0.10273E-05, 0.82800E+00, 0.23698E-06 });
  levels.push_back({ 5500000.0, 0.56800E-06, 0.44045E+00, 0.13104E-06 });
  levels.push_back({ 6000000.0, 0.30906E-06, 0.22771E+00, 0.71295E-07 });
  levels.push_back({ 6500000.0, 0.16285E-06, 0.11361E+00, 0.37569E-07 });
  levels.push_back({ 7000000.0, 0.82868E-07, 0.54414E-01, 0.19114E-07 });
  levels.push_back({ 7500000.0, 0.40145E-07, 0.24940E-01, 0.92604E-08 });
  levels.push_back({ 8000000.0, 0.18430E-07, 0.10993E-01, 0.42513E-08 });
  levels.push_back({ 8500000.0, 0.82291E-08, 0.46676E-02, 0.18985E-08 });
  levels.push_back({ 9000000.0, 0.34321E-08, 0.19250E-02, 0.79163E-09 });
  levels.push_back({ 9500000.0, 0.14063E-08, 0.78968E-03, 0.32437E-09 });
  levels.push_back({ 10000000.0, 0.57185E-09, 0.32602E-03, 0.13189E-09 });
  levels.push_back({ 10500000.0, 0.24206E-09, 0.13421E-03, 0.55841E-10 });
  levels.push_back({ 11000000.0, 0.10312E-09, 0.52792E-04, 0.23788E-10 });
  levels.push_back({ 11500000.0, 0.46595E-10, 0.17216E-04, 0.10748E-10 });
  levels.push_back({ 12000000.0, 0.24596E-10, 0.00000E+00, 0.56734E-11 });
  return new LayeredAtmosphere(levels);
}
