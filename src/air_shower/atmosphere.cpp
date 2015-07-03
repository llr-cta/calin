/* 

   calin/air_shower/atmosphere.cpp -- Stephen Fegan -- 2015-06-11

   Classes to model density and refractive index of atmosphere

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

#include"air_shower/atmosphere.hpp"
#include"air_shower/atmosphere.hpp"

using namespace calin::air_shower::atmosphere;

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
  m_zs(zs), m_zmax(zmax), m_nmo0(nmo0), temperature_(temperature)
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

inline std::vector<calin::air_shower::atmosphere::LayeredAtmosphereLayer>::const_iterator
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
