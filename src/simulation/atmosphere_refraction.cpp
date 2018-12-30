/*

   calin/simulation/atmosphere_refraction.cpp -- Stephen Fegan -- 2018-12-24

   Classes to model density and refractive index of atmosphere and calculate
   effect of refraction

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <cmath>

#include <calin_global_definitions.hpp>
#include <math/accumulator.hpp>
#include <math/special.hpp>
#include <simulation/atmosphere.hpp>

using namespace calin::simulation::atmosphere;
using calin::math::accumulator::RecommendedAccumulator;
using namespace calin::math::special;

LayeredRefractiveAtmosphere::
LayeredRefractiveAtmosphere(const std::string& filename):
  LayeredRefractiveAtmosphere(calin::simulation::atmosphere::load_levels(filename))
{
  // nothing to see here
}

LayeredRefractiveAtmosphere::
LayeredRefractiveAtmosphere(const std::vector<Level> levels):
  Atmosphere(), levels_(levels)
{
  std::vector<double> v(levels.size());
  std::transform(levels.begin(), levels.end(), v.begin(), [](const Level& l){return l.z;});
  s_ = new math::spline_interpolation::CubicMultiSpline(v);
  // Spline 0 - density
  std::transform(levels.begin(), levels.end(), v.begin(), [](const Level& l){return std::log(l.rho);});
  s_->add_spline(v, "density [ln(g/cm^3)]");

  // Spline 1 - integrated thickness to altitude
  std::transform(levels.begin(), levels.end(), v.begin(), [](const Level& l){return std::log(l.t);});
  s_->add_spline(v, "thickness to altitude [ln(g/cm^2)]");

  // Spline 2 - refractive index : n minus one
  std::transform(levels.begin(), levels.end(), v.begin(), [](const Level& l){return std::log(l.nmo);});
  s_->add_spline(v, "refractive index minus one [ln(1)]");

  // Calculate the effects of refraction by integrating the path of rays from
  // various levels in the atmosphere

  double dzn = 0.5 * M_PI/180.0/60.0/60.0/1000.0; // 0.1 milli-arcsecond integration step for ray at 45 degrees
  double dz_max = 1000; // 10m maximum integration steps

  std::vector<double> zn0 = { 0, 45, 10, 20, 30, 40, 50, 60, 70 };
  double tan_znref_inv = 1.0/std::tan(zn0[1]/180.0*M_PI);
  std::vector<double> sin2_zn0(zn0.size());
  std::transform(zn0.begin(), zn0.end(), sin2_zn0.begin(),
    [](double x){ return SQR(std::sin(x/180.0*M_PI)); });

  std::vector<double> n0sq(levels.size());
  std::transform(levels.begin(), levels.end(), n0sq.begin(),
    [](const Level& l){return SQR(1.0 + l.nmo);});

  std::vector<std::vector<RecommendedAccumulator> > x(levels.size());
  std::vector<std::vector<RecommendedAccumulator> > t(levels.size());

  double z = levels.back().z;

  double n;
  double dn_dz = this->dn_dz(z, n);
  n += 1.0;
  double dz = std::min(dz_max, dzn*n*tan_znref_inv/std::abs(dn_dz));

  for(unsigned ilevel=levels.size()-1;ilevel>0;ilevel--)
  {
    for(unsigned iangle=0;iangle<zn0.size();iangle++) {
      x[ilevel].emplace_back(0);
      t[ilevel].emplace_back(0);
    }

    unsigned nstep = 0;
    double zmin = levels[ilevel-1].z;
    while(z > zmin) {
      ++nstep;
      double zmid = z;
      double dzstep;
      z -= dz;
      if(z > zmin) {
        dzstep = dz;
      } else {
        z = zmin;
        dzstep = zmid - z;
        std::cout << z*1e-5 << ' ' <<  dz << ' ' << nstep << '\n';
        nstep = 0;
      }
      zmid -= 0.5*dzstep;

      dn_dz = this->dn_dz(zmid, n  /* n minus one */);
      n += 1.0;
      dz = std::min(dz_max, dzn*n*tan_znref_inv/std::abs(dn_dz));

      double n2_inv = 1.0/SQR(n);

      for(unsigned jlevel = ilevel; jlevel<levels.size(); jlevel++) {
        for(unsigned iangle = 0; iangle<sin2_zn0.size(); iangle++) {
          double sin2_zn = sin2_zn0[iangle] * n0sq[jlevel] * n2_inv;
          double sec2_zn = 1.0/(1.0-sin2_zn);
          double tan_zn = std::sqrt(sin2_zn*sec2_zn);
          x[jlevel][iangle].accumulate(dzstep*tan_zn);
          t[jlevel][iangle].accumulate(dzstep*sqrt(sec2_zn)*n);
        }
      }
    }
  }
  for(unsigned iangle=0;iangle<zn0.size();iangle++) {
    x[0].emplace_back(0);
    t[0].emplace_back(0);
  }

  for(unsigned ilevel = 1; ilevel<levels.size(); ilevel++) {
    for(unsigned iangle = 0; iangle<sin2_zn0.size(); iangle++) {
      x[ilevel][iangle].accumulate(-levels[ilevel].z * std::tan(zn0[iangle]/180.0*M_PI));
      t[ilevel][iangle].accumulate(-levels[ilevel].z / std::cos(zn0[iangle]/180.0*M_PI));
    }
  }

  test_ray_zne_ = calin::std_to_eigenvec(zn0);
  test_ray_ze_ = calin::std_to_eigenvec(s_->xknot());
  test_ray_xg_.resize(levels.size(),zn0.size());
  test_ray_ctg_.resize(levels.size(),zn0.size());
  for(unsigned ilevel = 0; ilevel<levels.size(); ilevel++) {
    for(unsigned iangle = 0; iangle<sin2_zn0.size(); iangle++) {
      test_ray_xg_(ilevel,iangle) = x[ilevel][iangle];
      test_ray_ctg_(ilevel,iangle) = t[ilevel][iangle];
    }
  }

  // Spline 3 - vertical propagation time correction
  std::transform(t.begin(), t.end(), v.begin(),
    [](const std::vector<RecommendedAccumulator>& v){return v[0];});
  s_->add_spline(v, "vertical propagation ct correction [cm]");

  // Spline 4 - refraction impact point scaling
  std::transform(x.begin(), x.end(), v.begin(),
    [](const std::vector<RecommendedAccumulator>& v){return v[1];});
  s_->add_spline(v, "refraction impact dist correction [cm]");

  // Spline 5 - refraction propagation time correction
  std::transform(t.begin(), t.end(), v.begin(),
    [zn0](const std::vector<RecommendedAccumulator>& v){return v[1]-v[0]/std::cos(zn0[1]/180.0*M_PI);});
  s_->add_spline(v, "refraction impact ct correction [cm]");

}

LayeredRefractiveAtmosphere::~LayeredRefractiveAtmosphere()
{
  delete s_;
}

double LayeredRefractiveAtmosphere::rho(double z)
{
  return exp(s_->value(z, 0));
}

double LayeredRefractiveAtmosphere::thickness(double z)
{
  return exp(s_->value(z, 1));
}

double LayeredRefractiveAtmosphere::n_minus_one(double z)
{
  return exp(s_->value(z, 2));
}

double LayeredRefractiveAtmosphere::dn_dz(double z, double& n_minus_one)
{
  double dlogn_dz = s_->derivative_and_value(z, 2, n_minus_one);
  n_minus_one = exp(n_minus_one);
  return n_minus_one * dlogn_dz;
}

double LayeredRefractiveAtmosphere::propagation_ct_correction(double z)
{
  return exp(s_->value(z, 3));
}

void LayeredRefractiveAtmosphere::cherenkov_parameters(double z,
  double& n_minus_one, double& propagation_ct_correction)
{
  n_minus_one = exp(s_->value(z, 2));
  propagation_ct_correction = exp(s_->value(z, 3));
}

void LayeredRefractiveAtmosphere::
propagate_ray_with_refraction(const calin::math::ray::Ray ray)
{

}

double LayeredRefractiveAtmosphere::z_for_thickness(double t)
{
  return 0;
}

double LayeredRefractiveAtmosphere::top_of_atmosphere()
{
  return s_->xknot().back();
}

LayeredRefractiveAtmosphere*
LayeredRefractiveAtmosphere::LayeredRefractiveAtmosphere::us76()
{
  return new LayeredRefractiveAtmosphere(us76_levels());
}

void LayeredRefractiveAtmosphere::initialize()
{

}
