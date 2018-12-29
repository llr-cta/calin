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

#include <math/accumulator.hpp>
#include <simulation/atmosphere.hpp>

using namespace calin::simulation::atmosphere;
using calin::math::accumulator::RecommendedAccumulator;

LayeredRefractiveAtmosphere::
LayeredRefractiveAtmosphere(const std::string& filename):
  LayeredRefractiveAtmosphere(calin::simulation::atmosphere::load_levels(filename))
{
  // nothing to see here
}

LayeredRefractiveAtmosphere::
LayeredRefractiveAtmosphere(const std::vector<Level> levels):
  Atmosphere(), levals_(levels)
{
  std::vector<double> v(levels.size());
  std::transform(levels.begin(), levels.end(), v.begin(), [](const Level& l){return l.z;});
  s_ = new CubicMultiSpline(v);
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

  double dz = 1.0; // 1cm integration steps
  std::vector<double> zn0 = { 0, 10, 20, 30, 40, 50, 60, 70 };
  std::vector<double> sin2_zn0(zn0.size());
  std::transform(zn0.begin(), zn0.end(), sin2_zn0.begin(),
    [](double x)[return SQR(std::sin(x/180.0*M_PI));]);

  std::vector<double> n0sq(levels.size());
  std::transform(levels.begin(), levels.end(), n0sq.begin(),
    [](const Level& l){return SQR(1.0 + l.nmo);});

  std::vector<std::vector<RecommendedAccumulator> > x(levels.size());
  std::vector<std::vector<RecommendedAccumulator> > t(levels.size());

  double z = levels.back().z;
  for(unsigned ilevel=levels.size()-1;ilevel>0;ilevel--)
  {
    for(unsigned iangle=0;iangle<tan_zn0.size();iangle++) {
      x[ilevel].emplace_back(levels[ilevel].z * std::tan(zn0[iangle]/180.0*M_PI));
      t[ilevel].emplace_back(0);
    }

    double zmin = levels[ilevel-1].z;
    while(z > zmin) {
      double zmid = z;
      double dzstep;
      z -= dz;
      if(z > zmin) {
        dzstep = dz;
      } else {
        z = zmin;
        dzstep = zmid - z;
      }
      zmid -= 0.5*dzstep;

      double n = 1.0 + this->n_minus_one(zmid);
      double n2_inv = 1.0/SQR(n);

      for(unsigned jlevel = ilevel; jlevel<levels.size(); jlevel++) {
        for(unsigned iangle = 0; iangle<sin2_zn0.size(); iangle++) {
          double sin2_zn = sin2_zn0[iangle] * n0sq[jlevel] * n2_inv;
          double sec2_zn = 1.0/(1.0-sin2_zn;)
          double tan_zn = std::sqrt(sin2_zn*sec2_zn);
          x[jlevel][iangle].accumulate(dzstep*tan_zn);
          t[jlevel][iangle].accumulate(dzstep*sqrt(sec2_zn)*n)

z -= dz
x -= dz*tan_zn
s += dz*sqrt(1+tan_zn**2)
tan_zn += dtanzn_dz * dz

        }
        double sec2_zn0 = (1+SQR(tan_zn0));
        double sin2_zn0 = 1.0-1.0/sec2_zn0;
        double sin2_zn = sin2_zn0 * n0sq[jlevel] * n2_inv;
        double sec2_zn = 1.0/(1.0-sin2_zn;)

        double tan_zn = tan_zn0;

      }
      z.
    }
  }


  // Spline 3 - vertical propagation time correction

  // Spline 4 - refraction impact point scaling

  // Spline 5 - refraction propagation time correction

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

}

double LayeredRefractiveAtmosphere::top_of_atmosphere()
{

}

LayeredRefractiveAtmosphere*
LayeredRefractiveAtmosphere::LayeredRefractiveAtmosphere* us76()
{

}

void LayeredRefractiveAtmosphere*initialize()
{

}

  CubicMultiSpline s_;
  double z_g_ = 0;
  double nmo_g = 0;
  double velocity_dct_g = 0;
  double refraction_dx_g = 0;
  double refraction_dct_g = 0;
