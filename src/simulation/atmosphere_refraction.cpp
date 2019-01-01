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
LayeredRefractiveAtmosphere(const std::string& filename,
    const std::vector<double>& obs_levels):
  LayeredRefractiveAtmosphere(calin::simulation::atmosphere::load_levels(filename), obs_levels)
{
  // nothing to see here
}

LayeredRefractiveAtmosphere::
LayeredRefractiveAtmosphere(const std::string& filename, double obs_level):
  LayeredRefractiveAtmosphere(filename, std::vector<double>{ obs_level })
{
  // nothing to see here
}

LayeredRefractiveAtmosphere::
LayeredRefractiveAtmosphere(const std::vector<Level>& levels, const std::vector<double>& obs_levels):
  Atmosphere(), levels_(levels), zobs_(obs_levels)
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
  // various levels in the atmosphere. We use a matrix of test rays that have
  // predefined zenith angles at each of the atmospheric layers, but we track
  // them all from the top of the atmosphere - we use Snell's law to arrange
  // that they all have correct zenith angle at their prescribed layer. We
  // record the position and time when each ray passes its origin layer in the
  // atmpsohere (test_ray_lev_x_ and test_ray_lev_ct_) and when each of them passes
  // each of the ground layers (test_ray_obs_x_ and test_ray_obs_ct_).

  std::vector<double> zn0 = { 0, 45, 10, 20, 30, 40, 50, 60, 70 };

  // Set up ground levels - must be monotonic and within atmosphere
  if(zobs_.empty()) {
    zobs_.emplace_back(levels_.front().z);
  } else {
    double zg = levels_.front().z;
    if(zobs_.front() < zg) {
      throw std::runtime_error(
        "Observation-level altitudes must be above bottom of the atmospehere");
    }
    for(double z : zobs_) {
      if(z<zg)throw std::runtime_error(
        "Observation-level altitudes must be monotonic increasing");
      zg = z;
    }
    if(levels_.back().z < zg) {
      throw std::runtime_error(
        "Observation-level altitudes must be below top of the atmosphere");
    }
  }

  // Storage for results of integration at atmospheric and observation levels
  test_ray_lev_x_.resize(levels.size(),zn0.size());
  test_ray_lev_ct_.resize(levels.size(),zn0.size());
  test_ray_obs_x_.resize(zobs_.size(), { levels.size(),zn0.size() });
  test_ray_obs_ct_.resize(zobs_.size(), { levels.size(),zn0.size() });

  // Integration step
  double dzn = 0.5 * M_PI/180.0/60.0/60.0/1000.0; // 0.5 milli-arcsecond integration step for ray at 45 degrees
  double dz_max = 1000; // 10m maximum integration steps

  // Various derived values to speed up calculation
  double tan_znref_inv = 1.0/std::tan(zn0[1]/180.0*M_PI);
  std::vector<double> sin2_zn0(zn0.size());
  std::transform(zn0.begin(), zn0.end(), sin2_zn0.begin(),
    [](double x){ return SQR(std::sin(x/180.0*M_PI)); });

  std::vector<double> n0sq(levels.size());
  std::transform(levels.begin(), levels.end(), n0sq.begin(),
    [](const Level& l){return SQR(1.0 + l.nmo);});

  // Track position (x) and c * time (t) of ray as it propagates .. initialize
  // all test rays to x=0 and ct=0 at top of atmosphere
  std::vector<std::vector<RecommendedAccumulator> > x(levels.size());
  std::vector<std::vector<RecommendedAccumulator> > t(levels.size());

  for(unsigned ilevel=0;ilevel<levels.size();ilevel++) {
    for(unsigned iangle=0;iangle<zn0.size();iangle++) {
      x[ilevel].emplace_back(0);
      t[ilevel].emplace_back(0);
    }
  }

  // Initialize z at top of atmosphere and calculate step in z
  double z = levels.back().z;

  double n;
  double dn_dz = this->dn_dz(z, n);
  n += 1.0;
  double dz = std::min(dz_max, dzn*n*tan_znref_inv/std::abs(dn_dz));

  int iobs = zobs_.size()-1;

  for(unsigned ilevel=levels.size()-1;ilevel>0;ilevel--)
  {
    for(unsigned iangle=0;iangle<zn0.size();iangle++) {
      test_ray_lev_x_(ilevel,iangle) = x[ilevel][iangle].total();
      test_ray_lev_ct_(ilevel,iangle) = t[ilevel][iangle].total();
    }

    unsigned nstep = 0;
    double zmin = levels[ilevel-1].z;
    while(z > zmin) {
      ++nstep;
      double zold = z;
      double dzstep;
      z -= dz;
      if(z > zmin) {
        dzstep = dz;
      } else {
        z = zmin;
        dzstep = zold - z;
        std::cout << z*1e-5 << ' ' <<  dz << ' ' << nstep << '\n';
      }
      double zmid = zold - 0.5*dzstep;

      dn_dz = this->dn_dz(zmid, n  /* n minus one */);
      n += 1.0;
      dz = std::min(dz_max, dzn*n*tan_znref_inv/std::abs(dn_dz));

      double n2_inv = 1.0/SQR(n);

      while(iobs>=0 and zold>zobs_[iobs] and z<=zobs_[iobs]) {
        // This integraion step intersects an observation level so store
        // relevant data. This rarely happens so there is no penalty from
        // recalculating the refracted angles here

        double dzobsstep = zold-zobs_[iobs];
        for(unsigned jlevel = 0; jlevel<levels.size(); jlevel++) {
          for(unsigned iangle = 0; iangle<sin2_zn0.size(); iangle++) {
            double sin2_zn = sin2_zn0[iangle] * n0sq[jlevel] * n2_inv;
            double sec2_zn = 1.0/(1.0-sin2_zn);
            double tan_zn = std::sqrt(sin2_zn*sec2_zn);
            test_ray_obs_x_[iobs](jlevel,iangle) = x[jlevel][iangle] + dzobsstep*tan_zn;
            test_ray_obs_ct_[iobs](jlevel,iangle) = t[jlevel][iangle] + dzobsstep*sqrt(sec2_zn)*n;
          }
        }
        iobs--;
      }

      for(unsigned jlevel = 0; jlevel<levels.size(); jlevel++) {
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
    test_ray_lev_x_(0,iangle) = x[0][iangle].total();
    test_ray_lev_ct_(0,iangle) = t[0][iangle].total();
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
LayeredRefractiveAtmosphere::LayeredRefractiveAtmosphere::us76(const std::vector<double>& obs_levels)
{
  return new LayeredRefractiveAtmosphere(us76_levels(), obs_levels);
}

void LayeredRefractiveAtmosphere::initialize()
{

}
