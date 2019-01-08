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

//                             _              __  __           _      _
//     /\                     | |            |  \/  |         | |    | |
//    /  \   _ __   __ _ _   _| | __ _ _ __  | \  / | ___   __| | ___| |
//   / /\ \ | '_ \ / _` | | | | |/ _` | '__| | |\/| |/ _ \ / _` |/ _ \ |
//  / ____ \| | | | (_| | |_| | | (_| | |    | |  | | (_) | (_| |  __/ |
// /_/    \_\_| |_|\__, |\__,_|_|\__,_|_|    |_|  |_|\___/ \__,_|\___|_|
//                  __/ |
//                 |___/

// The angular model consists of two terms in sin_i, cos_i, sin_r and cos2_r
// that can be computed for rays without using trignometric functions

namespace {
  inline void calculate_refraction_angular_terms(
    double sin_i, double cos_i, double n_i, double n_r_inv,
    double& sin_r, double& cos_r,
    double& t1_x, double& t2_x, double& t1_ct, double& t2_ct)
  {
    sin_r = sin_i * n_i * n_r_inv;
    const double cos2_r = 1.0 - SQR(sin_r);
    cos_r = std::sqrt(cos2_r);
    const double sec2_i = 1.0/SQR(cos_i);

    t1_x = sin_i * sec2_i / cos_r;
    t2_x = t1_x * SQR(sin_i) * sec2_i;

    t1_ct = t1_x * sin_i;
    t2_ct = t2_x * sin_i;
  }

  void fit_refraction_angular_coefficients(
    double z_i, double z_r, double n_i, double n_r_inv, const std::vector<double>& zn,
    Eigen::VectorXd x_obs, Eigen::VectorXd ct_obs,
    double& a_x, double& b_x, double& a_ct, double& b_ct,
    bool use_fixed_ba_ratio = false, double ba_ratio_x = 0, double ba_ratio_ct = 0)
  {
    const double z_prop = z_i-z_r;
    Eigen::Matrix2d M_x = Eigen::Matrix2d::Zero();
    Eigen::Vector2d V_x = Eigen::Vector2d::Zero();
    Eigen::Matrix2d M_ct = Eigen::Matrix2d::Zero();
    Eigen::Vector2d V_ct = Eigen::Vector2d::Zero();
    for(unsigned izn=0; izn<zn.size(); ++izn) {
      const double sin_i = std::sin(zn[izn]);
      const double cos_i = std::cos(zn[izn]);
      double sin_r;
      double cos_r;
      double t1_x;
      double t2_x;
      double t1_ct;
      double t2_ct;
      calculate_refraction_angular_terms(sin_i, cos_i, n_i, n_r_inv,
        cos_r, sin_r, t1_x, t2_x, t1_ct, t2_ct);
      double sec_i = 1.0/cos_i;

      double mc_x = z_prop*sin_i*sec_i - x_obs(izn);
      M_x(0,0) += t1_x*t1_x;
      M_x(1,0) += t1_x*t2_x;
      M_x(0,1) += t1_x*t2_x;
      M_x(1,1) += t2_x*t2_x;
      V_x(0) += mc_x*t1_x;
      V_x(1) += mc_x*t2_x;

      double mc_ct = ct_obs(0)*sec_i - ct_obs(izn);
      M_ct(0,0) += t1_ct*t1_ct;
      M_ct(1,0) += t1_ct*t2_ct;
      M_ct(0,1) += t1_ct*t2_ct;
      M_ct(1,1) += t2_ct*t2_ct;
      V_ct(0) += mc_ct*t1_ct;
      V_ct(1) += mc_ct*t2_ct;
    }

    if(use_fixed_ba_ratio) {
      a_x = (V_x(0) + ba_ratio_x * V_x(1)) /
        (M_x(0,0)*M_x(0,0) + 2*ba_ratio_x*M_x(0,0)*M_x(1,0) + ba_ratio_x*ba_ratio_x*M_x(1,1)*M_x(1,1));
      b_x = a_x * ba_ratio_x;
      a_ct = (V_ct(0) + ba_ratio_ct * V_ct(1)) /
        (M_ct(0,0)*M_ct(0,0) + 2*ba_ratio_ct*M_ct(0,0)*M_ct(1,0) + ba_ratio_ct*ba_ratio_ct*M_ct(1,1)*M_ct(1,1));
      b_ct = a_ct * ba_ratio_ct;
    } else {
      Eigen::Vector2d ab_x = M_x.colPivHouseholderQr().solve(V_x);
      Eigen::Vector2d ab_ct = M_ct.colPivHouseholderQr().solve(V_ct);
      a_x = ab_x(0);
      b_x = ab_x(1);
      a_ct = ab_ct(0);
      b_ct = ab_ct(1);
    }
    // std::cout << M_x << '\n' << V_x << '\n' << M_ct << '\n' << V_ct << '\n';
    // std::cout << a_x << ' ' << b_x << ' ' << a_ct << ' ' << b_ct << '\n';
  }
} // anonymous namespace

calin::ix::simulation::atmosphere::LayeredRefractiveAtmosphereConfig
LayeredRefractiveAtmosphere::default_config()
{
  calin::ix::simulation::atmosphere::LayeredRefractiveAtmosphereConfig config;
  config.set_zn_reference(45.0);
  for(double z=5.0; z<41.0; z+=5.0)config.add_zn_optimize(z);
  for(double z=50.0; z<71.0; z+=5.0)config.add_zn_optimize(z);
  config.set_step_delta_zn(0.5);
  config.set_step_delta_z_max(1000);
  config.set_high_accuracy_mode(false);
  config.set_angular_model_optimization_altitude(1800000.0);
  return config;
}

LayeredRefractiveAtmosphere::
LayeredRefractiveAtmosphere(const std::string& filename,
    const std::vector<double>& obs_levels,
    const calin::ix::simulation::atmosphere::LayeredRefractiveAtmosphereConfig& config):
  LayeredRefractiveAtmosphere(calin::simulation::atmosphere::load_levels(filename), obs_levels, config)
{
  // nothing to see here
}

LayeredRefractiveAtmosphere::
LayeredRefractiveAtmosphere(const std::string& filename, double obs_level,
const calin::ix::simulation::atmosphere::LayeredRefractiveAtmosphereConfig& config):
  LayeredRefractiveAtmosphere(filename, std::vector<double>{ obs_level }, config)
{
  // nothing to see here
}

LayeredRefractiveAtmosphere::
LayeredRefractiveAtmosphere(const std::vector<Level>& levels, const std::vector<double>& obs_levels,
const calin::ix::simulation::atmosphere::LayeredRefractiveAtmosphereConfig& config):
  Atmosphere(), levels_(levels), zobs_(obs_levels),
  high_accuracy_mode_(config.high_accuracy_mode())
{
  //          _                             _
  //     /\  | |                           | |
  //    /  \ | |_ _ __ ___   ___  ___ _ __ | |__   ___ _ __ ___
  //   / /\ \| __| '_ ` _ \ / _ \/ __| '_ \| '_ \ / _ \ '__/ _ \
  //  / ____ \ |_| | | | | | (_) \__ \ |_) | | | |  __/ | |  __/
  // /_/    \_\__|_| |_| |_|\___/|___/ .__/|_| |_|\___|_|  \___|
  //                                 | |
  //                                 |_|

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

  //  _____       _ _   _       _           _   _
  // |_   _|     (_) | (_)     (_)         | | (_)
  //   | |  _ __  _| |_ _  __ _ _ ___  __ _| |_ _  ___  _ __
  //   | | | '_ \| | __| |/ _` | / __|/ _` | __| |/ _ \| '_ \
  //  _| |_| | | | | |_| | (_| | \__ \ (_| | |_| | (_) | | | |
  // |_____|_| |_|_|\__|_|\__,_|_|___/\__,_|\__|_|\___/|_| |_|

  // Calculate the effects of refraction by integrating the path of rays from
  // various levels in the atmosphere. We use a matrix of test rays that have
  // predefined zenith angles at each of the atmospheric layers, but we track
  // them all from the top of the atmosphere - we use Snell's law to arrange
  // that they all have correct zenith angle at their prescribed layer. We
  // record the position and time when each ray passes its origin layer in the
  // atmpsohere (test_ray_emi_x_ and test_ray_emi_ct_) and when each of them
  // passes each of the ground layers (test_ray_obs_x_ and test_ray_obs_ct_).

  // We track all rays from the top of the atmosphere to the bottom so that
  // the resulting splines are smooth - otherwise we would potentially have
  // problems around each observation level.

  std::vector<double> zn0 = { 0, config.zn_reference() };
  for(auto zn : config.zn_optimize()) {
    if(zn>0 and zn<90 and zn!=config.zn_reference()) {
      zn0.push_back(zn);
    }
  }
  std::transform(zn0.begin(), zn0.end(), zn0.begin(),
    [](const double zn){return zn/180.0*M_PI;});

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
  nobs_inv_.resize(zobs_.size());
  std::transform(zobs_.begin(), zobs_.end(), nobs_inv_.begin(),
    [this](double z){return 1.0/(1.0+this->n_minus_one(z));});

  // Storage for results of integration at atmospheric and observation levels
  test_ray_emi_x_.resize(levels.size(),zn0.size());
  test_ray_emi_ct_.resize(levels.size(),zn0.size());
  test_ray_obs_x_.resize(zobs_.size(), { levels.size(),zn0.size() });
  test_ray_obs_ct_.resize(zobs_.size(), { levels.size(),zn0.size() });

  // Integration step
  double dzn = config.step_delta_zn() * M_PI/180.0/60.0/60.0/1000.0; // 0.5 milli-arcsecond integration step for ray at 45 degrees
  double dz_max = config.step_delta_z_max(); // 10m maximum integration steps

  // Various derived values to speed up calculation
  double tan_znref_inv = 1.0/std::tan(zn0[1]);
  std::vector<double> sin2_zn0(zn0.size());
  std::transform(zn0.begin(), zn0.end(), sin2_zn0.begin(),
    [](double x){ return SQR(std::sin(x)); });

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

  // Initialize z at top of atmosphere and calculate initial step in z
  double z = levels.back().z;

  double n;
  double dn_dz = this->dn_dz(z, n);
  n += 1.0;
  double dz = std::min(dz_max, dzn*n*tan_znref_inv/std::abs(dn_dz));

  //  _____       _                       _   _
  // |_   _|     | |                     | | (_)
  //   | |  _ __ | |_ ___  __ _ _ __ __ _| |_ _  ___  _ __
  //   | | | '_ \| __/ _ \/ _` | '__/ _` | __| |/ _ \| '_ \
  //  _| |_| | | | ||  __/ (_| | | | (_| | |_| | (_) | | | |
  // |_____|_| |_|\__\___|\__, |_|  \__,_|\__|_|\___/|_| |_|
  //                       __/ |
  //                      |___/

  // Integrate all test ray paths from top of atmosphere (x=ct=0) to bottom
  // recording positions at each atmospheric and observation level
  int iobs = zobs_.size()-1;
  for(unsigned ilevel=levels.size()-1;ilevel>0;ilevel--)
  {
    for(unsigned iangle=0;iangle<zn0.size();iangle++) {
      test_ray_emi_x_(ilevel,iangle) = x[ilevel][iangle].total();
      test_ray_emi_ct_(ilevel,iangle) = t[ilevel][iangle].total();
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
    test_ray_emi_x_(0,iangle) = x[0][iangle].total();
    test_ray_emi_ct_(0,iangle) = t[0][iangle].total();
  }

  //   _____ _                   _____                 _ _
  //  / ____| |                 |  __ \               | | |
  // | (___ | |_ ___  _ __ ___  | |__) |___  ___ _   _| | |_ ___
  //  \___ \| __/ _ \| '__/ _ \ |  _  // _ \/ __| | | | | __/ __|
  //  ____) | || (_) | | |  __/ | | \ \  __/\__ \ |_| | | |_\__ \
  // |_____/ \__\___/|_|  \___| |_|  \_\___||___/\__,_|_|\__|___/

  test_ray_emi_zn_ = calin::std_to_eigenvec(zn0);
  test_ray_emi_z_ = calin::std_to_eigenvec(s_->xknot());

  // Save the values at the bottom of the atmosphere
  test_ray_boa_x_.resize(levels.size(),zn0.size());
  test_ray_boa_ct_.resize(levels.size(),zn0.size());
  for(unsigned ilevel=0; ilevel<levels.size(); ilevel++) {
    for(unsigned iangle=0;iangle<zn0.size();iangle++) {
      test_ray_boa_x_(ilevel,iangle) = x[ilevel][iangle].total();
      test_ray_boa_ct_(ilevel,iangle) = t[ilevel][iangle].total();
    }
  }

  // Subtract off emission time and place to shift rays to origin of x & ct
  // at their emission points
  test_ray_boa_x_ -= test_ray_emi_x_;
  test_ray_boa_ct_ -= test_ray_emi_ct_;
  for(unsigned iobs=0; iobs<zobs_.size(); iobs++) {
    test_ray_obs_x_[iobs] -= test_ray_emi_x_;
    test_ray_obs_ct_[iobs] -= test_ray_emi_ct_;
  }

  // Spline 3 - vertical propagation time correction to bottom of atmosphere
  for(unsigned ilevel=0; ilevel<levels.size(); ilevel++) {
    v[ilevel] = t[ilevel][0].total() - levels[ilevel].z;
  }
  s_->add_spline(v, "vertical ct correction (boa) [cm]");

  //  ______ _ _                                _              __  __           _      _
  // |  ____(_) |       /\                     | |            |  \/  |         | |    | |
  // | |__   _| |_     /  \   _ __   __ _ _   _| | __ _ _ __  | \  / | ___   __| | ___| |
  // |  __| | | __|   / /\ \ | '_ \ / _` | | | | |/ _` | '__| | |\/| |/ _ \ / _` |/ _ \ |
  // | |    | | |_   / ____ \| | | | (_| | |_| | | (_| | |    | |  | | (_) | (_| |  __/ |
  // |_|    |_|\__| /_/    \_\_| |_|\__, |\__,_|_|\__,_|_|    |_|  |_|\___/ \__,_|\___|_|
  //                                 __/ |
  //                                |___/

  // Five splines per observarion level : (4,5,6,7,8), (9,10,11,12,13) etc...

  // 4 + iobs*5 : vertical time correction to observation level
  // 5 + iobs*5 : position (x) correction "a" coefficient - term in tan_i
  // 6 + iobs*5 : position (x) correction "b" coefficient - term in tan^3_i
  // 5 + iobs*5 : time (ct) correction "a" coefficient - term in tan_i
  // 6 + iobs*5 : time (ct) correction "b" coefficient - term in tan^3_i

  // In "normal" mode only the first two are used, as the ratio of b/a is fixed
  // and the ct and x corrections are both derived from the "a_x" spline

  for(unsigned iobs=0; iobs<zobs_.size(); iobs++) {
    // Spline 4,9,14 ... vertical propagation time correction to obs level
    for(unsigned ilevel=0; ilevel<levels.size(); ilevel++) {
      double zlev = (levels[ilevel].z-zobs_[iobs]);
      v[ilevel] = test_ray_obs_ct_[iobs](ilevel,0) - zlev;
    }

    std::vector<double> va_x(levels.size());
    std::vector<double> vb_x(levels.size());
    std::vector<double> va_ct(levels.size());
    std::vector<double> vb_ct(levels.size());

    double ab_ratio_x = 0;
    double ab_ratio_ct = 0;
    if(!high_accuracy_mode_) {
      double dzmin = std::abs(levels[0].z - config.angular_model_optimization_altitude());
      unsigned ilevel = 0;
      for(unsigned jlevel=1; jlevel<levels.size(); jlevel++) {
        double dz = std::abs(levels[jlevel].z - config.angular_model_optimization_altitude());
        if(dz < dzmin) {
          dzmin = dz;
          ilevel = jlevel;
        }
      }
      const double z_i = levels[ilevel].z;
      const double z_r = zobs_[iobs];
      const double n_i = 1.0 + this->n_minus_one(z_i);
      const double n_r_inv = 1.0 / (1.0 + this->n_minus_one(z_r));

      double a_x;
      double b_x;
      double a_ct;
      double b_ct;
      fit_refraction_angular_coefficients(z_i, z_r, n_i, n_r_inv, zn0,
        test_ray_obs_x_[iobs].row(ilevel), test_ray_obs_ct_[iobs].row(ilevel),
        a_x, b_x, a_ct, b_ct);

      ab_ratio_x = b_x/a_x;
      ab_ratio_ct = b_ct/a_ct;
    }

    for(unsigned ilevel=0; ilevel<levels.size(); ilevel++) {
      const double z_i = levels[ilevel].z;
      const double z_r = zobs_[iobs];
      const double n_i = 1.0 + this->n_minus_one(z_i);
      const double n_r_inv = 1.0 / (1.0 + this->n_minus_one(z_r));

      double a_x;
      double b_x;
      double a_ct;
      double b_ct;
      fit_refraction_angular_coefficients(z_i, z_r, n_i, n_r_inv, zn0,
        test_ray_obs_x_[iobs].row(ilevel), test_ray_obs_ct_[iobs].row(ilevel),
        a_x, b_x, a_ct, b_ct,
        !high_accuracy_mode_, ab_ratio_x, ab_ratio_ct);

      va_x[ilevel] = a_x;
      vb_x[ilevel] = b_x;
      va_ct[ilevel] = a_ct;
      vb_ct[ilevel] = b_ct;
    }

    //   _____ _                   __  __           _      _    _____            __  __
    //  / ____| |                 |  \/  |         | |    | |  / ____|          / _|/ _|
    // | (___ | |_ ___  _ __ ___  | \  / | ___   __| | ___| | | |     ___   ___| |_| |_ ___
    //  \___ \| __/ _ \| '__/ _ \ | |\/| |/ _ \ / _` |/ _ \ | | |    / _ \ / _ \  _|  _/ __|
    //  ____) | || (_) | | |  __/ | |  | | (_) | (_| |  __/ | | |___| (_) |  __/ | | | \__ \
    // |_____/ \__\___/|_|  \___| |_|  |_|\___/ \__,_|\___|_|  \_____\___/ \___|_| |_| |___/

    s_->add_spline(v, "vertical ct correction (iobs=" + std::to_string(iobs) + ") [cm]");
    s_->add_spline(va_x, "refraction a_x coeff (iobs=" + std::to_string(iobs) + ") [cm]");
    s_->add_spline(vb_x, "refraction b_x coeff (iobs=" + std::to_string(iobs) + ") [cm]");
    s_->add_spline(va_ct, "refraction a_ct coeff (iobs=" + std::to_string(iobs) + ") [cm]");
    s_->add_spline(vb_ct, "refraction b_ct coeff (iobs=" + std::to_string(iobs) + ") [cm]");
  }
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

bool LayeredRefractiveAtmosphere::
propagate_ray_with_refraction(calin::math::ray::Ray& ray, unsigned iobs)
{
  double dz = zobs_[iobs] - ray.z();
  if(ray.uz()>=0 or dz>0)return false;

  double n0;
  double dct_v;
  double dct_r;
  double dx_r;

  unsigned ispline = 4+iobs*3;
  s_->value(ray.z(), 2, n0, ispline, dct_v, ispline+1, dct_r, ispline+2, dx_r);
  n0 = std::exp(n0);

  double sec_i = 1.0/ray.uz();

  ray.propagate_dist(dz*sec_i);

  double sin2_i = SQR(ray.ux()) + SQR(ray.uy());
  double ni_over_nr = n0 * nobs_inv_[iobs];
  double sin2_r = sin2_i * SQR(ni_over_nr);
  double cos_i = ray.uz();
  double cos_r = -std::sqrt(1.0 - sin2_r);

  ray.uz() = cos_r;
  ray.ux() *= ni_over_nr;
  ray.uy() *= ni_over_nr;

  double sec2_i = SQR(sec_i);
  double sin_i = std::sqrt(sin2_i);

  // f3 = A*sin(zn1)/(cos(zn0)*cos(zn1)**2) + (1-A)*sin(zn1)**3/(cos(zn0)*cos(zn1)**4)

  double t1 = sin_i * sec2_i / cos_r;
  double t2 = t1 * sin2_i * sec2_i;

  // f3 = A*sin(zn0)/(cos(zn0)**2*cos(zn1)) + (1-A)*sin(zn0)**3/(cos(zn0)**4*cos(zn1))


  if(sin2_r > 0) {
    double sin_r_inv = 1.0/std::sqrt(sin2_r);

  }

  return true;
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
