/*

   calin/simulation/iact_array_tracker.cpp -- Stephen Fegan -- 2016-07-24

   Air shower track visitor to calcilate Cherenkov cone parameters.

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

   Certain portions are from code that is Copyright 2012, Stephen Fegan
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

// Certain fragments are from my EGS5 simulation code with following header:

// EGS5AtmosphericDetector.hpp - Detector class for EGS5 system that
// - implements a layered atmospheric detector
// Stephen Fegan - sfegan@llr.in2p3.fr - September 2012

#include <Eigen/Geometry>
#include <math/special.hpp>
#include <simulation/air_cherenkov_tracker.hpp>
#include <util/log.hpp>
#include <math/geometry.hpp>
#include <math/nspace.hpp>
#include <cmath>

using namespace calin::util::log;
using namespace calin::simulation::air_cherenkov_tracker;
using calin::math::special::SQR;

AirCherenkovTrackVisitor::~AirCherenkovTrackVisitor()
{
  // nothing to see here
}

void AirCherenkovTrackVisitor::set_atmosphere(calin::simulation::atmosphere::Atmosphere* atm)
{
  // default is to do nothing
}

void AirCherenkovTrackVisitor::visit_event(const Event& event, bool& kill_event)
{
  // default is to do nothing
}

void AirCherenkovTrackVisitor::
visit_cherenkov_track(const AirCherenkovTrack& track, bool& kill_track)
{
  // default is to do nothing
}

void AirCherenkovTrackVisitor::leave_event()
{
  // default is to do nothing
}

AirCherenkovParameterCalculatorTrackVisitor::
AirCherenkovParameterCalculatorTrackVisitor(AirCherenkovTrackVisitor* visitor,
    calin::simulation::atmosphere::Atmosphere* atm,
    const calin::ix::simulation::tracker::AirCherenkovParameterCalculatorTrackVisitorConfig& cfg,
    bool adopt_visitor, bool adopt_atm):
  calin::simulation::tracker::TrackVisitor(),
  visitor_(visitor), adopt_visitor_(adopt_visitor),
  atm_(atm), adopt_atm_(adopt_atm)
{
  if(cfg.enable_forced_cherenkov_angle_mode()) {
    enable_forced_cherenkov_angle_mode_ = true;
    set_forced_cherenkov_angle(cfg.forced_cherenkov_angle());
  }
  visitor->set_atmosphere(atm);
}

AirCherenkovParameterCalculatorTrackVisitor::
~AirCherenkovParameterCalculatorTrackVisitor()
{
  if(adopt_visitor_)delete visitor_;
  if(adopt_atm_)delete atm_;
}

void AirCherenkovParameterCalculatorTrackVisitor::
visit_event(const Event& event, bool& kill_event)
{
  visitor_->visit_event(event, kill_event);
}

void AirCherenkovParameterCalculatorTrackVisitor::
visit_track(const calin::simulation::tracker::Track& track, bool& kill_track)
{
  if(track.q == 0 or track.dx <= 0)return; // it's store policy: no charge, no track = no radiation

  AirCherenkovTrack cherenkov;

  cherenkov.particle_track = &track;

  cherenkov.x0             = track.x0;
  cherenkov.e0             = track.e0;
  cherenkov.t0             = track.t0;

  cherenkov.dx_hat         = track.dx_hat;
  cherenkov.dx             = track.dx;
  cherenkov.de             = track.de;
  cherenkov.dt             = track.dt;

  cherenkov.x_mid          = cherenkov.x0 + 0.5*cherenkov.dx*cherenkov.dx_hat;
  cherenkov.e_mid          = cherenkov.e0 + 0.5*cherenkov.de;
  cherenkov.t_mid          = cherenkov.t0 + 0.5*cherenkov.dt;

  atm_->cherenkov_parameters(cherenkov.x_mid(2),
    cherenkov.n, cherenkov.propagation_delay);
  if(cherenkov.n<0)return;
  cherenkov.n             += 1.0;

  // Watch out for gamma^2 slightly less than 1.0
  const double g2 = SQR(std::max(cherenkov.e_mid/track.mass,1.0)); // gamma^2
  cherenkov.gamma_sq       = g2;
  const double b2 = 1.0 - 1.0/g2;                    // beta^2
  if(enable_forced_cherenkov_angle_mode_)
    cherenkov.sin2_thetac  = forced_sin2_thetac_;
  else
    cherenkov.sin2_thetac  = 1.0 - 1.0/(b2*SQR(cherenkov.n));
  if(cherenkov.sin2_thetac <= 0.0)return;
  cherenkov.yield_density  =
    YIELD_CONST*SQR(track.q)*cherenkov.sin2_thetac*cherenkov.dx;
  cherenkov.cos_thetac     = std::sqrt(1.0 - cherenkov.sin2_thetac);
  cherenkov.sin_thetac     = std::sqrt(cherenkov.sin2_thetac);

  if(std::isnan(cherenkov.yield_density) or cherenkov.yield_density>1e6) {
    LOG(INFO) << "Excessive Cherenkov yield density: Y="
              << cherenkov.yield_density
              << "\n  x0=(" << cherenkov.particle_track->x0.transpose()
              << ")\n  x1=(" << cherenkov.particle_track->x1.transpose()
              << ") dx=" << cherenkov.dx
              << "\n  n=" << cherenkov.n
              << " q=" << track.q
              << " pdg=" << track.pdg_type
              << " (" << particle_type_to_string(track.type)
              << ") e="<< cherenkov.e_mid
              << " m=" << track.mass
              << "\n  g2=" << g2
              << " b2=" << b2
              << " sin2(tc)=" << cherenkov.sin2_thetac;
  }

  visitor_->visit_cherenkov_track(cherenkov, kill_track);
}

void AirCherenkovParameterCalculatorTrackVisitor::leave_event()
{
  visitor_->leave_event();
}

void AirCherenkovParameterCalculatorTrackVisitor::
set_forced_cherenkov_angle(double theta_c_deg)
{
  if(enable_forced_cherenkov_angle_mode_)
    forced_sin2_thetac_ = SQR(std::sin(theta_c_deg/180.0*M_PI));
  else
    throw std::runtime_error("Forced cherenkov angle mode not enabled.");
}

calin::ix::simulation::tracker::AirCherenkovParameterCalculatorTrackVisitorConfig
AirCherenkovParameterCalculatorTrackVisitor::default_config()
{
  calin::ix::simulation::tracker::AirCherenkovParameterCalculatorTrackVisitorConfig cfg;
  cfg.set_forced_cherenkov_angle(1.2);
  return cfg;
}

CherenkovPhotonVisitor::~CherenkovPhotonVisitor()
{
  // nothing to see here
}

void CherenkovPhotonVisitor::set_bandpass(double epsilon0, double bandwidth, bool do_color_photons)
{
  // nothing to see here
}

void CherenkovPhotonVisitor::visit_event(const Event& event, bool& kill_event)
{
  // nothing to see here
}

void CherenkovPhotonVisitor::
visit_cherenkov_photon(const CherenkovPhoton& cherenkov_photon)
{
  // nothing to see here
}

void CherenkovPhotonVisitor::leave_event()
{
  // nothing to see here
}

MCCherenkovPhotonGenerator::
MCCherenkovPhotonGenerator(CherenkovPhotonVisitor* visitor,
    double epsilon0, double bandwidth, bool do_color_photons,
    calin::math::rng::RNG* rng, bool adopt_visitor, bool adopt_rng):
  visitor_(visitor), adopt_visitor_(adopt_visitor),
  epsilon0_(epsilon0), bandwidth_(bandwidth), do_color_photons_(do_color_photons),
  weight_(1.0), weighted_bandwidth_(bandwidth_/weight_),
  rng_(rng ? rng : new calin::math::rng::RNG(__PRETTY_FUNCTION__)), adopt_rng_(rng ? adopt_rng : true)
{
  visitor->set_bandpass(epsilon0, bandwidth, do_color_photons);
  dX_emission_ = rng_->exponential();
}

MCCherenkovPhotonGenerator::~MCCherenkovPhotonGenerator()
{
  if(adopt_visitor_)delete visitor_;
  if(adopt_rng_)delete rng_;
}

void MCCherenkovPhotonGenerator::visit_event(const Event& event, bool& kill_event)
{
  return visitor_->visit_event(event, kill_event);
}

void MCCherenkovPhotonGenerator::leave_event()
{
  return visitor_->leave_event();
}

void MCCherenkovPhotonGenerator::
visit_cherenkov_track(const AirCherenkovTrack& cherenkov_track, bool& kill_track)
{
  CherenkovPhoton photon;
  bool no_photon_emitted = true;

  const double dX_track = cherenkov_track.yield_density * weighted_bandwidth_;

  Eigen::Matrix3d Mrot;
  Eigen::Vector3d u0;

  double dX_left = dX_track;
  while(dX_emission_ < dX_left)
  {
    dX_left -= dX_emission_;
    dX_emission_ = rng_->exponential();

    if(no_photon_emitted)
    {
      no_photon_emitted = false;
      photon.air_cherenkov_track = &cherenkov_track;
      calin::math::geometry::rotation_z_to_vec_Rzy(Mrot, cherenkov_track.dx_hat);
      photon.epsilon = 0;
      photon.weight = weight_;
    }

    double dX_frac = 1.0 - dX_left/dX_track;
    photon.x0 = cherenkov_track.x0 + cherenkov_track.dx_hat*cherenkov_track.dx*dX_frac;
    double phi = 2.0*M_PI*rng_->uniform();
    u0 << cherenkov_track.sin_thetac*std::cos(phi), cherenkov_track.sin_thetac*std::sin(phi), cherenkov_track.cos_thetac;
    photon.u0 = Mrot * u0;
    photon.t0 = cherenkov_track.t0 + cherenkov_track.dt*dX_frac;
    if(do_color_photons_)photon.epsilon = epsilon0_ + rng_->uniform()*bandwidth_;
    visitor_->visit_cherenkov_photon(photon);
  }
  dX_emission_ -= dX_left;
}

CherenkovTrackYieldLogger::CherenkovTrackYieldLogger(): AirCherenkovTrackVisitor()
{
  // nothing to see here
}

CherenkovTrackYieldLogger::~CherenkovTrackYieldLogger()
{
  // nothing to see here
}

void CherenkovTrackYieldLogger::visit_event(const Event& event, bool& kill_event)
{
  track_altitude_.clear();
  track_length_.clear();
  track_yield_.clear();
}

void CherenkovTrackYieldLogger::visit_cherenkov_track(
  const AirCherenkovTrack& cherenkov_track, bool& kill_track)
{
  track_altitude_.push_back(cherenkov_track.x_mid.z());
  track_length_.push_back(cherenkov_track.dx);
  track_yield_.push_back(cherenkov_track.yield_density);
}

CherenkovTrackYieldNSpaceVisitor::CherenkovTrackYieldNSpaceVisitor(
    const calin::ix::simulation::tracker::CherenkovTrackYieldNSpaceVisitorConfig& config):
  AirCherenkovTrackVisitor(),
  config_(config), space_(space_axes(config)), p_(space_.naxes())
{
  // nothing to see here
}

CherenkovTrackYieldNSpaceVisitor::~CherenkovTrackYieldNSpaceVisitor()
{
  // nothing to see here
}

void CherenkovTrackYieldNSpaceVisitor::clear()
{
  space_.clear();
}

void CherenkovTrackYieldNSpaceVisitor::set_atmosphere(calin::simulation::atmosphere::Atmosphere* atm)
{
  atm_ = atm;
}

void CherenkovTrackYieldNSpaceVisitor::visit_event(const Event& event, bool& kill_event)
{
  x0_ = event.x0;
  rot_ = calin::math::geometry::rotation_z_to_vec_Rzyz(-event.u0).transpose();
  uz0_inv_ = 1.0/std::abs(event.u0.z());
}

void CherenkovTrackYieldNSpaceVisitor::visit_cherenkov_track(
  const AirCherenkovTrack& cherenkov_track, bool& kill_track)
{
  double t = atm_->thickness(cherenkov_track.x_mid.z())*uz0_inv_;
  Eigen::Vector3d x = rot_ * (cherenkov_track.x_mid - x0_);
  Eigen::Vector3d u = rot_ * cherenkov_track.dx_hat;

  switch(config_.axis_variables()) {
  default:
  case calin::ix::simulation::tracker::DEPTH:
    p_ << t;
    break;
  case calin::ix::simulation::tracker::XY:
    p_ << x.x(), x.y();
    break;
  case calin::ix::simulation::tracker::UXUY:
    p_ << u.x(), u.y();
    break;
  case calin::ix::simulation::tracker::DEPTH_XY:
    p_ << t, x.x(), x.y();
    break;
  case calin::ix::simulation::tracker::DEPTH_UXUY:
    p_ << t, u.x(), u.y();
    break;
  case calin::ix::simulation::tracker::DEPTH_XY_UXUY:
    p_ << t, x.x(), x.y(), u.x(), u.y();
    break;
  }

  space_.accumulate(p_, cherenkov_track.yield_density);
}

namespace {

  calin::math::nspace::Axis make_axis(
    const calin::ix::simulation::tracker::AxisBinning& ax)
  {
    return { ax.lo_bound(), ax.hi_bound(), ax.num_bins() };
  }

} // anonymous namespace

std::vector<calin::math::nspace::Axis> CherenkovTrackYieldNSpaceVisitor::space_axes(
  const calin::ix::simulation::tracker::CherenkovTrackYieldNSpaceVisitorConfig& config)
{
  std::vector<calin::math::nspace::Axis> axes;
  switch(config.axis_variables()) {
  default:
  case calin::ix::simulation::tracker::DEPTH:
    axes.push_back(make_axis(config.depth_axis()));
    break;
  case calin::ix::simulation::tracker::XY:
    axes.push_back(make_axis(config.xy_axis()));
    axes.push_back(make_axis(config.xy_axis()));
    break;
  case calin::ix::simulation::tracker::UXUY:
    axes.push_back(make_axis(config.uxuy_axis()));
    axes.push_back(make_axis(config.uxuy_axis()));
    break;
  case calin::ix::simulation::tracker::DEPTH_XY:
    axes.push_back(make_axis(config.depth_axis()));
    axes.push_back(make_axis(config.xy_axis()));
    axes.push_back(make_axis(config.xy_axis()));
    break;
  case calin::ix::simulation::tracker::DEPTH_UXUY:
    axes.push_back(make_axis(config.depth_axis()));
    axes.push_back(make_axis(config.uxuy_axis()));
    axes.push_back(make_axis(config.uxuy_axis()));
    break;
  case calin::ix::simulation::tracker::DEPTH_XY_UXUY:
    axes.push_back(make_axis(config.depth_axis()));
    axes.push_back(make_axis(config.xy_axis()));
    axes.push_back(make_axis(config.xy_axis()));
    axes.push_back(make_axis(config.uxuy_axis()));
    axes.push_back(make_axis(config.uxuy_axis()));
    break;
  }
  return axes;
}

calin::ix::simulation::tracker::CherenkovTrackYieldNSpaceVisitorConfig
CherenkovTrackYieldNSpaceVisitor::default_config()
{
  calin::ix::simulation::tracker::CherenkovTrackYieldNSpaceVisitorConfig cfg;
  cfg.set_axis_variables(calin::ix::simulation::tracker::DEPTH_XY);

  cfg.mutable_depth_axis()->set_lo_bound(0);
  cfg.mutable_depth_axis()->set_hi_bound(1300);
  cfg.mutable_depth_axis()->set_num_bins(26);

  cfg.mutable_xy_axis()->set_lo_bound(-1050.0);
  cfg.mutable_xy_axis()->set_hi_bound(1050.0);
  cfg.mutable_xy_axis()->set_num_bins(21);

  cfg.mutable_uxuy_axis()->set_lo_bound(-2.0*M_PI/180.0);
  cfg.mutable_uxuy_axis()->set_hi_bound(2.0*M_PI/180.0);
  cfg.mutable_uxuy_axis()->set_num_bins(21);

  return cfg;
}

MultiDelegatingAirCherenkovTrackVisitor::
MultiDelegatingAirCherenkovTrackVisitor(): AirCherenkovTrackVisitor()
{
  // nothing to see here
}

MultiDelegatingAirCherenkovTrackVisitor::
~MultiDelegatingAirCherenkovTrackVisitor()
{
  for(auto* d : adopted_delegates_) {
    delete d;
  }
}

void MultiDelegatingAirCherenkovTrackVisitor::
set_atmosphere(calin::simulation::atmosphere::Atmosphere* atm)
{
  for(auto* d : delegates_) {
    d->set_atmosphere(atm);
  }
}

void MultiDelegatingAirCherenkovTrackVisitor::
visit_event(const Event& event, bool& kill_event)
{
  for(auto* d : delegates_) {
    bool d_kill_event = false;
    d->visit_event(event, d_kill_event);
    kill_event = kill_event or d_kill_event;
  }
}

void MultiDelegatingAirCherenkovTrackVisitor::
visit_cherenkov_track(const AirCherenkovTrack& cherenkov_track, bool& kill_track)
{
  for(auto* d : delegates_) {
    bool d_kill_track = false;
    d->visit_cherenkov_track(cherenkov_track, d_kill_track);
    kill_track = kill_track or d_kill_track;
  }
}

void MultiDelegatingAirCherenkovTrackVisitor::leave_event()
{
  for(auto* d : delegates_) {
    d->leave_event();
  }
}

void MultiDelegatingAirCherenkovTrackVisitor::
add_delegate(AirCherenkovTrackVisitor* delegate, bool adopt_delegate)
{
  delegates_.push_back(delegate);
  if(adopt_delegate) {
    adopted_delegates_.push_back(delegate);
  }
}
