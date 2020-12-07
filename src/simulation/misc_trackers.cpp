/*

   calin/simulation/tracker.cpp -- Stephen Fegan -- 2015-07-17

   Base class for all air shower track visitors

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <sys/time.h>
#include <iostream>
#include <cassert>
#include <simulation/tracker.hpp>
#include <simulation/misc_trackers.hpp>
#include <util/log.hpp>
#include <math/brent.hpp>
#include <math/constants.hpp>

using namespace calin::simulation::air_cherenkov_tracker;
using namespace calin::simulation::tracker;
using namespace calin::simulation::misc_trackers;
using namespace calin::util::log;

LengthLimitingTrackVisitor::LengthLimitingTrackVisitor(TrackVisitor* visitor,
    double dx_max, double z_max, bool adopt_visitor): TrackVisitor(),
  visitor_(visitor), adopt_visitor_(adopt_visitor),
  dx_max_(dx_max), z_max_(z_max)
{
  // nothing to see here
}

LengthLimitingTrackVisitor::~LengthLimitingTrackVisitor()
{
  if(adopt_visitor_)delete visitor_;
}

void LengthLimitingTrackVisitor::
visit_event(const Event& event, bool& kill_event)
{
  visitor_->visit_event(event, kill_event);
}

void LengthLimitingTrackVisitor::
visit_track(const Track& track, bool& kill_track)
{
  if(track.dx <= dx_max_ or track.x1[2] > z_max_)
    return visitor_->visit_track(track, kill_track);

  // Otherwise ...
  Track subtrack = track;

  double dbl_n = std::ceil(track.dx / dx_max_);
  unsigned n_subtrack = dbl_n;
  double norm = 1.0/dbl_n;

  subtrack.dx = track.dx * norm;
  subtrack.de = track.de * norm;
  subtrack.dt = track.dt * norm;

  Eigen::Vector3d du = (track.u1 - track.u0) * norm;

  subtrack.x0 = track.x0;
  subtrack.u0 = track.u0;
  subtrack.e0 = track.e0;
  subtrack.t0 = track.t0;

  for(unsigned i_subtrack=1; i_subtrack<n_subtrack; i_subtrack++)
  {
    subtrack.x1 = track.x0 + (i_subtrack*subtrack.dx) * subtrack.dx_hat;
    subtrack.u1 = track.u0 + i_subtrack*du;
    subtrack.u1.normalize();
    subtrack.e1 = track.e0 + i_subtrack*subtrack.de;
    subtrack.t1 = track.t0 + i_subtrack*subtrack.dt;

    visitor_->visit_track(subtrack, kill_track);
    if(kill_track)return;

    subtrack.x0 = subtrack.x1;
    subtrack.u0 = subtrack.u1;
    subtrack.e0 = subtrack.e1;
    subtrack.t0 = subtrack.t1;
  }

  subtrack.x1 = track.x1;
  subtrack.u1 = track.u1;
  subtrack.e1 = track.e1;
  subtrack.t1 = track.t1;
  return visitor_->visit_track(subtrack, kill_track);
}

void LengthLimitingTrackVisitor::leave_event()
{
  visitor_->leave_event();
}

TimeLimitingTrackVisitor::
TimeLimitingTrackVisitor(double t_max_sec): TrackVisitor(),
  t_event_max_(t_max_sec * 1000000000.0)
{
  // nothing to see here
}

TimeLimitingTrackVisitor::~TimeLimitingTrackVisitor()
{
  // nothing to see here
}

void TimeLimitingTrackVisitor::
visit_event(const Event& event, bool& kill_event)
{
  t_event_start_ = get_time();
  do_kill_track_ = false;
  num_track_since_check_ = 0;
}

void TimeLimitingTrackVisitor::
visit_track(const Track& track, bool& kill_track)
{
restart:
  if(do_kill_track_)
  {
    LOG(INFO)
      << int(track.type) << ' '
      << track.pdg_type << ' '
      << track.q << ' '
      << track.mass << ' '
      << "[ " << track.x0.transpose() << " ] "
      << "[ " << track.u0.transpose() << " ] "
      << track.e0 << ' '
      << track.t0 << ' '
      << "[ " << track.dx_hat.transpose() << " ] "
      << track.dx << ' '
      << track.de << ' '
      << track.dt << '\n';
    kill_track = true;
    return;
  }

  if(num_track_since_check_ >= num_track_before_check_)
  {
    uint64_t t_now = get_time();
    if(t_now - t_event_start_ > t_event_max_) {
      do_kill_track_ = true;
      LOG(WARNING) << "Event time limit exceeded, remaining tracks will be printed.";
      goto restart;
    }
    num_track_since_check_ = 0;
  }
  else
  {
    ++num_track_since_check_;
  }
}

uint64_t TimeLimitingTrackVisitor::get_time()
{
  struct timeval tv;
  ::gettimeofday(&tv, nullptr);
  uint64_t t = uint64_t(tv.tv_sec)*1000000000 + uint64_t(tv.tv_usec)*1000;
  return t;
};

DebugStatsTrackVisitor::DebugStatsTrackVisitor()
{
  // nothing to see here
}

DebugStatsTrackVisitor::~DebugStatsTrackVisitor()
{
  // nothing to see here
}

void DebugStatsTrackVisitor::visit_event(const Event& event, bool& kill_event)
{
  stats_.clear();
}

void DebugStatsTrackVisitor::visit_track(const Track& track, bool& kill_track)
{
  TrackStats& stats = stats_[track.pdg_type];
  if(track.dx == 0)++stats.track_length_zero;
  else stats.track_length_hist.insert(std::log10(track.dx));
}

std::vector<int> DebugStatsTrackVisitor::particle_pdg_types() const
{
  std::vector<int> pdg_types;
  for(const auto& istat : stats_)pdg_types.emplace_back(istat.first);
  return pdg_types;
}

unsigned DebugStatsTrackVisitor::track_length_zero(int pdg_type) const
{
  return find_stats(pdg_type).track_length_zero;
}

const calin::math::histogram::SimpleHist&
DebugStatsTrackVisitor::track_length_log_hist(int pdg_type) const
{
  return find_stats(pdg_type).track_length_hist;
}

const DebugStatsTrackVisitor::TrackStats&
DebugStatsTrackVisitor::find_stats(int pdg_type) const
{
  auto ifind = stats_.find(pdg_type);
  if(ifind == stats_.end())
    throw std::out_of_range("Statistics for particle type not found: " +
      std::to_string(pdg_type));
  return ifind->second;
}

GroundInterceptTrackVisitor::
GroundInterceptTrackVisitor(double ground_level_cm):
  TrackVisitor(), ground_level_cm_(ground_level_cm)
{
  // nothing to see here
}

GroundInterceptTrackVisitor::~GroundInterceptTrackVisitor()
{
  // nothing to see here
}

void GroundInterceptTrackVisitor::
add_particle_type_filter(calin::simulation::tracker::ParticleType pt)
{
  type_filter_.insert(pt);
}

void GroundInterceptTrackVisitor::
visit_event(const Event& event, bool& kill_event)
{
  ground_tracks_.clear();
}

void GroundInterceptTrackVisitor::
visit_track(const Track& track, bool& kill_track)
{
  //LOG(INFO) << track.x0.z() << ' ' <<
  if(track.x0.z() > ground_level_cm_ and track.x1.z() <= ground_level_cm_ and
      (type_filter_.empty() or type_filter_.find(track.type)!=type_filter_.end()))
    ground_tracks_.push_back(track);
}

// =============================================================================
// =============================================================================
//
// RecordingTrackVisitor
//
// =============================================================================
// =============================================================================

RecordingTrackVisitor::RecordingTrackVisitor(): TrackVisitor()
{
  // nothing to see here
}

RecordingTrackVisitor::~RecordingTrackVisitor()
{
  // nothing to see here
}

void RecordingTrackVisitor::
add_particle_type_filter(calin::simulation::tracker::ParticleType pt)
{
  type_filter_.insert(pt);
}

void RecordingTrackVisitor::
visit_event(const Event& event, bool& kill_event)
{
  tracks_.clear();
  event_ = event;
}

void RecordingTrackVisitor::
visit_track(const Track& track, bool& kill_track)
{
  if(type_filter_.empty() or type_filter_.find(track.type)!=type_filter_.end())
    tracks_.push_back(track);
}

bool RecordingTrackVisitor::
replay_event(calin::simulation::tracker::TrackVisitor* visitor) const
{
  bool kill_event = false;
  visitor->visit_event(event_, kill_event);
  if(!kill_event) {
    for(const auto& track : tracks_) {
      bool kill_track = false;
      visitor->visit_track(track, kill_track);
    }
  }
  visitor->leave_event();
  return kill_event;
}

// =============================================================================
// =============================================================================
//
// SubshowerTrackVisitor
//
// =============================================================================
// =============================================================================

SubshowerTrackVisitor::SubshowerTrackVisitor(double subshower_energy_mev): TrackVisitor()
{
  // nothing to see here
}

SubshowerTrackVisitor::~SubshowerTrackVisitor()
{
  // nothing to see here
}

void SubshowerTrackVisitor::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
  if(not replay_) {
    trunk_track_visitor_.visit_event(event, kill_event);
    subshowers_.clear();
  }
}

void SubshowerTrackVisitor::
visit_track(const calin::simulation::tracker::Track& track, bool& kill_track)
{
  if(replay_) {
    Track subshower_track = track;
    subshower_track.t0 += current_subshower_->t0;
    subshower_track.t1 += current_subshower_->t0;
    subshower_track.weight *= current_subshower_->weight;
    replay_visitor_->visit_track(subshower_track, kill_track);
  } else if(track.type==calin::simulation::tracker::ParticleType::GAMMA and
      track.e0<subshower_energy_mev_) {
    Event subshower_event;
    subshower_event.type = track.type;
    subshower_event.x0 = track.x0;
    subshower_event.u0 = track.u0;
    subshower_event.e0 = track.e0;
    subshower_event.t0 = track.t0;
    subshower_event.weight = track.weight;
    subshowers_.push_back(subshower_event);
    kill_track = true;
  } else {
    trunk_track_visitor_.visit_track(track, kill_track);
  }
}

void SubshowerTrackVisitor::leave_event()
{
  if(not replay_) {
    trunk_track_visitor_.leave_event();
  }
}

calin::simulation::tracker::Event SubshowerTrackVisitor::event() const
{
  return trunk_track_visitor_.event();
}

void SubshowerTrackVisitor::clear()
{
  trunk_track_visitor_.clear_tracks();
  subshowers_.clear();
}

bool SubshowerTrackVisitor::generate_shower(
  calin::simulation::tracker::ShowerGenerator* generator,
  calin::simulation::tracker::TrackVisitor* visitor)
{
  replay_ = true;
  bool kill_event = trunk_track_visitor_.replay_event(visitor);
  if(not kill_event) {
    for(const auto& sub_event : subshowers_) {
      current_subshower_ = &sub_event;
      generator->generate_showers(visitor, 1, sub_event.type, sub_event.e0,
        sub_event.x0, sub_event.u0);
    }
  }
  return kill_event;
}

// =============================================================================
// =============================================================================
//
// ShowerMovieProducerTrackVisitor
//
// =============================================================================
// =============================================================================

ShowerMovieProducerTrackVisitor::
ShowerMovieProducerTrackVisitor(calin::simulation::atmosphere::Atmosphere* atm,
  double zground, calin::simulation::detector_efficiency::AtmosphericAbsorption* atm_abs,
  const calin::ix::simulation::tracker::ShowerMovieProducerTrackVisitorConfig& config,
  bool adopt_atm, bool adopt_atm_abs): TrackVisitor(), config_(config), zground_(zground)
{
  if(config_.disable_cherenkov_light()) {
    if(atm and adopt_atm)delete atm;
  } else if(atm) {
    auto* rng = new calin::math::rng::RNG(__PRETTY_FUNCTION__);
    cherenkov_ = new calin::simulation::air_cherenkov_tracker::AirCherenkovParameterCalculatorTrackVisitor(
      new MCCherenkovPhotonGenerator(new ShowerMovieProducerCherenkovPhotonVisitor(
          /* parent = */ this, atm, rng, atm_abs, adopt_atm_abs),
        config_.cherenkov_epsilon0(),
        config_.cherenkov_bandwidth() * config_.cherenkov_yield_factor(),
        /* do_color_photons = */ false,
        /* rng = */ rng,
        /* adopt_visitor = */ true,
        /* adopt_rng = */ true),
      atm, config_.cherenkov_params(),
      /* adopt_visitor = */ true,
      /* adopt_atm = */ adopt_atm);
  }
}

ShowerMovieProducerTrackVisitor::~ShowerMovieProducerTrackVisitor()
{
  delete cherenkov_;
}

void ShowerMovieProducerTrackVisitor::
visit_event(const Event& event, bool& kill_event)
{
  frames_.clear();
  if(cherenkov_)cherenkov_->visit_event(event, kill_event);
}

void ShowerMovieProducerTrackVisitor::
visit_track(const Track& track, bool& kill_track)
{
  const double t_exposure = (config_.frame_exposure_time() > 0) ?
    config_.frame_exposure_time() : config_.frame_advance_time();
  const double t0 = track.t0;
  const double t1 = (config_.max_time() > 0) ?
    std::min(track.t0+track.dt,config_.max_time()) : track.t0+track.dt;

  double t=track.t0;
  int it = int(std::floor(t/config_.frame_advance_time()));
  while(it*config_.frame_advance_time()+t_exposure > t)--it;
  ++it;
  while(t<t1)
  {
    Frame& frame = frames_[it];
    double t_it = it*config_.frame_advance_time();
    double tseg0 = std::max(t,t_it);
    double tseg1 = std::min(t_it+t_exposure, t1);
    Eigen::Vector3d x0 = track.x0 + (tseg0-t0)/track.dt * track.dx * track.dx_hat;
    Eigen::Vector3d x1 = track.x0 + (tseg1-t0)/track.dt * track.dx * track.dx_hat;
    if(x0.z() < zground_)goto next_track;
    if(x1.z() < zground_)x1 -= (x1.z() - zground_)/track.dx_hat.z() * track.dx_hat;
    switch(track.type) {
    case calin::simulation::tracker::ParticleType::GAMMA:
      frame.gamma.emplace_back(x0,x1);
      break;
    case calin::simulation::tracker::ParticleType::ELECTRON:
    case calin::simulation::tracker::ParticleType::POSITRON:
      frame.electron.emplace_back(x0,x1);
      break;
    case calin::simulation::tracker::ParticleType::MUON:
    case calin::simulation::tracker::ParticleType::ANTI_MUON:
      frame.muon.emplace_back(x0,x1);
      break;
    case calin::simulation::tracker::ParticleType::PROTON:
    case calin::simulation::tracker::ParticleType::ANTI_PROTON:
    case calin::simulation::tracker::ParticleType::OTHER:
      frame.other.emplace_back(x0,x1);
      break;
    }
  next_track:
    ++it;
    t = std::max(track.t0, it*config_.frame_advance_time());
  }
  if(cherenkov_)cherenkov_->visit_track(track, kill_track);
}

void ShowerMovieProducerTrackVisitor::leave_event()
{
  if(cherenkov_)cherenkov_->leave_event();
}

calin::ix::simulation::tracker::ShowerMovieProducerTrackVisitorConfig
ShowerMovieProducerTrackVisitor::default_config()
{
  calin::ix::simulation::tracker::ShowerMovieProducerTrackVisitorConfig cfg;
  cfg.set_frame_advance_time(3);
  cfg.set_frame_exposure_time(0);
  cfg.set_max_time(300 * 1e6);
  *cfg.mutable_cherenkov_params() =
    calin::simulation::air_cherenkov_tracker::AirCherenkovParameterCalculatorTrackVisitor::default_config();
  cfg.set_cherenkov_epsilon0(1.5);
  cfg.set_cherenkov_bandwidth(3.0);
  cfg.set_cherenkov_yield_factor(1.0);
  return cfg;
}

ShowerMovieProducerTrackVisitor::ShowerMovieProducerCherenkovPhotonVisitor::
ShowerMovieProducerCherenkovPhotonVisitor(ShowerMovieProducerTrackVisitor* parent,
    calin::simulation::atmosphere::Atmosphere* atm, calin::math::rng::RNG* rng,
    calin::simulation::detector_efficiency::AtmosphericAbsorption* atm_abs,
    bool adopt_atm_abs):
  calin::simulation::air_cherenkov_tracker::CherenkovPhotonVisitor(),
  parent_(parent), atm_(atm), rng_(rng),
  atm_abs_(atm_abs), adopt_atm_abs_(adopt_atm_abs)
{
  // nothing to see here
}

ShowerMovieProducerTrackVisitor::ShowerMovieProducerCherenkovPhotonVisitor::
~ShowerMovieProducerCherenkovPhotonVisitor()
{
  if(adopt_atm_abs_)delete atm_abs_;
}

void ShowerMovieProducerTrackVisitor::ShowerMovieProducerCherenkovPhotonVisitor::
visit_cherenkov_photon(const calin::simulation::air_cherenkov_tracker::CherenkovPhoton& cherenkov_photon)
{
#if 0
  static unsigned nprint=100;
  if(nprint) {
    LOG(INFO) << cherenkov_photon.x0.transpose() << ' '
      << cherenkov_photon.t0 << ' '
      << cherenkov_photon.u0.transpose();
    --nprint;
  }
#endif

  const calin::ix::simulation::tracker::ShowerMovieProducerTrackVisitorConfig& config =
    parent_->config_;
  const double t_exposure = (config.frame_exposure_time() > 0) ?
    config.frame_exposure_time() : config.frame_advance_time();

  Eigen::Vector3d x = cherenkov_photon.x0;
  double t = cherenkov_photon.t0;
  double z1 = parent_->zground_;
  if(atm_abs_ and cherenkov_photon.u0.z()<0) {
    double w = std::abs(cherenkov_photon.u0.z());
    double prob = rng_->uniform();
    double survival_prob = bw_.bandwidth(x.z(), w);
    if(prob > survival_prob) {
      z1 = calin::math::brent::brent_zero(z1, x.z(),
        [this,w,prob](double z){ return prob - bw_.bandwidth(z, w); },
        prob-1.0,prob-survival_prob);
    }

#if 0
    static unsigned nprint=100;
    if(nprint) {
      LOG(INFO) << x.transpose() << ' ' << w << ' '
        << prob << ' ' << survival_prob << ' ' << parent_->zground_ << ' '
        << z1;
      --nprint;
    }
#endif
  }

  if(cherenkov_photon.u0.z()>=0 and config.max_time()<=0)
    return; // upward going photon will travel forever, so don't track it

  int it = int(std::floor(t/config.frame_advance_time()));
  while(it*config.frame_advance_time()+t_exposure > t)--it;
  ++it;
  while(x.z()>z1 and (config.max_time()<=0 or t<config.max_time()))
  {
    double n = atm_->n_minus_one(x.z()) + 1.0;
    double t_it = it*config.frame_advance_time();
    double t1 = t_it+t_exposure;
    Eigen::Vector3d x1 = x + cherenkov_photon.u0*(t1-t)*calin::math::constants::cgs_c*1e-9/n;
    if(x1.z() < z1)x1 -= (x1.z() - z1)/cherenkov_photon.u0.z() * cherenkov_photon.u0;
    Frame& frame = parent_->frames_[it];
    frame.cherenkov.emplace_back(x,x1);
    ++it;
    double t_new = std::max(cherenkov_photon.t0, it*config.frame_advance_time());
    double dt = t_new - t;
    x += cherenkov_photon.u0*dt*calin::math::constants::cgs_c*1e-9/n;
    t = t_new;
  }
}

void ShowerMovieProducerTrackVisitor::ShowerMovieProducerCherenkovPhotonVisitor::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
  const calin::ix::simulation::tracker::ShowerMovieProducerTrackVisitorConfig& config =
    parent_->config_;
  const double deps = 0.00001;
  const double eff0 = 1.0/config.cherenkov_bandwidth();
  if(atm_abs_)
  {
    calin::simulation::detector_efficiency::DetectionEfficiency eff(eff0);
    eff.insert(config.cherenkov_epsilon0() - 2*deps, 0.0);
    eff.insert(config.cherenkov_epsilon0() - deps, 0.0);
    eff.insert(config.cherenkov_epsilon0() + deps, eff0);
    for(auto epsilon : atm_abs_->energy_ev())
      if(epsilon > config.cherenkov_epsilon0() + deps and
          epsilon < config.cherenkov_epsilon0() + config.cherenkov_bandwidth() - deps)
        eff.insert(epsilon, eff0);
    eff.insert(config.cherenkov_epsilon0() + config.cherenkov_bandwidth() - deps, eff0);
    eff.insert(config.cherenkov_epsilon0() + config.cherenkov_bandwidth() + deps, 0.0);
    eff.insert(config.cherenkov_epsilon0() + config.cherenkov_bandwidth() + 2*deps, 0.0);
    bw_ = atm_abs_->integrateBandwidth(parent_->zground_, std::abs(event.u0.z()), eff);
  }
}
