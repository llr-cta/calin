/*

   calin/simulation/tracker.cpp -- Stephen Fegan -- 2015-07-17

   Base class for all air shower track visitors

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

#include <sys/time.h>
#include <iostream>
#include <cassert>
#include <simulation/tracker.hpp>
#include <io/log.hpp>

using namespace calin::simulation::tracker;
using namespace calin::io::log;

calin::simulation::tracker::ParticleType
calin::simulation::tracker::pdg_type_to_particle_type(int pdg_type)
{
  switch(pdg_type)
  {
    case 22:    return ParticleType::GAMMA;
    case 1:     return ParticleType::ELECTRON;
    case -1:    return ParticleType::POSITRON;
    case 13:    return ParticleType::MUON;
    case -13:   return ParticleType::ANTI_MUON;
    case 2212:  return ParticleType::PROTON;
    case -2212: return ParticleType::ANTI_PROTON;
    default:    return ParticleType::OTHER;
  };
  assert(0);
  return ParticleType::OTHER;
}

int calin::simulation::tracker::
particle_type_to_pdg_type(calin::simulation::tracker::ParticleType track_type)
{
  switch(track_type)
  {
    case ParticleType::GAMMA:       return 22;
    case ParticleType::ELECTRON:    return 1;
    case ParticleType::POSITRON:    return -1;
    case ParticleType::MUON:        return 13;
    case ParticleType::ANTI_MUON:   return -13;
    case ParticleType::PROTON:      return 2212;
    case ParticleType::ANTI_PROTON: return -2212;
    case ParticleType::OTHER:
      throw(std::runtime_error("ParticleType::OTHER has no PDG type code"));
  };
  assert(0);
  return 0;
}

double calin::simulation::tracker::
particle_type_to_mass(ParticleType track_type)
{
  switch(track_type)
  {
    case ParticleType::GAMMA:
      return 0.0;
    case ParticleType::ELECTRON:
    case ParticleType::POSITRON:
      return 0.510998928;
    case ParticleType::MUON:
    case ParticleType::ANTI_MUON:
      return 105.6583715;
    case ParticleType::PROTON:
    case ParticleType::ANTI_PROTON:
      return 938.262046;
    case ParticleType::OTHER:
      throw(std::runtime_error("ParticleType::OTHER has no mass"));
  };
  assert(0);
  return 0;
}

double calin::simulation::tracker::
particle_type_to_charge(ParticleType track_type)
{
  switch(track_type)
  {
  case ParticleType::GAMMA:
    return 0.0;
  case ParticleType::ELECTRON:
  case ParticleType::MUON:
  case ParticleType::ANTI_PROTON:
    return -1.0;
  case ParticleType::POSITRON:
  case ParticleType::ANTI_MUON:
  case ParticleType::PROTON:
    return 1.0;
  case ParticleType::OTHER:
    throw(std::runtime_error("ParticleType::OTHER has no charge"));
  };
  assert(0);
  return 0;
}

TrackVisitor::~TrackVisitor()
{
  // nothing to see here
}

void TrackVisitor::visit_event(const Event& event, bool& kill_event)
{
  // default is to do nothing
}

void TrackVisitor::visit_track(const Track& track, bool& kill_track)
{
  // default is to do nothing
}

void TrackVisitor::leave_event()
{
  // default is to do nothing
}

VectorTrackVisitor::~VectorTrackVisitor()
{
  for(auto* ivisitor : adopted_visitors_)delete ivisitor;
}

void VectorTrackVisitor::visit_event(const Event& event, bool& kill_event)
{
  for(auto* ivisitor : visitors_)
  {
    bool ikill_event = false;
    ivisitor->visit_event(event, ikill_event);
    kill_event |= ikill_event;
  }
}

void VectorTrackVisitor::visit_track(const Track& track, bool& kill_track)
{
  for(auto* ivisitor : visitors_)
  {
    bool ikill_track = false;
    ivisitor->visit_track(track, ikill_track);
    kill_track |= ikill_track;
  }
}

void VectorTrackVisitor::leave_event()
{
  for(auto* ivisitor : visitors_)ivisitor->leave_event();
}

void VectorTrackVisitor::add_visitor(TrackVisitor* visitor, bool adopt_visitor)
{
  visitors_.emplace_back(visitor);
  if(adopt_visitor)adopted_visitors_.emplace_back(visitor);
}

void VectorTrackVisitor::
add_visitor_at_front(TrackVisitor* visitor, bool adopt_visitor)
{
  visitors_.insert(visitors_.begin(), visitor);
  if(adopt_visitor)adopted_visitors_.insert(adopted_visitors_.begin(), visitor);
}

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
