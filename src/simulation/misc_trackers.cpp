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
#include <simulation/misc_trackers.hpp>
#include <io/log.hpp>

using namespace calin::simulation::tracker;
using namespace calin::simulation::misc_trackers;
using namespace calin::io::log;

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

void RecordingTrackVisitor::
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
}
