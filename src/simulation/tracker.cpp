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

#include <cassert>
#include <simulation/tracker.hpp>

using namespace calin::simulation::tracker;

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


LengthLimitingTrackVisitor::LengthLimitingTrackVisitor(TrackVisitor* visitor,
    double dx_max, bool adopt_visitor): TrackVisitor(),
  visitor_(visitor), adopt_visitor_(adopt_visitor), dx_max_(dx_max)
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
  if(track.dx <= dx_max_)return visitor_->visit_track(track, kill_track);

  // Otherwise ...
  Track subtrack = track;
  double dx_sum = 0;

  subtrack.dx = dx_max_;
  subtrack.de = track.de * dx_max_/track.dx;
  subtrack.dt = track.dt * dx_max_/track.dx;

  Eigen::Vector3d du = (track.u1 - track.u0) * dx_max_/track.dx;

  while(dx_sum<track.dx)
  {
    dx_sum += subtrack.dx;
    if(dx_sum >= track.dx)
    {
      subtrack.x1 = track.x1;
      subtrack.u1 = track.u1;
      subtrack.e1 = track.e1;
      subtrack.t1 = track.t1;
      subtrack.dx = (subtrack.x1 - subtrack.x0).norm();
      subtrack.de = subtrack.e1 - subtrack.e0;
      subtrack.dt = subtrack.t1 - subtrack.t0;

      return visitor_->visit_track(subtrack, kill_track);
    }
    else
    {
      subtrack.x1 = subtrack.x0 + subtrack.dx * subtrack.dx_hat;
      subtrack.u1 = subtrack.u0 + du;
      subtrack.u1.normalize();
      subtrack.e1 = subtrack.e0 + subtrack.de;
      subtrack.t1 = subtrack.t0 + subtrack.dt;

      visitor_->visit_track(subtrack, kill_track);
      if(kill_track)return;

      subtrack.x0 = subtrack.x1;
      subtrack.u0 = subtrack.u1;
      subtrack.e0 = subtrack.e1;
      subtrack.t0 = subtrack.t1;
    }
  }
  assert(0);
}

void LengthLimitingTrackVisitor::leave_event()
{
  visitor_->leave_event();
}
