/*

   calin/simulation/tracker.cpp -- Stephen Fegan -- 2015-07-17

   Base class for all air shower track visitors

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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
    case 11:    return ParticleType::ELECTRON;
    case -11:   return ParticleType::POSITRON;
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
    case ParticleType::ELECTRON:    return 11;
    case ParticleType::POSITRON:    return -11;
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

MultiDelegatingTrackVisitor::~MultiDelegatingTrackVisitor()
{
  for(auto* ivisitor : adopted_visitors_)delete ivisitor;
}

void MultiDelegatingTrackVisitor::visit_event(const Event& event, bool& kill_event)
{
  for(auto* ivisitor : visitors_)
  {
    bool ikill_event = false;
    ivisitor->visit_event(event, ikill_event);
    kill_event |= ikill_event;
  }
}

void MultiDelegatingTrackVisitor::visit_track(const Track& track, bool& kill_track)
{
  for(auto* ivisitor : visitors_)
  {
    bool ikill_track = false;
    ivisitor->visit_track(track, ikill_track);
    kill_track |= ikill_track;
  }
}

void MultiDelegatingTrackVisitor::leave_event()
{
  for(auto* ivisitor : visitors_)ivisitor->leave_event();
}

void MultiDelegatingTrackVisitor::add_visitor(TrackVisitor* visitor, bool adopt_visitor)
{
  visitors_.emplace_back(visitor);
  if(adopt_visitor)adopted_visitors_.emplace_back(visitor);
}

void MultiDelegatingTrackVisitor::
add_visitor_at_front(TrackVisitor* visitor, bool adopt_visitor)
{
  visitors_.insert(visitors_.begin(), visitor);
  if(adopt_visitor)adopted_visitors_.insert(adopted_visitors_.begin(), visitor);
}
