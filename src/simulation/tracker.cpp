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

#include <simulation/tracker.hpp>

using namespace calin::simulation::tracker;

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
  Track subtrack { track };
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

      visitor_->visit_track(subtrack, kill_track);
      return;
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
