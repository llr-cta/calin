/*

   calin/simulation/iact_array_tracker.cpp -- Stephen Fegan -- 2016-07-24

   Base class for all air shower track visitors

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <simulation/air_cherenkov_tracker.hpp>

using namespace calin::simulation::air_cherenkov_tracker;

AirCherenkovTrackVisitor::~AirCherenkovTrackVisitor()
{
  // nothing to see here
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
