/*

   calin/simulation/shower_movie.hpp -- Stephen Fegan -- 2017-10-13

   Shower movie maker visitor

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#pragma once

#include<vector>
#include<map>

#include"Eigen/Core"
#include"simulation/tracker.hpp"
#include"simulation/air_cherenkov_tracker.hpp"

namespace calin { namespace simulation { namespace shower_movie {

class ShowerMovieProducerTrackVisitor: public TrackVisitor
{
public:
  ShowerMovieProducerTrackVisitor(
    calin::simulation::atmosphere::Atmosphere* atm,
    const calin::ix::simulation::tracker::AirCherenkovParameterCalculatorTrackVisitorConfig& cfg = calin::ix::simulation::tracker::AirCherenkovParameterCalculatorTrackVisitor::default_config(),
    bool adopt_atm = false);
  virtual ~ShowerMovieProducerTrackVisitor();
  virtual void visit_event(const Event& event, bool& kill_event);
  virtual void visit_track(const Track& track, bool& kill_track);
  virtual void leave_event();
private:

};


} } } // namespace calin::simulation::shower_movie
