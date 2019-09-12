/*

   calin/simulation/straight_track_generator.hpp -- Stephen Fegan -- 2016-07-26

   Class to genereate straight tracks

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#pragma once

#include<vector>
#include<cassert>

#include"simulation/tracker.hpp"

#include"calin_global_definitions.hpp"

// Units:
// Height, z:    cm
// Energy, e:    MeV

namespace calin { namespace simulation { namespace straight_track_generator {

class StraightTrackGenerator
{
 public:
  StraightTrackGenerator(calin::simulation::tracker::TrackVisitor* visitor,
                        double zground, bool adopt_visitor = false);
  virtual ~StraightTrackGenerator();

  void generate_showers(unsigned num_events,
                        calin::simulation::tracker::ParticleType type,
                        double total_energy,
                        const Eigen::Vector3d& x0 = Eigen::Vector3d(0,0,0),
                        const Eigen::Vector3d& u0 = Eigen::Vector3d(0,0,-1),
                        double weight=1.0);

 protected:
  calin::simulation::tracker::TrackVisitor* visitor_ = nullptr;
  bool adopt_visitor_ = false;
  double zground_ = 0;
  int event_id_ = 0;
};

} } } // namespace calin::simulation::straight_track_generator
