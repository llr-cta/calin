/*

   calin/simulation/bfield_track_generator.hpp -- Stephen Fegan -- 2017-06-16

   Class to genereate curved tracks in a magnetic field (neglecting
   synchrotron emission)

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
#include<cassert>

#include"simulation/tracker.hpp"

#include"calin_global_definitions.hpp"

// Units:
// Height, z:    cm
// Energy, e:    MeV

namespace calin { namespace simulation { namespace bfield_track_generator {

enum PropagationMode { FWD_TO_GROUND, BWD_FIXED_DISTANCE };

class BFieldTrackGenerator
{
 public:
  BFieldTrackGenerator(calin::simulation::tracker::TrackVisitor* visitor,
    const Eigen::Vector3d& bfield_nT, double zground_or_dist, double step_size,
    PropagationMode propagation_mode = FWD_TO_GROUND,
    bool adopt_visitor = false);
  BFieldTrackGenerator(calin::simulation::tracker::TrackVisitor* visitor,
    double zground_or_dist, double step_size,
    PropagationMode propagation_mode = FWD_TO_GROUND,
    bool adopt_visitor = false);
  virtual ~BFieldTrackGenerator();

  void generate_showers(unsigned num_events,
                        calin::simulation::tracker::ParticleType type,
                        double total_energy,
                        const Eigen::Vector3d& x0 = Eigen::Vector3d(0,0,0),
                        const Eigen::Vector3d& u0 = Eigen::Vector3d(0,0,-1),
                        double weight=1.0);

 protected:
  calin::simulation::tracker::TrackVisitor* visitor_ = nullptr;
  bool adopt_visitor_ = false;
  Eigen::Vector3d bfield_;
  PropagationMode propagation_mode_ = FWD_TO_GROUND;
  double zground_or_dist_ = 0.0;
  double step_size_ = 1.0;
  int event_id_ = 0;
};

} } } // namespace calin::simulation::straight_track_generator
