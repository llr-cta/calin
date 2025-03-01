/*

   calin/simulation/corsika8_shower_generator.hpp -- Stephen Fegan -- 2024-08-23

   Class to generate extensive air showers using Corsika 8

   Copyright 2024, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include"simulation/atmosphere.hpp"
#include"simulation/tracker.hpp"
#include"simulation/world_magnetic_model.hpp"
#include"simulation/corsika8_shower_generator.pb.h"

#include"calin_global_definitions.hpp"

// Units:
// Height, z:    cm
// Energy, e:    MeV

namespace calin { namespace simulation { namespace corsika8_shower_generator {

  class CORSIKA8ShowerGenerator: public calin::simulation::tracker::ShowerGenerator
  {
  public:
    CALIN_TYPEALIAS(config_type,
      calin::ix::simulation::corsika8_shower_generator::CORSIKA8ShowerGeneratorConfiguration);

    virtual ~CORSIKA8ShowerGenerator();

    void generate_showers(calin::simulation::tracker::TrackVisitor* visitor,
                          unsigned num_events,
                          calin::simulation::tracker::ParticleType type,
                          double total_energy,
                          const Eigen::Vector3d& x0 = Eigen::Vector3d(0,0,0),
                          const Eigen::Vector3d& u0 = Eigen::Vector3d(0,0,-1),
                          double weight=1.0) override = 0;

    static config_type default_config();

    static CORSIKA8ShowerGenerator* new_instance(const config_type& config);
  };

} } } // namespace calin::simulation::corsika8_shower_generator
