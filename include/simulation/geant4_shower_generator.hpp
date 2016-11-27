/*

   calin/simulation/geant4_shower_generator.hpp -- Stephen Fegan -- 2015-06-23

   Class to generate extensive air showers using Geant-4

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

#pragma once

#include<vector>
#include<cassert>

#include"simulation/atmosphere.hpp"
#include"simulation/tracker.hpp"
#include"simulation/world_magnetic_model.hpp"

#include"calin_global_definitions.hpp"

// Units:
// Height, z:    cm
// Energy, e:    MeV

#ifndef SWIG
class G4UImanager;
class G4UIsession;
class G4RunManager;
class G4VExceptionHandler;
#endif

namespace calin { namespace simulation { namespace geant4_shower_generator {

#ifndef SWIG
class EAS_SteppingAction;
class EAS_PrimaryGeneratorAction;
class EAS_StackingAction;
class EAS_UserEventAction;
#endif

enum class VerbosityLevel {
  SUPPRESSED_ALL, SUPRESSED_STDOUT, NORMAL, VERBOSE_EVENT, VERBOSE_TRACKING,
    VERBOSE_EVERYTHING };

class Geant4ShowerGenerator
{
 public:
  Geant4ShowerGenerator(calin::simulation::tracker::TrackVisitor* visitor,
                        calin::simulation::atmosphere::Atmosphere* atm,
                        unsigned num_atm_layers, double zground, double ztop,
                        calin::simulation::world_magnetic_model::FieldVsElevation* bfield = nullptr,
                        VerbosityLevel verbose_level =
                        VerbosityLevel::SUPRESSED_STDOUT, uint32_t seed = 0,
                        bool adopt_visitor = false,
                        bool adopt_atm = false,
                        bool adopt_bfield = false);
  Geant4ShowerGenerator(calin::simulation::tracker::TrackVisitor* visitor,
                        calin::simulation::atmosphere::Atmosphere* atm,
                        unsigned num_atm_layers, double zground, double ztop,
                        VerbosityLevel verbose_level,
                        calin::simulation::world_magnetic_model::FieldVsElevation* bfield = nullptr,
                        uint32_t seed = 0,
                        bool adopt_visitor = false,
                        bool adopt_atm = false,
                        bool adopt_bfield = false):
    Geant4ShowerGenerator(visitor, atm, num_atm_layers, zground, ztop, bfield,
      verbose_level, seed, adopt_visitor, adopt_atm, adopt_bfield) { }
  virtual ~Geant4ShowerGenerator();

  void set_minimum_energy_cut(double emin_mev);
  void generate_showers(unsigned num_events,
                        calin::simulation::tracker::ParticleType type,
                        double total_energy,
                        const Eigen::Vector3d& x0 = Eigen::Vector3d(0,0,0),
                        const Eigen::Vector3d& u0 = Eigen::Vector3d(0,0,-1),
                        double weight=1.0);

  uint32_t random_seed() const { return seed_; }

 protected:
  calin::simulation::tracker::TrackVisitor* visitor_ = nullptr;
  bool adopt_visitor_ = false;
  calin::simulation::atmosphere::Atmosphere* atm_ = nullptr;
  bool adopt_atm_ = false;
  double ztop_of_atm_ = 0;
  double zground_ = 0;
  calin::simulation::world_magnetic_model::FieldVsElevation* bfield_ = nullptr;
  bool adopt_bfield_ = false;
  uint32_t seed_ = 0;

  G4UImanager* ui_manager_                = nullptr;    // don't delete me
  G4UIsession* ui_session_                = nullptr;    // delete me
  G4RunManager* run_manager_              = nullptr;    // delete me
  EAS_SteppingAction* step_action_        = nullptr;    // don't delete me
  EAS_PrimaryGeneratorAction* gen_action_ = nullptr;    // don't delete me
  EAS_StackingAction* stack_action_       = nullptr;    // don't delete me
  EAS_UserEventAction* event_action_      = nullptr;    // don't delete me
  G4VExceptionHandler* exception_handler_ = nullptr;    // delete me
};

} } } // namespace calin::simulation::generator
