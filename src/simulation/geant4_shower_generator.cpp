/*

   calin/simulation/geant4_shower_generator.cpp -- Stephen Fegan -- 2015-07-02

   Class to generate extensive air showers using Geant-4

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <G4StateManager.hh>

#include <math/rng.hpp>
#include <simulation/geant4_shower_generator.hpp>
#include <simulation/geant4_shower_generator_internals.hpp>
#include <provenance/chronicle.hpp>

using namespace calin::simulation::geant4_shower_generator;
using namespace calin::util::log;

Geant4ShowerGenerator::
Geant4ShowerGenerator(calin::simulation::tracker::TrackVisitor* visitor,
                      calin::simulation::atmosphere::Atmosphere* atm,
                      config_type config,
                      calin::simulation::world_magnetic_model::FieldVsElevation* bfield,
                      bool adopt_visitor, bool adopt_atm, bool adopt_bfield):
  visitor_(visitor), adopt_visitor_(adopt_visitor),
  atm_(atm), adopt_atm_(adopt_atm), ztop_of_atm_(config.ztop()), zground_(config.zground()),
  bfield_(bfield), adopt_bfield_(adopt_bfield), seed_(config.seed())
{
  VerbosityLevel verbose_level;
  switch(config.verbosity()) {
    case calin::ix::simulation::geant4_shower_generator::SUPPRESSED_ALL:
      verbose_level = VerbosityLevel::SUPPRESSED_ALL; break;
    case calin::ix::simulation::geant4_shower_generator::SUPRESSED_STDOUT:
      verbose_level = VerbosityLevel::SUPRESSED_STDOUT; break;
    case calin::ix::simulation::geant4_shower_generator::NORMAL:
    default:
      verbose_level = VerbosityLevel::NORMAL; break;
    case calin::ix::simulation::geant4_shower_generator::VERBOSE_EVENT:
      verbose_level = VerbosityLevel::VERBOSE_EVENT; break;
    case calin::ix::simulation::geant4_shower_generator::VERBOSE_TRACKING:
      verbose_level = VerbosityLevel::VERBOSE_TRACKING; break;
    case calin::ix::simulation::geant4_shower_generator::VERBOSE_EVERYTHING:
      verbose_level = VerbosityLevel::VERBOSE_EVERYTHING; break;
  }
  construct(config.num_atm_layers(), verbose_level, config.tracking_cut_scale(),
    config.detector_box_size(), config.material());
}

Geant4ShowerGenerator::
Geant4ShowerGenerator(calin::simulation::tracker::TrackVisitor* visitor,
                      calin::simulation::atmosphere::Atmosphere* atm,
                      unsigned num_atm_layers, double zground, double ztop,
                      calin::simulation::world_magnetic_model::FieldVsElevation* bfield,
                      VerbosityLevel verbose_level, uint32_t seed,
                      double default_cut_value_cm,
                      bool adopt_visitor, bool adopt_atm, bool adopt_bfield):
    visitor_(visitor), adopt_visitor_(adopt_visitor),
    atm_(atm), adopt_atm_(adopt_atm), ztop_of_atm_(ztop), zground_(zground),
    bfield_(bfield), adopt_bfield_(adopt_bfield), seed_(seed)
{
  auto config = default_config();
  construct(num_atm_layers, verbose_level, default_cut_value_cm,
    config.detector_box_size(), config.material());
}

void Geant4ShowerGenerator::construct(unsigned num_atm_layers,
  VerbosityLevel verbose_level, double default_cut_value_cm,
  double detector_size, const std::string& material_name)
{
  while(seed_ == 0)seed_ = calin::math::rng::RNG::uint32_from_random_device();
  CLHEP::HepRandom::setTheSeed(seed_);
  calin::provenance::chronicle::register_external_rng(seed_, "CLHEP::HepRandom",
    __PRETTY_FUNCTION__);

  // get the pointer to the User Interface manager
  ui_manager_ = G4UImanager::GetUIpointer();

  // construct a session which receives G4cout/G4cerr

  calin::util::log::Level cout_level = calin::util::log::VERBOSE;
  calin::util::log::Level cerr_level = calin::util::log::WARNING;
  G4int verbose_everything = 0;
  G4int verbose_event = 0;
  G4int verbose_track = 0;
  switch(verbose_level)
  {
    case VerbosityLevel::SUPPRESSED_ALL:
      cerr_level = calin::util::log::DISCARD;
      // fall through
    case VerbosityLevel::SUPRESSED_STDOUT:
      cout_level = calin::util::log::DISCARD;
      break;
    case VerbosityLevel::NORMAL:
      break;
    case VerbosityLevel::VERBOSE_EVERYTHING:
      verbose_everything = 1;
      // fall through
    case VerbosityLevel::VERBOSE_TRACKING:
      verbose_track = 1;
      // fall through
    case VerbosityLevel::VERBOSE_EVENT:
      verbose_event = 1;
      break;
  }

  // construct the default run manager
  run_manager_ = new G4RunManager;

  // ---------------------------------------------------------------------------
  // This should be done before the run manager is constrcuted but we must wait
  // until the run manager doesn't trash our choice of Exception handler

  // set the UI through the logger
  ui_session_ = new CoutCerrLogger(cout_level, cerr_level);
  ui_manager_->SetCoutDestination(ui_session_);

  // construct exception handler that avoids abort
  exception_handler_ = new EAS_ExceptionHandler();
  // ---------------------------------------------------------------------------

  // set mandatory initialization classes
  FTFP_BERT* physlist = new FTFP_BERT(verbose_everything);
  physlist->SetDefaultCutValue(default_cut_value_cm*CLHEP::cm);
  physlist->SetVerboseLevel(verbose_everything);
  run_manager_->SetUserInitialization(physlist);

  EAS_FlatDetectorConstruction* detector_constructor =
      new EAS_FlatDetectorConstruction(atm_, num_atm_layers, zground_, ztop_of_atm_, bfield_,
        detector_size, material_name);
  run_manager_->SetUserInitialization(detector_constructor);

  event_action_ = new EAS_UserEventAction(visitor_);
  run_manager_->SetUserAction(event_action_);

  step_action_ = new EAS_SteppingAction(visitor_, detector_constructor);
  run_manager_->SetUserAction(step_action_);

  gen_action_ = new EAS_PrimaryGeneratorAction();
  run_manager_->SetUserAction(gen_action_);

  //run_manager_->SetUserInitialization(new MyUserActionInitialization);

  // initialize G4 kernel
  run_manager_->Initialize();

  // set verbosities
  ui_manager_->
      ApplyCommand(verbose_everything?"/run/verbose 1":"/run/verbose 0");
  ui_manager_->
      ApplyCommand(verbose_event?"/event/verbose 1":"/event/verbose 0");
  ui_manager_->
      ApplyCommand(verbose_track?"/tracking/verbose 1":"/tracking/verbose 0");
}

Geant4ShowerGenerator::~Geant4ShowerGenerator()
{
  delete run_manager_;
  delete ui_session_;
  delete exception_handler_;
}

void Geant4ShowerGenerator::set_minimum_energy_cut(double emin_mev)
{
  step_action_->setEminCut(emin_mev);
  if(stack_action_ == nullptr)
  {
    stack_action_ = new EAS_StackingAction();
    run_manager_->SetUserAction(stack_action_);
  }
  stack_action_->setEminCut(emin_mev);
}

void Geant4ShowerGenerator::
generate_showers(unsigned num_events,
                calin::simulation::tracker::ParticleType type,
                double total_energy,
                const Eigen::Vector3d& x0,
                const Eigen::Vector3d& u0,
                double weight)
{
  G4GeneralParticleSource* gps = new G4GeneralParticleSource();

  // default particle kinematic
  G4ParticleTable* particle_table = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle
      = particle_table->FindParticle(particle_type_to_pdg_type(type));

  if(total_energy * CLHEP::MeV < particle->GetPDGMass())
    throw std::invalid_argument(
      "Total energy must be larger than particle rest mass ("
      + std::to_string(particle->GetPDGMass()/CLHEP::MeV) + " MeV)");

  G4ThreeVector position;
  eigen_to_g4vec(position, x0, CLHEP::cm);

  G4ThreeVector momentum_direction;
  eigen_to_g4vec(momentum_direction, u0);

  G4SingleParticleSource* sps = gps->GetCurrentSource();
  sps->SetParticleDefinition(particle);
  sps->GetPosDist()->SetPosDisType("Point");
  sps->GetPosDist()->SetCentreCoords(position);
  sps->GetAngDist()->SetAngDistType("planar");
  sps->GetAngDist()->SetParticleMomentumDirection(momentum_direction);
  sps->GetEneDist()->SetEnergyDisType("Mono");
  sps->GetEneDist()->SetMonoEnergy(total_energy * CLHEP::MeV - particle->GetPDGMass());

  gen_action_->setGPS(gps);

  // start a run
  run_manager_->BeamOn(num_events);
}

double Geant4ShowerGenerator::
get_material_density(const std::string& material_name)
{
  G4NistManager* nist = G4NistManager::Instance();
  if(nist == nullptr)throw std::runtime_error("Could not construct NIST manager");
  G4Material* material = nist->FindOrBuildMaterial(material_name);
  if(material == nullptr)throw std::runtime_error("Could not find material " + material_name);
  return material->GetDensity() / (CLHEP::g/CLHEP::cm3);
}

Geant4ShowerGenerator::config_type Geant4ShowerGenerator::default_config()
{
  config_type config;
  config.set_num_atm_layers(1000);
  config.set_zground(0);
  config.set_ztop(100E5);
  config.set_tracking_cut_scale(10);
  config.set_detector_box_size(1000E5);
  config.set_material("G4_AIR");
  config.set_seed(0);
  config.set_verbosity(calin::ix::simulation::geant4_shower_generator::SUPRESSED_STDOUT);
  return config;
}
