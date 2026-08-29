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
#include <G4IonTable.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4EmExtraPhysics.hh>
#include <G4DecayPhysics.hh>
#include <G4RadioactiveDecayPhysics.hh>
#include <G4HadronElasticPhysics.hh>
#include <G4StoppingPhysics.hh>
#include <G4IonPhysics.hh>
#include <G4NeutronTrackingCut.hh>
#include <G4PhysListFactory.hh>

#include <math/rng.hpp>
#include <simulation/geant4_shower_generator.hpp>
#include <simulation/geant4_shower_generator_internals.hpp>
#include <provenance/chronicle.hpp>

using calin::math::constants::g4_1_c;

using namespace calin::simulation::geant4_shower_generator;
using namespace calin::util::log;

Geant4ShowerGenerator::
Geant4ShowerGenerator(calin::simulation::atmosphere::Atmosphere* atm,
                      config_type config,
                      calin::simulation::world_magnetic_model::FieldVsElevation* bfield,
                      bool adopt_atm, bool adopt_bfield):
  calin::simulation::tracker::ShowerGenerator(),
  atm_(atm), adopt_atm_(adopt_atm), config_(config),
  bfield_(bfield), adopt_bfield_(adopt_bfield)
{
  construct();
}

Geant4ShowerGenerator::
Geant4ShowerGenerator(calin::simulation::atmosphere::Atmosphere* atm,
                      unsigned num_atm_layers, double zground, double ztop_of_atmosphere,
                      calin::simulation::world_magnetic_model::FieldVsElevation* bfield,
                      VerbosityLevel verbose_level, uint32_t seed,
                      double default_cut_value_cm,
                      bool adopt_atm, bool adopt_bfield):
  Geant4ShowerGenerator(atm, 
    customized_config(num_atm_layers, zground, ztop_of_atmosphere, verbose_level, seed, default_cut_value_cm), 
    bfield, adopt_atm, adopt_bfield)
{
  // nothing to see here
}

void Geant4ShowerGenerator::construct()
{
  uint32_t seed = 0;
  while(seed == 0) {
    seed = calin::math::rng::RNG::uint32_from_random_device(); 
  }
  config_.set_seed(seed);
  CLHEP::HepRandom::setTheSeed(seed);
  calin::provenance::chronicle::register_external_rng_open(seed, "CLHEP::HepRandom",
    __PRETTY_FUNCTION__);

  // get the pointer to the User Interface manager
  ui_manager_ = G4UImanager::GetUIpointer();

  // construct a session which receives G4cout/G4cerr
  calin::util::log::Level cout_level = calin::util::log::VERBOSE;
  calin::util::log::Level cerr_level = calin::util::log::WARNING;
  G4int verbose_everything = 0;
  G4int verbose_event = 0;
  G4int verbose_track = 0;

  switch(config_.verbosity())
  {
  case calin::ix::simulation::geant4_shower_generator::SUPPRESSED_ALL:
    cerr_level = calin::util::log::DISCARD;
    // fall through
  case calin::ix::simulation::geant4_shower_generator::SUPRESSED_STDOUT:
    cout_level = calin::util::log::DISCARD;
    break;
  case calin::ix::simulation::geant4_shower_generator::NORMAL:
  default:
    break;
  case calin::ix::simulation::geant4_shower_generator::VERBOSE_EVERYTHING:
    verbose_everything = 1;
    // fall through
  case calin::ix::simulation::geant4_shower_generator::VERBOSE_TRACKING:
    verbose_track = 1;
    // fall through
  case calin::ix::simulation::geant4_shower_generator::VERBOSE_EVENT:
    verbose_event = 1;
    break;
  }

  // construct the default run manager
  run_manager_ = new G4RunManager;

  // ---------------------------------------------------------------------------
  // This should be done before the run manager is constructed but we must wait
  // until the run manager doesn't trash our choice of Exception handler

  // set the UI through the logger
  ui_session_ = new CoutCerrLogger(cout_level, cerr_level);
  ui_manager_->SetCoutDestination(ui_session_);

  // set verbosities
  ui_manager_->
    ApplyCommand(verbose_everything?"/run/verbose 1":"/run/verbose 0");
  ui_manager_->
    ApplyCommand(verbose_event?"/event/verbose 1":"/event/verbose 0");
  ui_manager_->
    ApplyCommand(verbose_track?"/tracking/verbose 1":"/tracking/verbose 0");

  // construct exception handler that avoids abort
  exception_handler_ = new EAS_ExceptionHandler();
  // ---------------------------------------------------------------------------

  G4PhysListFactory physListFactory;

  if(config_.physics_list() != "" and !physListFactory.IsReferencePhysList(config_.physics_list())) {
    LOG(WARNING) << "Physics list \"" << config_.physics_list() << "\" is not recognized. Defaulting to FTFP_BERT.";
    config_.set_physics_list("FTFP_BERT");
  } else if(config_.physics_list() == "") {
    config_.set_physics_list("FTFP_BERT");
  } 

  G4VModularPhysicsList* physlist = physListFactory.GetReferencePhysList(config_.physics_list());
  physlist->SetDefaultCutValue(config_.tracking_cut_scale()*CLHEP::cm);
  physlist->SetVerboseLevel(verbose_everything);
  run_manager_->SetUserInitialization(physlist);

  // Run pre-initialisation commands 
  G4cout.flush();
  G4cerr.flush();
  for(const auto& c : config_.pre_init_commands()) {
    auto retval = apply_command(c);
    if(retval != 0) {
      LOG(WARNING) << "Command: \"" << c << "\" returned " << retval;
    }
  }
  
  // Construct detector
  EAS_FlatDetectorConstruction* detector_constructor =
      new EAS_FlatDetectorConstruction(atm_, config_.num_atm_layers(), 
        config_.zground(), config_.ztop_of_atmosphere(), bfield_,
        config_.detector_box_size(), config_.material());
  run_manager_->SetUserInitialization(detector_constructor);

  event_action_ = new EAS_UserEventAction();
  run_manager_->SetUserAction(event_action_);

  step_action_ = new EAS_SteppingAction(detector_constructor);
  run_manager_->SetUserAction(step_action_);

  gen_action_ = new EAS_PrimaryGeneratorAction();
  run_manager_->SetUserAction(gen_action_);

  // Set minimum energy cut if requested
  if(config_.minimum_energy_cut() > 0) {
    set_minimum_energy_cut(config_.minimum_energy_cut());
  }

  // initialize G4 kernel
  run_manager_->Initialize();

  if(config_.enable_ions()) {
    // Make the ions we support
    // G4IonTable::GetIonTable()->GetIon(2, 4, 0.0); // Helium - already made by G4 as "alpha"
    G4IonTable::GetIonTable()->GetIon(6, 12, 0.0); // Carbon
    G4IonTable::GetIonTable()->GetIon(8, 16, 0.0); // Oxygen
    G4IonTable::GetIonTable()->GetIon(12, 24, 0.0); // Magnesium
    G4IonTable::GetIonTable()->GetIon(14, 28, 0.0); // Silicon
    G4IonTable::GetIonTable()->GetIon(26, 56, 0.0); // Iron
  }
}

Geant4ShowerGenerator::~Geant4ShowerGenerator()
{
  delete run_manager_;
  delete ui_session_;
  delete exception_handler_;
}

int Geant4ShowerGenerator::apply_command(const std::string command)
{
  if(command == "" or command == "help" or command == "?" or command == "/help") {
    LOG(INFO) 
      << "Usage: apply_command(\"/path/to/command [parameter...]\")\n"
      << "For example: apply_command(\"/control/manual /process\")\n"
      << "Note that STDOUT must not be suppressed to see the output of the command.\n";
    return 0;
  }
  auto retval = ui_manager_->ApplyCommand(command);
  G4cout.flush();
  G4cerr.flush();
  return retval;
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
generate_showers(calin::simulation::tracker::TrackVisitor* visitor,
                unsigned num_events,
                calin::simulation::tracker::ParticleType type,
                double total_energy,
                const Eigen::Vector3d& x0,
                const Eigen::Vector3d& u0,
                double ct0, double weight)
{
  event_action_->set_visitor(visitor);
  step_action_->set_visitor(visitor);

  G4ParticleTable* particle_table = G4ParticleTable::GetParticleTable();
  auto pdg_type = particle_type_to_pdg_type(type);
  if(std::abs(pdg_type) > 1000020040 and config_.enable_ions() == false) {
    throw std::invalid_argument(
      "Cannot simulate ion primaries unless \"enable_ions\" option is selected");
  }
  G4ParticleDefinition* particle = particle_table->FindParticle(pdg_type);

  double kinetic_energy = total_energy * CLHEP::MeV - particle->GetPDGMass();
  if(kinetic_energy < 0) {
    throw std::invalid_argument(
      "Total energy must be larger than particle rest mass ("
      + std::to_string(particle->GetPDGMass()/CLHEP::MeV) + " MeV)");
  }

  G4ThreeVector position;
  eigen_to_g4vec(position, x0, CLHEP::cm);

  G4ThreeVector momentum_direction;
  eigen_to_g4vec(momentum_direction, u0);

  auto* particle_generator = gen_action_->particleGenerator();
  particle_generator->SetParticleDefinition(particle);
  particle_generator->SetParticleEnergy(kinetic_energy);
  particle_generator->SetParticlePosition(position);
  particle_generator->SetParticleMomentumDirection(momentum_direction);
  particle_generator->SetParticleTime(ct0 * g4_1_c * CLHEP::ns);
  particle_generator->SetParticleWeight(weight);

  // start a run
  run_manager_->BeamOn(num_events);

  event_action_->set_visitor(nullptr);
  step_action_->set_visitor(nullptr);

  G4cout.flush();
  G4cerr.flush();
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

std::string Geant4ShowerGenerator::pdg_type_to_string(int pdg_type)
{
  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  if(table) {
    const G4ParticleDefinition* particle = table->FindParticle(pdg_type);
    if(particle) {
      return particle->GetParticleName();
    }
  }
  return std::to_string(pdg_type);
}

double Geant4ShowerGenerator::pdg_type_to_mass(int pdg_type)
{
  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  if(table) {
    const G4ParticleDefinition* particle = table->FindParticle(pdg_type);
    if(particle) {
      return particle->GetPDGMass() / CLHEP::MeV;
    }
  }
  return -1;
}

Geant4ShowerGenerator::config_type Geant4ShowerGenerator::default_config()
{
  config_type config;
  config.set_num_atm_layers(1000);
  config.set_zground(0);
  config.set_ztop_of_atmosphere(100E5);
  config.set_tracking_cut_scale(10);
  config.set_minimum_energy_cut(20);
  config.set_physics_list("FTFP_BERT");
  config.set_detector_box_size(1000E5);
  config.set_material("G4_AIR");
  config.set_seed(0);
  config.set_verbosity(calin::ix::simulation::geant4_shower_generator::SUPRESSED_STDOUT);
  return config;
}

Geant4ShowerGenerator::config_type Geant4ShowerGenerator::customized_config(
  unsigned num_atm_layers, double zground, double ztop_of_atmosphere,
  VerbosityLevel verbose_level, uint32_t seed, double default_cut_value_cm,
  double default_set_minimum_energy_cut_mev)
{
  config_type config = default_config();
  config.set_num_atm_layers(num_atm_layers);
  config.set_zground(zground);
  config.set_ztop_of_atmosphere(ztop_of_atmosphere);
  config.set_tracking_cut_scale(default_cut_value_cm);
  config.set_minimum_energy_cut(20);
  config.set_seed(seed);
  switch(verbose_level) {
    case VerbosityLevel::SUPPRESSED_ALL:
      config.set_verbosity(calin::ix::simulation::geant4_shower_generator::SUPPRESSED_ALL);
      break;
    case VerbosityLevel::SUPRESSED_STDOUT:
      config.set_verbosity(calin::ix::simulation::geant4_shower_generator::SUPRESSED_STDOUT);
      break;
    case VerbosityLevel::NORMAL:
    default:
      config.set_verbosity(calin::ix::simulation::geant4_shower_generator::NORMAL);
      break;
    case VerbosityLevel::VERBOSE_EVENT:
      config.set_verbosity(calin::ix::simulation::geant4_shower_generator::VERBOSE_EVENT);
      break;
    case VerbosityLevel::VERBOSE_TRACKING:
      config.set_verbosity(calin::ix::simulation::geant4_shower_generator::VERBOSE_TRACKING);
      break;
    case VerbosityLevel::VERBOSE_EVERYTHING:
      config.set_verbosity(calin::ix::simulation::geant4_shower_generator::VERBOSE_EVERYTHING);
      break;
  }
  return config;
}
