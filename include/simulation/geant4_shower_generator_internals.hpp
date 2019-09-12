/*

   calin/simulation/geant4_shower_generator_internals.hpp
   Stephen Fegan -- 2015-07-03

   Internals of Geant-4 extensive air shower generator

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

#pragma once

#include<vector>
#include<cassert>
#include<limits>
#include<map>

#include"simulation/atmosphere.hpp"
#include"simulation/tracker.hpp"
#include"math/special.hpp"
#include"util/log.hpp"
#include"simulation/world_magnetic_model.hpp"

#include"calin_global_definitions.hpp"

using calin::math::special::SQR;

// Units:
// Height, z:    cm
// Energy, e:    MeV

#include <G4RunManager.hh>
#include <G4UImanager.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4UserStackingAction.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4LogicalVolume.hh>
#include <G4Box.hh>
#include <G4Material.hh>
#include <G4PVPlacement.hh>
#include <G4NistManager.hh>
#include <G4VModularPhysicsList.hh>
#include <G4PhysListFactory.hh>
#include <G4GeneralParticleSource.hh>
#include <G4VUserActionInitialization.hh>
#include <G4UserEventAction.hh>
#include <FTFP_BERT.hh>
#include <G4UIsession.hh>
#include <G4TrackStatus.hh>
#include <G4ExceptionHandler.hh>

namespace calin { namespace simulation { namespace geant4_shower_generator {

void g4vec_to_eigen(Eigen::Vector3d& evec, const G4ThreeVector& g4vec);
void g4vec_to_eigen(Eigen::Vector3d& evec, const G4ThreeVector& g4vec,
                    double to_units);
void eigen_to_g4vec(G4ThreeVector& g4vec, const Eigen::Vector3d& evec);
void eigen_to_g4vec(G4ThreeVector& g4vec, const Eigen::Vector3d& evec,
                    double from_units);

inline bool apply_kinetic_energy_cut(int pdg_type)
{
  return pdg_type==11 // electron
      or pdg_type==2212 // proton
      or pdg_type==2112 // neutron
      ;
}

class EAS_StackingAction: public G4UserStackingAction
{
public:
  EAS_StackingAction();
  virtual ~EAS_StackingAction();

  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*) override;

  void setEminCut(double emin_MeV) { ecut_ = emin_MeV/CLHEP::MeV; }

protected:
  double ecut_ = 0;
};

class EAS_DetectorConstruction;

class EAS_SteppingAction: public G4UserSteppingAction
{
public:
  EAS_SteppingAction(calin::simulation::tracker::TrackVisitor* visitor,
    EAS_DetectorConstruction* detector_geometry = nullptr);
  virtual ~EAS_SteppingAction();

  void UserSteppingAction(const G4Step*) override;

  void setEminCut(double emin_MeV) { ecut_ = emin_MeV/CLHEP::MeV; }

protected:
  double ecut_  = 0;
  calin::simulation::tracker::TrackVisitor* visitor_ = nullptr;
  EAS_DetectorConstruction* detector_geometry_ = nullptr;
};

class EAS_PrimaryGeneratorAction: public G4VUserPrimaryGeneratorAction
{
public:
  EAS_PrimaryGeneratorAction();
  virtual ~EAS_PrimaryGeneratorAction();

  void GeneratePrimaries(G4Event* the_event) override;

  void setGPS(G4GeneralParticleSource* particle_source) {
    delete particle_source_; particle_source_ = particle_source; }
protected:
  G4GeneralParticleSource* particle_source_ = nullptr;
};

class EAS_DetectorConstruction: public G4VUserDetectorConstruction
{
public:
  EAS_DetectorConstruction();
  virtual ~EAS_DetectorConstruction();
  G4VPhysicalVolume* Construct() override = 0;
  void ConstructSDandField() override = 0;
  virtual bool ray_intersects_detector(const Eigen::Vector3d& pos,
    const Eigen::Vector3d& dir) = 0;
};

class EAS_FlatDetectorConstruction: public EAS_DetectorConstruction
{
public:
  EAS_FlatDetectorConstruction(calin::simulation::atmosphere::Atmosphere* atm,
      unsigned num_atm_layers, double zground_cm, double ztop_of_atm_cm,
      double layer_side_cm = 100E5,
      calin::simulation::world_magnetic_model::FieldVsElevation* bfield=nullptr);
  EAS_FlatDetectorConstruction(calin::simulation::atmosphere::Atmosphere* atm,
      unsigned num_atm_layers, double zground_cm, double ztop_of_atm_cm,
      calin::simulation::world_magnetic_model::FieldVsElevation* bfield,
      double layer_side_cm = 100E5):
    EAS_FlatDetectorConstruction(atm, num_atm_layers, zground_cm, ztop_of_atm_cm,
      layer_side_cm, bfield) { /* nothing to see here */ }

  virtual ~EAS_FlatDetectorConstruction();
  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;
  bool ray_intersects_detector(const Eigen::Vector3d& pos,
    const Eigen::Vector3d& dir) override;
private:
  calin::simulation::atmosphere::Atmosphere* atm_;
  unsigned num_atm_layers_;
  double zground_cm_;
  double ztop_of_atm_cm_;
  double layer_side_cm_;
  Eigen::Vector3d min_corner_;
  Eigen::Vector3d max_corner_;
  calin::simulation::world_magnetic_model::FieldVsElevation* bfield_;
  std::map<G4LogicalVolume*, G4ThreeVector> logical_bfield_;
};

class EAS_UserEventAction: public G4UserEventAction
{
public:
  EAS_UserEventAction(calin::simulation::tracker::TrackVisitor* visitor);
  virtual ~EAS_UserEventAction();
  void BeginOfEventAction(const G4Event *anEvent) override;
  void EndOfEventAction(const G4Event *anEvent) override;
private:
  calin::simulation::tracker::TrackVisitor* visitor_;
};

class EAS_ExceptionHandler: public G4ExceptionHandler
{
public:
  EAS_ExceptionHandler();
  virtual ~EAS_ExceptionHandler();
  G4bool Notify(const char* originOfException, const char *exceptionCode,
    G4ExceptionSeverity severity, const char * description) override;
};

class CoutCerrLogger: public G4UIsession
{
public:
  CoutCerrLogger(calin::util::log::Level cout_level = calin::util::log::VERBOSE,
                 calin::util::log::Level cerr_level = calin::util::log::WARNING,
                 calin::util::log::Logger* cout_logger =
                 calin::util::log::default_logger(),
                 calin::util::log::Logger* cerr_logger =
                 calin::util::log::default_logger(),
                 bool adopt_cout_logger = false,
                 bool adopt_cerr_logger = false);
  virtual ~CoutCerrLogger();
  G4int ReceiveG4cout(const G4String& cout_string) override;
  G4int ReceiveG4cerr(const G4String& cerr_string) override;
protected:
  calin::util::log::Level cout_level_;
  calin::util::log::Level cerr_level_;
  calin::util::log::Logger* cout_logger_;
  calin::util::log::Logger* cerr_logger_;
  bool adopt_cout_logger_;
  bool adopt_cerr_logger_;
};


} } } // namespace calin::simulation::generator
