/*

   calin/simulation/geant4_shower_generator_internals.cpp
   Stephen Fegan -- 2015-07-07

   Internals of Geant-4 extensive air shower generator

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

#include<stdexcept>
#include<G4GeometryManager.hh>
#include<G4UniformMagField.hh>
#include<G4AutoDelete.hh>
#include<G4FieldManager.hh>

#include <simulation/geant4_shower_generator_internals.hpp>
#include <io/log.hpp>
#include <math/geometry.hpp>

using namespace calin::simulation::geant4_shower_generator;
using namespace calin::io::log;

void calin::simulation::geant4_shower_generator::
g4vec_to_eigen(Eigen::Vector3d& evec, const G4ThreeVector& g4vec)
{
  evec[0] = g4vec[0];
  evec[1] = g4vec[1];
  evec[2] = g4vec[2];
}

void calin::simulation::geant4_shower_generator::
g4vec_to_eigen(Eigen::Vector3d& evec, const G4ThreeVector& g4vec,
                    double to_units)
{
  evec[0] = g4vec[0]/to_units;
  evec[1] = g4vec[1]/to_units;
  evec[2] = g4vec[2]/to_units;
}

void calin::simulation::geant4_shower_generator::
eigen_to_g4vec(G4ThreeVector& g4vec, const Eigen::Vector3d& evec)
{
  g4vec[0] = evec[0];
  g4vec[1] = evec[1];
  g4vec[2] = evec[2];
}

void calin::simulation::geant4_shower_generator::
eigen_to_g4vec(G4ThreeVector& g4vec, const Eigen::Vector3d& evec,
               double from_units)
{
  g4vec[0] = evec[0]*from_units;
  g4vec[1] = evec[1]*from_units;
  g4vec[2] = evec[2]*from_units;
}

// ============================================================================
//
// EAS_SteppingAction - Stacking action - kill low energy particles
//
// ============================================================================

EAS_StackingAction::EAS_StackingAction(): G4UserStackingAction()
{
  // nothing to see here
}

EAS_StackingAction::~EAS_StackingAction()
{
  // nothing to see here
}

G4ClassificationOfNewTrack EAS_StackingAction::
ClassifyNewTrack(const G4Track* track)
{
  const G4ParticleDefinition* pdg_info = track->GetParticleDefinition();
  if(apply_kinetic_energy_cut(pdg_info->GetPDGEncoding()))
    return track->GetKineticEnergy() < ecut_ ? fKill : fUrgent;
  else
    return track->GetTotalEnergy() < ecut_ ? fKill : fUrgent;
}

// ============================================================================
//
// EAS_SteppingAction - Stepping action
//
// ============================================================================

EAS_SteppingAction::
EAS_SteppingAction(calin::simulation::tracker::TrackVisitor* visitor,
    EAS_DetectorConstruction* detector_geometry):
  G4UserSteppingAction(), visitor_(visitor),
  detector_geometry_(detector_geometry)
{
  /* nothing to see here */
}

EAS_SteppingAction::~EAS_SteppingAction()
{
  // nothing to see here
}

void EAS_SteppingAction::UserSteppingAction(const G4Step* the_step)
{
  const G4StepPoint* pre_step_pt = the_step->GetPreStepPoint();

  double pre_step_pt_etot = pre_step_pt->GetTotalEnergy();

  const G4ParticleDefinition* pdg_info =
      the_step->GetTrack()->GetParticleDefinition();

  if(apply_kinetic_energy_cut(pdg_info->GetPDGEncoding()))
  {
    if(pre_step_pt->GetKineticEnergy() < ecut_) {
      the_step->GetTrack()->SetTrackStatus(fStopAndKill);
      return;
    }
  }
  else
  {
    if(pre_step_pt_etot < ecut_)
    {
      the_step->GetTrack()->SetTrackStatus(fStopAndKill);
      return;
    }
  }

  calin::simulation::tracker::Track track;

  track.pdg_type = pdg_info->GetPDGEncoding();
  track.q        = pdg_info->GetPDGCharge();
  track.mass     = pdg_info->GetPDGMass()/CLHEP::MeV;
  track.type     = tracker::pdg_type_to_particle_type(track.pdg_type);

  const G4ThreeVector& pre_step_pt_posn = pre_step_pt->GetPosition();
  track.e0       = pre_step_pt_etot/CLHEP::MeV;
  g4vec_to_eigen(track.x0, pre_step_pt_posn, CLHEP::cm);
  g4vec_to_eigen(track.u0, pre_step_pt->GetMomentumDirection());
  track.t0       = pre_step_pt->GetGlobalTime()/CLHEP::ns;

  const G4StepPoint* post_step_pt = the_step->GetPostStepPoint();
  const G4ThreeVector& post_step_pt_posn = post_step_pt->GetPosition();
  track.e1       = post_step_pt->GetTotalEnergy()/CLHEP::MeV;
  g4vec_to_eigen(track.x1, post_step_pt_posn, CLHEP::cm);
  g4vec_to_eigen(track.u1, post_step_pt->GetMomentumDirection());
  track.t1       = post_step_pt->GetGlobalTime()/CLHEP::ns;

  track.dx_hat   = track.x1 - track.x0;
  track.dx       = track.dx_hat.norm();
  track.dx_hat/=track.dx;
  track.de       = track.e1 - track.e0;
  track.dt       = track.t1 - track.t0;
  track.weight   = pre_step_pt->GetWeight();

#if 0
  static int nprint=0;
  if(nprint<100 and pdg_info->GetPDGEncoding()==2112) {
    LOG(INFO) << pre_step_pt->GetKineticEnergy() << ' ' << track.dx;
    ++nprint;
  }
#endif

  bool kill_track = false;
  visitor_->visit_track(track, kill_track);
  if(kill_track)the_step->GetTrack()->SetTrackStatus(fStopAndKill);

  if(detector_geometry_ and
    not detector_geometry_->ray_intersects_detector(track.x1, track.u1))
  {
    the_step->GetTrack()->SetTrackStatus(fStopAndKill);
    return;
  }
}

// ============================================================================
//
// EAS_PrimaryGeneratorAction - create primary particles on demand using
//                              supplied particle source
//
// ============================================================================

EAS_PrimaryGeneratorAction::EAS_PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction()
{
  // nothing to see here
}

EAS_PrimaryGeneratorAction::~EAS_PrimaryGeneratorAction()
{
  delete(particle_source_);
}

void EAS_PrimaryGeneratorAction::GeneratePrimaries(G4Event* the_event)
{
  // this function is called at the begining of event
  particle_source_->GeneratePrimaryVertex(the_event);
}

// ============================================================================
//
// Detector construction
//
// ============================================================================

EAS_DetectorConstruction::EAS_DetectorConstruction()
    : G4VUserDetectorConstruction()
{
  // nothing to see here
}

EAS_DetectorConstruction::~EAS_DetectorConstruction()
{
  // nothing to see here
}

EAS_FlatDetectorConstruction::
EAS_FlatDetectorConstruction(calin::simulation::atmosphere::Atmosphere* atm,
  unsigned num_atm_layers, double zground_cm, double ztop_of_atm_cm,
  double layer_side_cm,
  calin::simulation::world_magnetic_model::FieldVsElevation* bfield)
    : EAS_DetectorConstruction(), atm_(atm), num_atm_layers_(num_atm_layers),
      zground_cm_(zground_cm), ztop_of_atm_cm_(ztop_of_atm_cm),
      layer_side_cm_(layer_side_cm),
      min_corner_(-layer_side_cm, -layer_side_cm, zground_cm),
      max_corner_( layer_side_cm,  layer_side_cm, ztop_of_atm_cm),
      bfield_(bfield)
{
  // nothing to see here
}

EAS_FlatDetectorConstruction::~EAS_FlatDetectorConstruction()
{
  // nothing to see here
}

G4VPhysicalVolume* EAS_FlatDetectorConstruction::Construct()
{
  constexpr double eps = std::numeric_limits<double>::epsilon();

  std::vector<atmosphere::AtmSlice> slices =
      atm_->make_atm_slices(num_atm_layers_, ztop_of_atm_cm_, zground_cm_);

  G4NistManager* nist = G4NistManager::Instance();

  G4double world_hx = layer_side_cm_*CLHEP::cm;
  G4double world_hy = layer_side_cm_*CLHEP::cm;
  G4double world_hz = ztop_of_atm_cm_*CLHEP::cm;

  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(
    std::max({world_hx, world_hy, world_hz}));

  G4Box* world_box
      = new G4Box("BOX_World",
                  world_hx*(1+eps), world_hy*(1+eps), world_hz*(1+eps));

  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  G4LogicalVolume* world_logical
      = new G4LogicalVolume(world_box, world_mat, std::string("VOL_World"));

  G4VPhysicalVolume* world_physical
      = new G4PVPlacement(0,                           // no rotation
                          G4ThreeVector(0, 0, 0), // translation
                          world_logical    ,           // its logical volume
                          std::string("PHY_WORLD"),    // its name
                          0,                           // its mother volume
                          false,                       // boolean operations
                          0);                          // its copy number}

  for(const auto& islice : slices)
  {
    G4double density = islice.rho * CLHEP::g/CLHEP::cm3;
    std::string name {"LAYER_"};
    name += std::to_string(int(std::floor(islice.zb)));
    name += "_";
    name += std::to_string(int(std::floor(islice.zt)));

    G4Material* air =
        nist->BuildMaterialWithNewDensity(std::string("AIR_")+name,
                                          "G4_AIR",density);

    G4double layer_hz = 0.5*(islice.zt-islice.zb)*CLHEP::cm;
    G4double pos_z = 0.5*(islice.zt+islice.zb)*CLHEP::cm;
#if 0
    if(islice.zb-zground_cm_ < eps)
    {
      // Add 1mm guard to lowest layer simplify cut on z
      layer_hz += 0.5*CLHEP::mm;
      pos_z -= 0.5*CLHEP::mm;
    }
    layer_hz *= (1-eps);
#else
    //layer_hz -= 0.5*CLHEP::mm;
#endif
    G4Box* layer_box
        = new G4Box(std::string("BOX_")+name, world_hx, world_hy, layer_hz);

    G4LogicalVolume* layer_logical
        = new G4LogicalVolume(layer_box, air, std::string("VOL_")+name);

    if(bfield_)
      eigen_to_g4vec(logical_bfield_[layer_logical],
        bfield_->field_nT(pos_z), 1e-9*CLHEP::tesla);

    G4VPhysicalVolume* layer_physical __attribute__((__unused__))
        = new G4PVPlacement(0,                          // no rotation
                            G4ThreeVector(0, 0, pos_z), // translation
                            layer_logical    ,          // its logical volume
                            std::string("PHY_")+name,   // its name
                            world_logical,              // its mother volume
                            false,                      // boolean operations
                            0);                         // its copy number}
  }

  return world_physical;
}

void EAS_FlatDetectorConstruction::ConstructSDandField()
{
  for(auto& ilogical_bfield : logical_bfield_)
  {
    auto* mag_field = new G4UniformMagField(ilogical_bfield.second);
    auto* field_mgr = new G4FieldManager();
    field_mgr->SetDetectorField(mag_field);
    field_mgr->CreateChordFinder(mag_field);
    G4bool force_to_all_daughters = true;
    ilogical_bfield.first->SetFieldManager(field_mgr, force_to_all_daughters);
    // Register the field and its manager for deleting
    G4AutoDelete::Register(mag_field);
    G4AutoDelete::Register(field_mgr);
  }
}

bool EAS_FlatDetectorConstruction::
ray_intersects_detector(const Eigen::Vector3d& pos, const Eigen::Vector3d& dir)
{
  return calin::math::geometry::
    box_has_future_intersection(min_corner_,max_corner_,pos,dir);
}

// ============================================================================
//
// EAS_UserEventAction - send visit_event and leave_event messages
//
// ============================================================================


EAS_UserEventAction::
EAS_UserEventAction(calin::simulation::tracker::TrackVisitor* visitor):
  G4UserEventAction(), visitor_(visitor)
{
  // nothing to see here
}

EAS_UserEventAction::~EAS_UserEventAction()
{
  // nothing to see here
}

void EAS_UserEventAction::BeginOfEventAction(const G4Event *anEvent)
{
  calin::simulation::tracker::Event event;
  event.event_id = anEvent->GetEventID();

  assert(anEvent->GetNumberOfPrimaryVertex() == 1);
  const auto* vertex = anEvent->GetPrimaryVertex(0);
  assert(vertex->GetNumberOfParticle() == 1);
  const auto* primary = vertex->GetPrimary(0);

  const auto* pdg_info = primary->GetG4code();
  event.pdg_type = pdg_info->GetPDGEncoding();
  event.q        = pdg_info->GetPDGCharge();
  event.mass     = pdg_info->GetPDGMass()/CLHEP::MeV;
  event.type     = tracker::pdg_type_to_particle_type(event.pdg_type);

  const auto posn = vertex->GetPosition();
  event.e0       = primary->GetTotalEnergy()/CLHEP::MeV;
  g4vec_to_eigen(event.x0, posn, CLHEP::cm);
  g4vec_to_eigen(event.u0, primary->GetMomentumDirection());
  event.t0       = 0; // pre_step_pt->GetGlobalTime()/CLHEP::ns;

  event.weight   = vertex->GetWeight();

  bool kill_event = false;
  visitor_->visit_event(event, kill_event);
  if(kill_event)fpEventManager->AbortCurrentEvent();
}

void EAS_UserEventAction::EndOfEventAction(const G4Event *anEvent)
{
  visitor_->leave_event();
}

// ============================================================================
//
// EAS_ExceptionHandler - stop Geant4 from aborting, rather send C++ exception
//
// ============================================================================

EAS_ExceptionHandler::EAS_ExceptionHandler(): G4ExceptionHandler()
{
  // nothing to see here
}

EAS_ExceptionHandler::~EAS_ExceptionHandler()
{
  // nothing to see here
}

G4bool EAS_ExceptionHandler::Notify(const char* originOfException,
  const char *exceptionCode, G4ExceptionSeverity severity,
  const char * description)
{
  G4bool abort = G4ExceptionHandler::Notify(originOfException, exceptionCode,
    severity, description);
  if(abort)throw std::runtime_error("Geant4 exception");
  return abort;
}

// ============================================================================
//
// CoutCerrLogger - send Geant4 messages to the log system
//
// ============================================================================

CoutCerrLogger::
CoutCerrLogger(calin::io::log::Level cout_level,
               calin::io::log::Level cerr_level,
               calin::io::log::Logger* cout_logger,
               calin::io::log::Logger* cerr_logger,
               bool adopt_cout_logger,
               bool adopt_cerr_logger):
    G4UIsession(), cout_level_(cout_level), cerr_level_(cerr_level),
    cout_logger_(cout_logger), cerr_logger_(cerr_logger),
    adopt_cout_logger_(adopt_cout_logger),
    adopt_cerr_logger_(adopt_cerr_logger)
{
  // nothing to see here
}

CoutCerrLogger::~CoutCerrLogger()
{
  // nothing to see here
}

G4int CoutCerrLogger::ReceiveG4cout(const G4String& cout_string)
{
  cout_logger_->log_message(cout_level_, cout_string);
  return 0;
}

G4int CoutCerrLogger::ReceiveG4cerr(const G4String& cerr_string)
{
  cerr_logger_->log_message(cerr_level_, cerr_string);
  return 0;
}
