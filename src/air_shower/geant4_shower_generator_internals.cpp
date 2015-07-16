/* 

   calin/air_shower/geant4_shower_generator_internals.cpp 
   Stephen Fegan -- 2015-07-07

   Internals of Geant-4 extensive air shower generator

*/

#include"air_shower/geant4_shower_generator_internals.hpp"

using namespace calin::air_shower::shower_generator;

namespace {

void g4vec_to_eigen(Eigen::Vector3d& evec, const G4ThreeVector& g4vec)
{
  evec[0] = g4vec[0];
  evec[1] = g4vec[1];
  evec[2] = g4vec[2];
}

void g4vec_to_eigen(Eigen::Vector3d& evec, const G4ThreeVector& g4vec,
                    double to_units)
{
  evec[0] = g4vec[0]/to_units;
  evec[1] = g4vec[1]/to_units;
  evec[2] = g4vec[2]/to_units;
}

} // anonymous namespace

// ============================================================================
//
// EAS_Actions - Run / Track / Step actions
//
// ============================================================================

EAS_Actions::EAS_Actions(calin::air_shower::tracker::TrackVisitor* visitor)
    : G4UserSteppingAction(), visitor_(visitor)
{
  /* nothing to see here */
}

EAS_Actions::~EAS_Actions()
{
  // nothing to see here
}

void EAS_Actions::UserSteppingAction(const G4Step* the_step)
{
  const G4StepPoint* pre_step_pt = the_step->GetPreStepPoint();

  double pre_step_pt_etot = pre_step_pt->GetTotalEnergy();
  const G4ThreeVector& pre_step_pt_posn = pre_step_pt->GetPosition();
  if(pre_step_pt_etot < ecut_)
  {
    the_step->GetTrack()->SetTrackStatus(fStopAndKill);
    return;
  }

  if(r2cut_ > 0)
  {
    if(pre_step_pt_posn.x()*pre_step_pt_posn.x()
       + pre_step_pt_posn.y()*pre_step_pt_posn.y()
       + pre_step_pt_posn.z()*(pre_step_pt_posn.z() + 2.0*rzero_)
       + rzero_*rzero_ < r2cut_)
    {
      the_step->GetTrack()->SetTrackStatus(fStopAndKill);
      return;     
    }
  }
  else if(pre_step_pt_posn.z() < zcut_)
  {
    the_step->GetTrack()->SetTrackStatus(fStopAndKill);
    return;     
  }
  
  calin::air_shower::tracker::Track track;

  const G4ParticleDefinition* pdg_info =
      the_step->GetTrack()->GetParticleDefinition();
  track.pdg_type = pdg_info->GetPDGEncoding();
  track.q        = pdg_info->GetPDGCharge();
  track.mass     = pdg_info->GetPDGMass()/CLHEP::MeV;

  using calin::air_shower::tracker::ParticleType;
  switch(track.pdg_type)
  {
    case 22:    track.type = ParticleType::GAMMA; break;
    case 1:     track.type = ParticleType::ELECTRON; break;
    case -1:    track.type = ParticleType::POSITRON; break;
    case 13:    track.type = ParticleType::MUON; break;
    case -13:   track.type = ParticleType::ANTI_MUON; break;
    case 2212:  track.type = ParticleType::PROTON; break;
    case -2212: track.type = ParticleType::ANTI_PROTON; break;
    default:    track.type = ParticleType::OTHER; break;
  };

  track.e0       = pre_step_pt_etot/CLHEP::MeV;
  g4vec_to_eigen(track.x0, pre_step_pt_posn, CLHEP::cm);
  g4vec_to_eigen(track.u0, pre_step_pt->GetMomentumDirection());
  track.t0       = pre_step_pt->GetGlobalTime()/CLHEP::ns;

  const G4StepPoint* post_step_pt = the_step->GetPostStepPoint();
  track.e1       = post_step_pt->GetTotalEnergy()/CLHEP::MeV;
  g4vec_to_eigen(track.x1, post_step_pt->GetPosition(), CLHEP::cm);
  g4vec_to_eigen(track.u1, post_step_pt->GetMomentumDirection());  
  track.t1       = post_step_pt->GetGlobalTime()/CLHEP::ns;

  track.dx       = the_step->GetStepLength()/CLHEP::cm;
  track.weight   = pre_step_pt->GetWeight();
      
  bool kill_track = false;
  visitor_->visitTrack(track, kill_track);
  if(kill_track)the_step->GetTrack()->SetTrackStatus(fStopAndKill);
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
EAS_FlatDetectorConstruction(calin::air_shower::atmosphere::Atmosphere* atm,
                             unsigned num_atm_layers,
                             double zground_cm, double ztop_of_atm_cm,
                             double layer_side_cm)
    : EAS_DetectorConstruction(), atm_(atm), num_atm_layers_(num_atm_layers),
      zground_cm_(zground_cm), ztop_of_atm_cm_(ztop_of_atm_cm),
      layer_side_cm_(layer_side_cm)
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
  G4double world_hz = 0.5*(ztop_of_atm_cm_-zground_cm_)*CLHEP::cm;
  G4Box* world_box
      = new G4Box("BOX_World",
                  world_hx*(1+eps), world_hy*(1+eps), world_hz*(1+eps));

  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  G4LogicalVolume* world_logical
      = new G4LogicalVolume(world_box, world_mat, std::string("VOL_World"));

  G4double wpos_z = 0.5*(ztop_of_atm_cm_-zground_cm_)*CLHEP::cm;
                      
  G4VPhysicalVolume* world_physical
      = new G4PVPlacement(0,                           // no rotation
                          G4ThreeVector(0, 0, wpos_z), // translation
                          world_logical    ,           // its logical volume
                          std::string("PHY_WORLD"),    // its name
                          0,                           // its mother volume
                          false,                       // boolean operations
                          0);                          // its copy number}

  for(const auto& islice : slices)
  {
    G4double density = islice.rho * CLHEP::g/CLHEP::cm3;
    std::string name {"LAYER_"};
    name += std::to_string(islice.zb);
    name += "_";
    name += std::to_string(islice.zt);
    
    G4Material* air =
        nist->BuildMaterialWithNewDensity(std::string("AIR_")+name,
                                          "G4_AIR",density);

    G4double layer_hz = 0.5*(islice.zt-islice.zb)*CLHEP::cm*(1-eps);

    G4Box* layer_box
        = new G4Box(std::string("BOX_")+name, world_hx, world_hy, layer_hz);
    
    G4LogicalVolume* layer_logical
        = new G4LogicalVolume(layer_box, air, std::string("VOL_")+name);

    G4double pos_z = 0.5*(islice.zt+islice.zb)*CLHEP::cm;
  
    G4VPhysicalVolume* layer_physical
        = new G4PVPlacement(0,                          // no rotation
                            G4ThreeVector(0, 0, pos_z), // translation
                            world_logical    ,          // its logical volume
                            std::string("PHY_")+name,   // its name
                            world_logical,              // its mother volume
                            false,                      // boolean operations
                            0);                         // its copy number}
  }
  
  return world_physical;
}

void EAS_FlatDetectorConstruction::ConstructSDandField()
{
  
};


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
