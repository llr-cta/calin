/* 

   calin/air_shower/geant4_shower_generator.hpp -- Stephen Fegan -- 2015-06-23

   Class to genereate extensive air showers using Geant-4

*/

#pragma once

#include<vector>
#include<cassert>

#include"air_shower/atmosphere.hpp"
#include"air_shower/tracker.hpp"

#include"package_wide_definitions.hpp"

// Units:
// Height, z:    cm
// Energy, e:    MeV

class G4UImanager;
class G4UIsession;
class G4RunManager;

namespace calin { namespace air_shower { namespace shower_generator {

class EAS_Actions;

enum class VerbosityLevel {
  SUPPRESSED_ALL, SUPRESSED_STDOUT, NORMAL, VERBOSE_EVENT, VERBOSE_TRACKING,
    VERBOSE_EVERYTHING };

class Geant4ShowerGenerator
{
 public:
  Geant4ShowerGenerator(calin::air_shower::tracker::TrackVisitor* visitor,
                        calin::air_shower::atmosphere::Atmosphere* atm,
                        unsigned num_atm_layers, double zground, double ztop,
                        VerbosityLevel verbose_level =
                        VerbosityLevel::SUPRESSED_STDOUT,
                        bool adopt_visitor = false,
                        bool adopt_atm = false);
  virtual ~Geant4ShowerGenerator();

  void setMinimumEnergyCut(double emin_mev);
  void generateShower(calin::air_shower::tracker::ParticleType type,
                      double total_energy,
                      const Eigen::Vector3d& x0 = Eigen::Vector3d(0,0,0),
                      const Eigen::Vector3d& u0 = Eigen::Vector3d(0,0,-1),
                      double weight=1.0);

 protected:
  calin::air_shower::tracker::TrackVisitor* visitor_ = nullptr;
  bool adopt_visitor_ = false;
  calin::air_shower::atmosphere::Atmosphere* atm_ = nullptr;
  bool adopt_atm_ = false;
  double ztop_of_atm_ = 0;
  double zground_ = 0;
  
  G4UImanager* ui_manager_   = nullptr;    // singleton - do not delete!
  G4UIsession* ui_session_   = nullptr;    // delete me
  G4RunManager* run_manager_ = nullptr;    // delete me
  EAS_Actions* action_       = nullptr;    // don't delete me
}; 

} } } // namespace calin::air_shower::generator
