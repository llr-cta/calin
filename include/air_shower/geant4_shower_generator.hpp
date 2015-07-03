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
// Energy, e:    eV

class G4UImanager;
class G4UIsession;
class G4RunManager;

namespace calin { namespace air_shower { namespace shower_generator {

class ThresholdEnergyKiller;

class Geant4ShowerGenerator
{
 public:
  Geant4ShowerGenerator(calin::air_shower::tracker::TrackVisitor* visitor,
                        calin::air_shower::atmosphere::Atmosphere& atm,
                        unsigned num_atm_layers, double zground, double ztop,
                        bool adopt_visitor = false);
  virtual ~Geant4ShowerGenerator();

  void setMinimumEnergyCut(double emin);
  void generateShower(calin::air_shower::tracker::ParticleType type,
                      double total_energy,
                      const Eigen::Vector3d& x0 = Eigen::Vector3d(0,0,0),
                      const Eigen::Vector3d& u0 = Eigen::Vector3d(0,0,-1),
                      double weight=1.0);

 protected:
  calin::air_shower::tracker::TrackVisitor* visitor_ = nullptr;
  bool adopt_visitor_ = false;

  G4UImanager* ui_manager_;   // singleton - do not delete!
  G4UIsession* ui_session_;   // delete me
  G4RunManager* run_manager_; // delete me
  ThresholdEnergyKiller* killer_; // don't delete me
}; 

} } } // namespace calin::air_shower::generator
