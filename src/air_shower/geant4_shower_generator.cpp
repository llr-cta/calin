/* 

   calin/air_shower/geant4_shower_generator.cpp -- Stephen Fegan -- 2015-07-02

   Class to genereate extensive air showers using Geant-4

*/

#include "air_shower/geant4_shower_generator.hpp"

#include <G4RunManager.hh>
#include <G4UImanager.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4LogicalVolume.hh>
#include <G4Box.hh>
#include <G4Material.hh>
#include <G4PVPlacement.hh>
#include <G4NistManager.hh>
#include <G4VModularPhysicsList.hh>
#include <G4PhysListFactory.hh>
#include <G4ParticleGun.hh>
#include <G4VUserActionInitialization.hh>
#include <FTFP_BERT.hh>
#include <G4UIsession.hh>
#include <G4TrackStatus.hh>

using namespace calin::air_shower::shower_generator;

namespace {



}
