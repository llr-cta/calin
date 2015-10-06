#include <gtest/gtest.h>

#include "simulation/atmosphere.hpp"
#include "simulation/tracker.hpp"
#include "simulation/geant4_shower_generator.hpp"

using namespace calin::simulation::atmosphere;
using namespace calin::simulation::tracker;
using namespace calin::simulation::shower_generator;

TEST(TestGeant4, MakeUS76Atmosphere) {
  Atmosphere* atm = LayeredAtmosphere::us76();
  ASSERT_EQ(atm->thickness(0),1035);
  delete(atm);
}

TEST(TestGeant4, MakeGeant4Simulator) {
  //CLHEP::HepRandom::setTheSeed(time(0));
  Atmosphere* atm = LayeredAtmosphere::us76();
  TrackVisitor* visitor = new TrackVisitor;
  Geant4ShowerGenerator sim(visitor, atm,
                            100, 0, atm->top_of_atmosphere(),
                            //VerbosityLevel::SUPRESSED_STDOUT);
                            //VerbosityLevel::NORMAL);
                            VerbosityLevel::VERBOSE_EVERYTHING);
  sim.setMinimumEnergyCut(20); // MeV
  sim.generateShowers(100, ParticleType::MUON, 1e6,
                      Eigen::Vector3d(0, 0, atm->top_of_atmosphere()));
  delete(visitor);
  delete(atm);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
