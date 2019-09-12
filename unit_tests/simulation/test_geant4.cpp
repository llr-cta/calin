/*

   calin/unit_tests/simulation/test_geant4.cpp -- Stephen Fegan -- 2015-10-06

   Unit tests for Geant4 simulation code

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

#include <gtest/gtest.h>

#include "simulation/atmosphere.hpp"
#include "simulation/tracker.hpp"
#include "simulation/geant4_shower_generator.hpp"

using namespace calin::simulation::atmosphere;
using namespace calin::simulation::tracker;
using namespace calin::simulation::geant4_shower_generator;

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
                            100, 0, atm->top_of_atmosphere(), nullptr,
                            //VerbosityLevel::SUPRESSED_STDOUT);
                            //VerbosityLevel::NORMAL);
                            VerbosityLevel::VERBOSE_EVERYTHING);
  sim.set_minimum_energy_cut(20); // MeV
  sim.generate_showers(100, ParticleType::MUON, 1e6,
                      Eigen::Vector3d(0, 0, atm->top_of_atmosphere()));
  delete(visitor);
  delete(atm);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
