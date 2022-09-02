/*

   calin/unit_tests/simulation/speedtest_vcl_iact.cpp -- Stephen Fegan -- 2019-02-24

   Performance tests for VCL IACT code

   Copyright 2019, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include "simulation/geant4_shower_generator.hpp"
#include "simulation/vcl_iact_ground_map.hpp"
#include "simulation/atmosphere.hpp"
#include "provenance/system_info.hpp"
#include "util/vcl.hpp"
#include "util/string.hpp"

char** global_argv;
int global_argc;

using namespace calin::simulation::atmosphere;
using namespace calin::simulation::tracker;
using namespace calin::simulation::geant4_shower_generator;

TEST(SpeedTestVCLIACT, Generate100Protons1TeV) {
  //CLHEP::HepRandom::setTheSeed(time(0));
  std::string datadir = calin::provenance::system_info::build_info()->data_install_dir();
  auto* atm = new LayeredRefractiveAtmosphere(datadir + "/simulation/atmprof36.dat");

  auto config = calin::simulation::vcl_iact::VCLIACTGroundMap<
    calin::util::vcl::VCL256Architecture>::default_config();
  if(global_argc > 1) {
    config.set_bandwidth(calin::util::string::double_from_string(global_argv[1]));
    std::cerr << "HELLO: " << global_argv[1] << ' ' << config.bandwidth() << '\n';
  }
  auto* act = new calin::simulation::vcl_iact::VCLIACTGroundMap<
    calin::util::vcl::VCL256Architecture>(atm, config);

  Geant4ShowerGenerator sim(atm,
                            1000, 0, atm->top_of_atmosphere(), nullptr,
                            //VerbosityLevel::SUPRESSED_STDOUT);
                            VerbosityLevel::NORMAL);
                            //VerbosityLevel::VERBOSE_EVERYTHING, 0, 10.0);
  sim.set_minimum_energy_cut(20); // MeV

  uint64_t ntrack = 0;
  uint64_t nstep = 0;
  uint64_t nray = 0;

  for(unsigned i=0; i<100; i++) {
    sim.generate_showers(act, 1, ParticleType::PROTON, 1000E3,
                        Eigen::Vector3d(0, 0, atm->top_of_atmosphere()),
                        Eigen::Vector3d(0, 0, -1.0));
    ntrack += act->num_tracks();
    nstep += act->num_steps();
    nray += act->num_rays();
  }

  std::cout
    << "Ntrack: " << ntrack << '\n'
    << "Nsteps: " << nstep << "  (" << double(nstep)/double(ntrack) << "/track)\n"
    << "Nrays:  " << nray << "  (" << double(nray)/double(ntrack) << "/track, "
      << double(nray)/double(nstep) << "/step)\n";

  delete(act);
  delete(atm);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  global_argc = argc;
  global_argv = argv;
  return RUN_ALL_TESTS();
}
