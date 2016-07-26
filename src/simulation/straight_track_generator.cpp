/*

   calin/simulation/straight_track_generator.cpp -- Stephen Fegan -- 2016-07-26

   Class to genereate straight tracks

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <simulation/straight_track_generator.hpp>
#include <io/log.hpp>

using namespace calin::simulation::straight_track_generator;

StraightTrackGenerator::
StraightTrackGenerator(calin::simulation::tracker::TrackVisitor* visitor,
    double zground, bool adopt_visitor):
  visitor_(visitor), adopt_visitor_(adopt_visitor), zground_(zground)
{
  // nothing to see here
}

StraightTrackGenerator::~StraightTrackGenerator()
{
  if(adopt_visitor_)delete visitor_;
}

void StraightTrackGenerator::generate_showers(unsigned num_events,
  calin::simulation::tracker::ParticleType type, double total_energy,
  const Eigen::Vector3d& x0, const Eigen::Vector3d& u0, double weight)
{

}
