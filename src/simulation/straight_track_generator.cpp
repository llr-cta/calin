/*

   calin/simulation/straight_track_generator.cpp -- Stephen Fegan -- 2016-07-26

   Class to genereate straight tracks

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <simulation/straight_track_generator.hpp>
#include <util/log.hpp>
#include <math/constants.hpp>
#include <math/special.hpp>

using namespace calin::simulation::straight_track_generator;
using calin::math::special::SQR;

StraightTrackGenerator::
StraightTrackGenerator(double zground):
  calin::simulation::tracker::ShowerGenerator(), zground_(zground)
{
  // nothing to see here
}

StraightTrackGenerator::~StraightTrackGenerator()
{
  // nothing to see here
}

void StraightTrackGenerator::generate_showers(
  calin::simulation::tracker::TrackVisitor* visitor, unsigned num_events,
  calin::simulation::tracker::ParticleType type, double total_energy,
  const Eigen::Vector3d& x0, const Eigen::Vector3d& u0, double weight)
{
  tracker::Event event;
  event.type     = type;
  event.pdg_type = calin::simulation::tracker::particle_type_to_pdg_type(type);
  event.q        = calin::simulation::tracker::particle_type_to_charge(type);
  event.mass     = calin::simulation::tracker::particle_type_to_mass(type);
  event.x0       = x0;
  event.u0       = u0;
  event.u0.normalize();
  event.e0       = total_energy;
  event.t0       = 0.0;
  event.weight   = weight;

  tracker::Track track;
  track.type     = event.type;
  track.pdg_type = event.pdg_type;
  track.q        = event.q;
  track.mass     = event.mass;

  track.x0       = event.x0;
  track.u0       = event.u0;
  track.e0       = event.e0;
  track.t0       = event.t0;

  const double dz = zground_ - event.x0(2);
  const double dx = dz/track.u0(2);
  const double gamma = track.e0/track.mass;
  assert(gamma >= 1);
  const double v = std::sqrt(1.0-1.0/SQR(gamma))*calin::math::constants::cgs_c;
  const double dt = dx/v;

  track.x1       = track.x0 + track.u0 * dx;
  track.u1       = track.u0;
  track.e1       = track.e0;
  track.t1       = track.t0 + dt;

  track.dx_hat   = track.u0;
  track.dx       = dx;
  track.de       = 0.0;
  track.dt       = dt;

  track.weight   = event.weight;

  while(num_events--)
  {
    bool kill_event = false;
    event.event_id = event_id_++;
    visitor->visit_event(event, kill_event);
    if(kill_event)continue;
    bool kill_track = false;
    visitor->visit_track(track, kill_track);
    visitor->leave_event();
  }
}
