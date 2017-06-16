/*

   calin/simulation/bfield_track_generator.cpp -- Stephen Fegan -- 2017-06-16

   Class to genereate curved tracks in a magnetic field (neglecting
   synchrotron emission)

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <simulation/bfield_track_generator.hpp>
#include <io/log.hpp>
#include <math/constants.hpp>
#include <math/special.hpp>

using namespace calin::simulation::bfield_track_generator;
using calin::math::special::SQR;

namespace {
  constexpr double GYRO_CONST = 89.8755178736818; // e*nanotesla/(MeV/c^2) [1/s]
}

BFieldTrackGenerator::
BFieldTrackGenerator(calin::simulation::tracker::TrackVisitor* visitor,
    const Eigen::Vector3d& bfield_nT, double zground_or_dist, double step_size,
    PropagationMode propagation_mode, bool adopt_visitor):
  visitor_(visitor), adopt_visitor_(adopt_visitor), bfield_(bfield_nT),
  propagation_mode_(propagation_mode), zground_or_dist_(zground_or_dist),
  step_size_(step_size)
{
  // nothing to see here
}

BFieldTrackGenerator::~BFieldTrackGenerator()
{
  if(adopt_visitor_)delete visitor_;
}

void BFieldTrackGenerator::generate_showers(unsigned num_events,
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

  const double gamma = track.e0/track.mass;
  assert(gamma >= 1);
  const double v = std::sqrt(1.0-1.0/SQR(gamma))*calin::math::constants::cgs_c;
  const double dt = step_size_/v;
  const double gyro = GYRO_CONST * track.q / gamma / track.mass * dt;

  tracker::Track track;
  track.type     = event.type;
  track.pdg_type = event.pdg_type;
  track.q        = event.q;
  track.mass     = event.mass;
  track.weight   = event.weight;
  track.e0       = event.e0;
  track.e1       = event.e0;
  track.dx       = step_size_;
  track.de       = 0.0;
  track.dt       = dt;

  while(num_events--)
  {
    bool kill_event = false;
    event.event_id = event_id_++;
    visitor_->visit_event(event, kill_event);
    if(kill_event)continue;

    Eigen::Vector3d x = x0;
    if(propagation_mode_ == FWD_TO_GROUND)
    {
      track.x1       = event.x0;
      track.u1       = event.u0;
      track.t1       = event.t0;

      Eigen::Vector3d umid = u0 + 0.5*gyro*u0.cross(bfield_))
      umid.normalize();

      while(x.z() > zground_or_dist_)
      {
        track.x0     = track.x1
        track.u0     = track.u1;
        track.t0     = track.t1;

        track.x1     = track.x0 + umid * step_size_;
        track.t1     = track.t0 + dt;

        track.dx_hat = umid;

        Eigen::Vector3d umid_new = umid + gyro*umid.cross(bfield_));
        umid_new.normalize();

        track.u1     = 0.5*(umid + umid_new);
        track.u1.normalize();

        umid = umid_new;
      }
    }
    else
    {
      track.x0       = event.x0;
      track.u0       = event.u0;
      track.t0       = event.t0;

      Eigen::Vector3d umid = u0 - 0.5*gyro*u0.cross(bfield_))
      umid.normalize();

      double dist = 0;
      while(dist < zground_or_dist_)
      {
        track.x1     = track.x0
        track.u1     = track.u0;
        track.t1     = track.t0;

        track.x0     = track.x1 - umid * step_size_;
        track.t0     = track.t1 - dt;

        track.dx_hat = umid;

        Eigen::Vector3d umid_new = umid - gyro*umid.cross(bfield_));
        umid_new.normalize();

        track.u0     = 0.5*(umid + umid_new);
        track.u0.normalize();

        umid = umid_new;

        dist += step_size;
      }

    }

    bool kill_track = false;
    visitor_->visit_track(track, kill_track);
    visitor_->leave_event();
  }
}
