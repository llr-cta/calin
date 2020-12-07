/*

   calin/simulation/bfield_track_generator.cpp -- Stephen Fegan -- 2017-06-16

   Class to genereate curved tracks in a magnetic field (neglecting
   synchrotron emission)

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <Eigen/Geometry>
#include <simulation/bfield_track_generator.hpp>
#include <util/log.hpp>
#include <math/constants.hpp>
#include <math/special.hpp>

using namespace calin::simulation::bfield_track_generator;
using calin::math::special::SQR;

namespace {
  constexpr double GYRO_CONST = 89.8755178736818; // e*nanotesla/(MeV/c^2) [1/s]
}

BFieldTrackGenerator::
BFieldTrackGenerator(const Eigen::Vector3d& bfield_nT,
    double zground_or_dist, double step_size,
    PropagationMode propagation_mode):
  bfield_(bfield_nT),
  propagation_mode_(propagation_mode), zground_or_dist_(zground_or_dist),
  step_size_(std::abs(step_size))
{
  // nothing to see here
}

BFieldTrackGenerator::
BFieldTrackGenerator(double zground_or_dist, double step_size,
    PropagationMode propagation_mode):
  bfield_(Eigen::Vector3d::Zero()),
  propagation_mode_(propagation_mode), zground_or_dist_(zground_or_dist),
  step_size_(std::abs(step_size))
{
  // nothing to see here
}

BFieldTrackGenerator::~BFieldTrackGenerator()
{
  // nothing to see here
}

void BFieldTrackGenerator::generate_showers(calin::simulation::tracker::TrackVisitor* visitor,
  unsigned num_events, calin::simulation::tracker::ParticleType type, double total_energy,
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
  track.weight   = event.weight;
  track.e0       = event.e0;
  track.e1       = event.e0;
  track.de       = 0.0;

  const double gamma = track.e0/track.mass;
  assert(gamma >= 1);
  const double v = std::sqrt(1.0-1.0/SQR(gamma))*calin::math::constants::cgs_c;
  const double dt = step_size_/v;
  const double gyro = GYRO_CONST * track.q / track.e0 * dt;

  track.dx       = step_size_;
  track.dt       = dt;

  while(num_events--)
  {
    bool kill_event = false;
    event.event_id = event_id_++;
    visitor->visit_event(event, kill_event);
    if(kill_event)continue;

    if(propagation_mode_ == FWD_TO_GROUND)
    {
      track.x1       = event.x0;
      track.u1       = event.u0;
      track.t1       = event.t0;

      Eigen::Vector3d umid = u0 + 0.5*gyro*u0.cross(bfield_);
      umid.normalize();

      bool kill_track = false;
      while(track.x1.z() > zground_or_dist_ and not kill_track)
      {
        track.x0     = track.x1;
        track.u0     = track.u1;
        track.t0     = track.t1;

        track.x1     = track.x0 + umid * step_size_;

        if(track.x1.z() < zground_or_dist_)
        {
          double step_size = (zground_or_dist_ - track.x0.z())/umid.z();
          track.x1.x() = track.x0.x() + umid.x() * step_size;
          track.x1.y() = track.x0.y() + umid.y() * step_size;
          track.x1.z() = zground_or_dist_;
          track.dx     = step_size;
          track.dt     = step_size/v;
        } else {
          track.dx     = step_size_;
          track.dt     = dt;
        }


        track.t1     = track.t0 + dt;

        track.dx_hat = umid;

        Eigen::Vector3d umid_new = umid + gyro*umid.cross(bfield_);
        umid_new.normalize();

        track.u1     = 0.5*(umid + umid_new);
        track.u1.normalize();

        umid = umid_new;

        visitor->visit_track(track, kill_track);
      }
    }
    else // propagation_mode_ == BWD_FIXED_DISTANCE
    {
      track.x0       = event.x0;
      track.u0       = event.u0;
      track.t0       = event.t0;

      Eigen::Vector3d umid = u0 - 0.5*gyro*u0.cross(bfield_);
      umid.normalize();

      double dist = 0;
      bool kill_track = false;
      while(dist < zground_or_dist_ and not kill_track)
      {
        track.x1     = track.x0;
        track.u1     = track.u0;
        track.t1     = track.t0;

        track.x0     = track.x1 - umid * step_size_;
        track.t0     = track.t1 - dt;

        track.dx_hat = umid;

        Eigen::Vector3d umid_new = umid - gyro*umid.cross(bfield_);
        umid_new.normalize();

        track.u0     = 0.5*(umid + umid_new);
        track.u0.normalize();

        umid = umid_new;

        dist += step_size_;

        visitor->visit_track(track, kill_track);
      }

    }

    visitor->leave_event();
  }
}
