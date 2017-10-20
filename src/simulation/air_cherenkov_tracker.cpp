/*

   calin/simulation/iact_array_tracker.cpp -- Stephen Fegan -- 2016-07-24

   Air shower track visitor to calcilate Cherenkov cone parameters.

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

   Certain portions are from code that is Copyright 2012, Stephen Fegan
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

// Certain fragments are from my EGS5 simulation code with following header:

// EGS5AtmosphericDetector.hpp - Detector class for EGS5 system that
// - implements a layered atmospheric detector
// Stephen Fegan - sfegan@llr.in2p3.fr - September 2012

#include <Eigen/Geometry>
#include <math/special.hpp>
#include <simulation/air_cherenkov_tracker.hpp>
#include <io/log.hpp>
#include <math/geometry.hpp>

using namespace calin::io::log;
using namespace calin::simulation::air_cherenkov_tracker;
using calin::math::special::SQR;

AirCherenkovTrackVisitor::~AirCherenkovTrackVisitor()
{
  // nothing to see here
}

void AirCherenkovTrackVisitor::visit_event(const Event& event, bool& kill_event)
{
  // default is to do nothing
}

void AirCherenkovTrackVisitor::
visit_cherenkov_track(const AirCherenkovTrack& track, bool& kill_track)
{
  // default is to do nothing
}

void AirCherenkovTrackVisitor::leave_event()
{
  // default is to do nothing
}

AirCherenkovParameterCalculatorTrackVisitor::
~AirCherenkovParameterCalculatorTrackVisitor()
{
  if(adopt_visitor_)delete visitor_;
  if(adopt_atm_)delete atm_;
}

void AirCherenkovParameterCalculatorTrackVisitor::
visit_event(const Event& event, bool& kill_event)
{
  visitor_->visit_event(event, kill_event);
}

void AirCherenkovParameterCalculatorTrackVisitor::
visit_track(const calin::simulation::tracker::Track& track, bool& kill_track)
{
  if(track.q == 0)return; // it's store policy: no charge = no radiation

  AirCherenkovTrack cherenkov;

  cherenkov.particle_track = &track;

  cherenkov.x0             = track.x0;
  cherenkov.e0             = track.e0;
  cherenkov.t0             = track.t0;

  cherenkov.dx_hat         = track.dx_hat;
  cherenkov.dx             = track.dx;
  cherenkov.de             = track.de;
  cherenkov.dt             = track.dt;

  cherenkov.x_mid          = cherenkov.x0 + 0.5*cherenkov.dx*cherenkov.dx_hat;
  cherenkov.e_mid          = cherenkov.e0 + 0.5*cherenkov.de;
  cherenkov.t_mid          = cherenkov.t0 + 0.5*cherenkov.dt;

  atm_->cherenkov_parameters(cherenkov.x_mid(2),
    cherenkov.n, cherenkov.propagation_delay);
  cherenkov.n             += 1.0;

  const double g2 = SQR(cherenkov.e_mid/track.mass); // gamma^2
  cherenkov.gamma_sq       = g2;
  const double b2 = 1.0 - 1.0/g2;                    // beta^2
  if(enable_forced_cherenkov_angle_mode_)
    cherenkov.sin2_thetac  = forced_sin2_thetac_;
  else
    cherenkov.sin2_thetac  = 1.0 - 1.0/(b2*SQR(cherenkov.n));
  if(cherenkov.sin2_thetac <= 0.0)return;
  cherenkov.yield_density  =
    YIELD_CONST*SQR(track.q)*cherenkov.sin2_thetac*cherenkov.dx;
  cherenkov.cos_thetac     = std::sqrt(1.0 - cherenkov.sin2_thetac);
  cherenkov.sin_thetac     = std::sqrt(cherenkov.sin2_thetac);

  visitor_->visit_cherenkov_track(cherenkov, kill_track);
}

void AirCherenkovParameterCalculatorTrackVisitor::leave_event()
{
  visitor_->leave_event();
}

void AirCherenkovParameterCalculatorTrackVisitor::
set_forced_cherenkov_angle(double theta_c_deg)
{
  if(enable_forced_cherenkov_angle_mode_)
    forced_sin2_thetac_ = SQR(std::sin(theta_c_deg/180.0*M_PI));
  else
    throw std::runtime_error("Forced cherenkov angle mode not enabled.");
}

calin::ix::simulation::tracker::AirCherenkovParameterCalculatorTrackVisitorConfig
AirCherenkovParameterCalculatorTrackVisitor::default_config()
{
  calin::ix::simulation::tracker::AirCherenkovParameterCalculatorTrackVisitorConfig cfg;
  cfg.set_forced_cherenkov_angle(1.2);
  return cfg;
}

CherenkovPhotonVisitor::~CherenkovPhotonVisitor()
{
  // nothing to see here
}

void CherenkovPhotonVisitor::visit_event(const Event& event, bool& kill_event)
{
  // nothing to see here
}

void CherenkovPhotonVisitor::
visit_cherenkov_photon(const CherenkovPhoton& cherenkov_photon)
{
  // nothing to see here
}

void CherenkovPhotonVisitor::leave_event()
{
  // nothing to see here
}

MCCherenkovPhotonGenerator::
MCCherenkovPhotonGenerator(CherenkovPhotonVisitor* visitor,
    double epsilon0, double bandwidth, bool do_color_photons,
    calin::math::rng::RNG* rng, bool adopt_visitor, bool adopt_rng):
  visitor_(visitor), adopt_visitor_(adopt_visitor),
  epsilon0_(epsilon0), bandwidth_(bandwidth), do_color_photons_(do_color_photons),
  rng_(rng ? rng : new calin::math::rng::RNG()), adopt_rng_(rng ? adopt_rng : true)
{
  dX_emission_ = rng_->exponential();
}

MCCherenkovPhotonGenerator::~MCCherenkovPhotonGenerator()
{
  if(adopt_visitor_)delete visitor_;
  if(adopt_rng_)delete rng_;
}

void MCCherenkovPhotonGenerator::visit_event(const Event& event, bool& kill_event)
{
  return visitor_->visit_event(event, kill_event);
}

void MCCherenkovPhotonGenerator::leave_event()
{
  return visitor_->leave_event();
}

void MCCherenkovPhotonGenerator::
visit_cherenkov_track(const AirCherenkovTrack& cherenkov_track, bool& kill_track)
{
  CherenkovPhoton photon;
  bool no_photon_emitted = true;

  const double dX_track = cherenkov_track.yield_density * bandwidth_;

  Eigen::Matrix3d Mrot;
  Eigen::Vector3d u0;

  double dX_left = dX_track;
  while(dX_emission_ < dX_left)
  {
    dX_left -= dX_emission_;
    dX_emission_ = rng_->exponential();

    if(no_photon_emitted)
    {
      no_photon_emitted = false;
      photon.air_cherenkov_track = &cherenkov_track;
      calin::math::geometry::rotation_z_to_vec(Mrot, cherenkov_track.dx_hat);
      photon.epsilon = 0;
    }

    double dX_frac = 1.0 - dX_left/dX_track;
    photon.x0 = cherenkov_track.x0 + cherenkov_track.dx_hat*cherenkov_track.dx*dX_frac;
    double phi = 2.0*M_PI*rng_->uniform();
    u0 << cherenkov_track.sin_thetac*std::cos(phi), cherenkov_track.sin_thetac*std::sin(phi), cherenkov_track.cos_thetac;
    photon.u0 = Mrot * u0;
    photon.t0 = cherenkov_track.t0 + cherenkov_track.dt*dX_frac;
    visitor_->visit_cherenkov_photon(photon);
    if(do_color_photons_)photon.epsilon = epsilon0_ + rng_->uniform()*bandwidth_;
  }
  dX_emission_ -= dX_left;
}
