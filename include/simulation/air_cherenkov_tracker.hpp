/*

   calin/simulation/air_cherenkov_tracker.hpp -- Stephen Fegan -- 2016-07-24

   Air shower track visitor to calcilate Cherenkov cone parameters

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

#pragma once

// Units:
// Time, t:      ns
// Height, z:    cm
// Energy, e:    MeV

#include"Eigen/Core"

#include"simulation/atmosphere.hpp"
#include"simulation/tracker.hpp"

namespace calin { namespace simulation { namespace air_cherenkov_tracker {

using calin::simulation::tracker::Event;

// Cherenkov yield in photons per track length per unit bandwidth [eV^-1 cm^-1]
// alpha^2 / r_e / m_e / c^2 -- see PDG 2014, eqn. 32.45
constexpr double YIELD_CONST = 369.81020849958;

struct AirCherenkovTrack
{
  double sin2_thetac;       // Sine^2 of emission angle          [1]
  double sin_thetac;        // Sine of emssion angle             [1]
  double cos_thetac;        // Cosine of emission angle          [1]
  double yield_density;     // Cherenkov photon density          [ph/eV]
  double n;                 // Refracive index                   [1]
  double propagation_delay; // Vertical propagation delay in ct  [cm]
  double gamma_sq;          // Particle gamma squared            [1]

  Eigen::Vector3d x0;       // Position of start of track        [cm]
  double e0;                // Total energy at start of trk      [MeV]
  double t0;                // Time at start of track            [ns]

  Eigen::Vector3d x_mid;    // Position of "middle" of track     [cm]
  double e_mid;             // Total energy at "middle" of trk   [MeV]
  double t_mid;             // Time at "middle" of track         [ns]

  Eigen::Vector3d dx_hat;   // Unit vector from x0 to x1         [1]
  double dx;                // Step length                       [cm]
  double de;                // Change in energy                  [MeV]
  double dt;                // Time step                         [ns]

  // Track from underlying simulator - this could be split by AirCherenkov
  // code if it is too long, i.e. if theta_c changes significantly along
  // the track
  const calin::simulation::tracker::Track* particle_track;
};

class AirCherenkovTrackVisitor
{
public:
  virtual ~AirCherenkovTrackVisitor();
  virtual void visit_event(const Event& event, bool& kill_event);
  virtual void visit_cherenkov_track(const AirCherenkovTrack& cherenkov_track,
    bool& kill_track);
  virtual void leave_event();
};

class AirCherenkovParameterCalculatorTrackVisitor:
  public calin::simulation::tracker::TrackVisitor
{
public:
  AirCherenkovParameterCalculatorTrackVisitor(AirCherenkovTrackVisitor* visitor,
      calin::simulation::atmosphere::Atmosphere* atm,
      bool adopt_visitor = false, bool adopt_atm = false):
    calin::simulation::tracker::TrackVisitor(),
    visitor_(visitor), adopt_visitor_(adopt_visitor),
    atm_(atm), adopt_atm_(adopt_atm) { /* nothing to see here */ }
  virtual ~AirCherenkovParameterCalculatorTrackVisitor();
  void visit_event(const Event& event, bool& kill_event) override;
  void visit_track(const calin::simulation::tracker::Track& track,
    bool& kill_track) override;
  void leave_event() override;
private:
  AirCherenkovTrackVisitor* visitor_ = nullptr;
  bool adopt_visitor_ = false;
  calin::simulation::atmosphere::Atmosphere* atm_ = nullptr;
  bool adopt_atm_ = false;
};

} } } // namespace calin::simulation::air_cherenkov_tracker
