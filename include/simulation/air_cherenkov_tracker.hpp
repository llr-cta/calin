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

#include"simulation/tracker.hpp"

namespace calin { namespace simulation { namespace air_cherenkov_tracker {

using calin::simulation::tracker::Event;

struct AirCherenkovTrack
{
  double sin2_thetac;      // Sine^2 of emission angle     [1]
  double sin_thetac;       // Sine of emssion angle        [1]
  double cos_thetac;       // Cosine of emission angle     [1]
  double yield_density;    // Cherenkov photon density     [ph/eV]
  double n;                // Refracive index              [1]
  double gamma;            // Particle gamma               [1]

  Eigen::Vector3d x0;      // Position of start of track   [cm]
  double e0;               // Total energy at start of trk [MeV]
  double t0;               // Time at start of track       [ns]

  Eigen::Vector3d dx_hat;  // Unit vector from x0 to x1    [1]
  double dx;               // Step length                  [cm]
  double de;               // Change in energy             [MeV]
  double dt;               // Time step                    [ns]

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
  virtual void visit_cherenkov_track(const AirCherenkovTrack& track, bool& kill_track);
  virtual void leave_event();
};

class AirCherenkovParametersTrackVisitor:
  public calin::simulation::tracker::TrackVisitor
{
public:
  virtual ~AirCherenkovParametersTrackVisitor();
  void visit_event(const Event& event, bool& kill_event) override;
  void visit_track(const calin::simulation::tracker::Track& track, bool& kill_track) override;
  void leave_event() override;
private:
  AirCherenkovTrackVisitor* visitor_ = nullptr;
  bool adopt_visitor_ = false;

};

} } } // namespace calin::simulation::air_cherenkov_tracker
