/*

   calin/simulation/tracker.hpp -- Stephen Fegan -- 2015-06-23

   Base class for all air shower track visitors

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

namespace calin { namespace simulation { namespace tracker {

enum class ParticleType { GAMMA, ELECTRON, POSITRON, MUON, ANTI_MUON,
    PROTON, ANTI_PROTON, OTHER };

struct Event
{
  int event_id;
  
  ParticleType type;     // Simplified particle type
  int pdg_type;          // PDG particle type code
  double q;              // PDG particle charge          [e]
  double mass;           // PDG particle rest mass       [MeV]

  Eigen::Vector3d x0;    // Position of start of event   [cm]
  Eigen::Vector3d u0;    // Direction of start of event  [1]
  double e0;             // Total energy at start of evt [MeV]
  double t0;             // Time at start of event       [ns]

  double weight;         // Track weighting for thinning
};

struct Track
{
  ParticleType type;     // Simplified particle type
  int pdg_type;          // PDG particle type code
  double q;              // PDG particle charge          [e]
  double mass;           // PDG particle rest mass       [MeV]

  Eigen::Vector3d x0;    // Position of start of track   [cm]
  Eigen::Vector3d u0;    // Direction of start of track  [1]
  double e0;             // Total energy at start of trk [MeV]
  double t0;             // Time at start of track       [ns]

  Eigen::Vector3d x1;    // Position of end of track     [cm]
  Eigen::Vector3d u1;    // Direction of end of track    [1]
  double e1;             // Total energy at end of track [MeV]
  double t1;             // Time at end of track         [ns]

  double dx;             // Step length                  [cm]

  double weight;         // Track weighting for thinning
};

class TrackVisitor
{
 public:
  virtual ~TrackVisitor();
  virtual void visit_event(const Event& event, bool& kill_event);
  virtual void visit_track(const Track& track, bool& kill_track);
  virtual void leave_event();
};

} } } // namespace calin::simulation::tracker
