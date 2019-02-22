/*

   calin/simulation/tracker.hpp -- Stephen Fegan -- 2015-06-23

   Base class for all air shower track visitors

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

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

#include<vector>
#include<ostream>
#include<map>

#include"Eigen/Core"
#include"math/histogram.hpp"

namespace calin { namespace simulation { namespace tracker {

enum class ParticleType { GAMMA, ELECTRON, POSITRON, MUON, ANTI_MUON,
    PROTON, ANTI_PROTON, OTHER };

ParticleType pdg_type_to_particle_type(int pdg_type);
int particle_type_to_pdg_type(ParticleType track_type);
double particle_type_to_mass(ParticleType track_type);
double particle_type_to_charge(ParticleType track_type);
std::string particle_type_to_string(ParticleType track_type);

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

  double weight;         // Event weighting for thinning
};

struct Track
{
  ParticleType type;       // Simplified particle type
  int pdg_type;            // PDG particle type code
  double q;                // PDG particle charge          [e]
  double mass;             // PDG particle rest mass       [MeV]

  Eigen::Vector3d x0;      // Position of start of track   [cm]
  Eigen::Vector3d u0;      // Direction at start of track  [1]
  double e0;               // Total energy at start of trk [MeV]
  double t0;               // Time at start of track       [ns]

  Eigen::Vector3d x1;      // Position of end of track     [cm]
  Eigen::Vector3d u1;      // Direction at end of track    [1]
  double e1;               // Total energy at end of track [MeV]
  double t1;               // Time at end of track         [ns]

  Eigen::Vector3d dx_hat;  // Unit vector from x0 to x1    [1]
  double dx;               // Step length                  [cm]
  double de;               // Change in energy             [MeV]
  double dt;               // Time step                    [ns]

  double weight;           // Track weighting for thinning
};

std::ostream& operator<<(std::ostream& stream, const Track& t);

class TrackVisitor
{
public:
  virtual ~TrackVisitor();
  virtual void visit_event(const Event& event, bool& kill_event);
  virtual void visit_track(const Track& track, bool& kill_track);
  virtual void leave_event();
};

class MultiDelegatingTrackVisitor: public TrackVisitor
{
public:
  MultiDelegatingTrackVisitor() : TrackVisitor() { /* nothing to see here */ }
  virtual ~MultiDelegatingTrackVisitor();
  virtual void visit_event(const Event& event, bool& kill_event);
  virtual void visit_track(const Track& track, bool& kill_track);
  virtual void leave_event();
  void add_visitor(TrackVisitor* visitor, bool adopt_visitor = false);
  void add_visitor_at_front(TrackVisitor* visitor, bool adopt_visitor = false);
protected:
  std::vector<TrackVisitor*> visitors_;
  std::vector<TrackVisitor*> adopted_visitors_;
};

} } } // namespace calin::simulation::tracker
