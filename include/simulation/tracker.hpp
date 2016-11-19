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

#include<vector>
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

class TrackVisitor
{
public:
  virtual ~TrackVisitor();
  virtual void visit_event(const Event& event, bool& kill_event);
  virtual void visit_track(const Track& track, bool& kill_track);
  virtual void leave_event();
};

class VectorTrackVisitor: public TrackVisitor
{
public:
  VectorTrackVisitor() : TrackVisitor() { /* nothing to see here */ }
  virtual ~VectorTrackVisitor();
  virtual void visit_event(const Event& event, bool& kill_event);
  virtual void visit_track(const Track& track, bool& kill_track);
  virtual void leave_event();
  void add_visitor(TrackVisitor* visitor, bool adopt_visitor = false);
  void add_visitor_at_front(TrackVisitor* visitor, bool adopt_visitor = false);
protected:
  std::vector<TrackVisitor*> visitors_;
  std::vector<TrackVisitor*> adopted_visitors_;
};

class LengthLimitingTrackVisitor: public TrackVisitor
{
public:
  LengthLimitingTrackVisitor(TrackVisitor* visitor, double dx_max,
    double z_max = std::numeric_limits<double>::infinity(),
    bool adopt_visitor = false);
  virtual ~LengthLimitingTrackVisitor();
  void visit_event(const Event& event, bool& kill_event) override;
  void visit_track(const Track& track, bool& kill_track) override;
  void leave_event() override;
private:
  TrackVisitor* visitor_ = nullptr;
  bool adopt_visitor_ = false;
  double dx_max_ = 0;
  double z_max_ = std::numeric_limits<double>::infinity();
};

class TimeLimitingTrackVisitor: public TrackVisitor
{
public:
  TimeLimitingTrackVisitor(double t_max_sec = 60.0);
  virtual ~TimeLimitingTrackVisitor();
  void visit_event(const Event& event, bool& kill_event) override;
  void visit_track(const Track& track, bool& kill_track) override;
private:
  static uint64_t get_time();
  uint64_t t_event_max_ = 60000000000;
  uint64_t t_event_start_ = 0;
  bool do_kill_track_ = false;
  unsigned num_track_since_check_ = 0;
  unsigned num_track_before_check_ = 1000000;
};

class DebugStatsTrackVisitor: public TrackVisitor
{
public:
  DebugStatsTrackVisitor();
  virtual ~DebugStatsTrackVisitor();
  virtual void visit_event(const Event& event, bool& kill_event);
  virtual void visit_track(const Track& track, bool& kill_track);
  std::vector<int> particle_pdg_types() const;
  unsigned track_length_zero(int pdg_type) const;
  const calin::math::histogram::SimpleHist& track_length_log_hist(int pdg_type) const;
private:
  struct TrackStats
  {
    unsigned track_length_zero = { 0 };
    calin::math::histogram::SimpleHist track_length_hist = { 0.1 };
  };
  const TrackStats& find_stats(int pdg_type) const;

  std::map<int, TrackStats> stats_;
};

} } } // namespace calin::simulation::tracker
