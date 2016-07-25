/*

   calin/simulation/iact_array_tracker.hpp -- Stephen Fegan -- 2016-07-22

   Air shower track visitor with array of multiple IACT detector spheres

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

#include"simulation/air_cherenkov_tracker.hpp"

namespace calin { namespace simulation { namespace iact_array_tracker {


struct IACTDetectorSphereHit
{
  Eigen::Vector3d x0;    // Position of particle at emisson
  Eigen::Vector3d u;     // Unit vector along motion of particle at emisson
  //Eigen::Vector3d v;     // Unit vector prependicular to u from x0 to r0
  double rxu;            //
  double rxv;            //
  double dmin;           //
  double cos_phimax;     //
  double phimax;         //
};

class IACTDetectorSphereHitProcessor
{
public:
  virtual ~IACTDetectorSphereHitProcessor();
  virtual void process_hit(const IACTDetectorSphereHit& hit);
};

struct IACTDetectorSphere
{
  Eigen::Vector3d r0;                        // Center of detector sphere [cm]
  double radius_sq;                          // Squared radius of sphere  [cm^2]
  IACTDetectorSphereHitProcessor* processor; // Hit processor for this detector
};

class HitIACTVisitor:
  public calin::simulation::tracker::TrackVisitor
{
public:
  virtual ~HitIACTVisitor();
  virtual std::vector<IACTDetectorSphere> spheres();
  virtual void visit_event(const calin::simulation::tracker::Event& event,
    bool& kill_event);
  virtual void visit_cherenkov_track(
    const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack& cherenkov_track,
    bool& kill_track);
  virtual void leave_cherenkov_track();
  virtual void leave_event();
};

class IACTDetectorSphereAirCherenkovTrackVisitor:
  public calin::simulation::air_cherenkov_tracker::AirCherenkovTrackVisitor
{
public:
  IACTDetectorSphereAirCherenkovTrackVisitor(HitIACTVisitor* visitor_,
      bool adopt_visitor = false);
  virtual ~IACTDetectorSphereAirCherenkovTrackVisitor();
  virtual void visit_event(const calin::simulation::tracker::Event& event,
    bool& kill_event);
  virtual void visit_cherenkov_track(
    const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack& cherenkov_track,
    bool& kill_track);
  virtual void leave_event();
private:
  HitIACTVisitor* visitor_;
  bool adopt_visitor_;
  std::vector<IACTDetectorSphere> spheres_;
};

} } } // namespace calin::simulation::iact_array_tracker
