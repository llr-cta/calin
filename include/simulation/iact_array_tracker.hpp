/*

   calin/simulation/iact_array_tracker.hpp -- Stephen Fegan -- 2016-07-22

   Air shower track visitor with array of multiple IACT detector spheres

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include"Eigen/Core"

#include"simulation/air_cherenkov_tracker.hpp"

namespace calin { namespace simulation { namespace iact_array_tracker {

struct IACTDetectorSphereHit
{
  Eigen::Vector3d x0;    // Position of particle at emisson
  Eigen::Vector3d rx;    // Vector from x0 to r0 (r0-x0)
  Eigen::Vector3d u;     // Unit vector along motion of particle at emisson
  Eigen::Vector3d v;     // Unit vector perpendicular to u from x0 to r0
  Eigen::Vector3d w;     // Third unit vector such that u cross v = w
  Eigen::Matrix3d rot;   // Matrix that maps x,y,z onto u,v,w
  double rx2;            // (r0-x0)^2
  double rxu;            // Component of r0-x0 along u
  double rxv;            // Component of r0-x0 along v
  double dmin;           // Distance of closest approach between cone and sphere
  double cos_phimax;     // Cosine of half-angle between cone and sphere
  double phimax;         // Half-angle between cone and sphere

  const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack* cherenkov_track;
};

class IACTDetectorSphereHitProcessor
{
public:
  virtual ~IACTDetectorSphereHitProcessor();
  virtual void process_hit(const IACTDetectorSphereHit& hit) = 0;
};

struct IACTDetectorSphere
{
  IACTDetectorSphere(const Eigen::Vector3d& r0_, double radius_sq_,
    IACTDetectorSphereHitProcessor* processor_):
      r0(r0_), radius_sq(radius_sq_), processor(processor_) {
    /* nothing to see here */ }
  Eigen::Vector3d r0;                        // Center of detector sphere [cm]
  double radius_sq;                          // Squared radius of sphere  [cm^2]
  IACTDetectorSphereHitProcessor* processor; // Hit processor for this detector
};

class HitIACTVisitor
{
public:
  virtual ~HitIACTVisitor();
  virtual std::vector<IACTDetectorSphere> spheres() = 0;
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
