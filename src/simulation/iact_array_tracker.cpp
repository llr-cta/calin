/*

   calin/simulation/iact_array_tracker.cpp -- Stephen Fegan -- 2016-07-24

   Base class for all air shower track visitors

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

#include <limits>

#include <simulation/iact_array_tracker.hpp>
#include <math/special.hpp>

using namespace calin::simulation::iact_array_tracker;
using calin::math::special::SQR;

IACTDetectorSphereAirCherenkovTrackVisitor::
IACTDetectorSphereAirCherenkovTrackVisitor(HitIACTVisitor* visitor,
    bool adopt_visitor):
  calin::simulation::air_cherenkov_tracker::AirCherenkovTrackVisitor(),
  visitor_(visitor), adopt_visitor_(adopt_visitor),
  spheres_(visitor->spheres())
{
  // nothing to see here
}

IACTDetectorSphereAirCherenkovTrackVisitor::
~IACTDetectorSphereAirCherenkovTrackVisitor()
{
  for(auto& isphere : spheres_)delete isphere.processor;
  if(adopt_visitor_)delete visitor_;
}

void IACTDetectorSphereAirCherenkovTrackVisitor::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
  visitor_->visit_event(event, kill_event);
}

void IACTDetectorSphereAirCherenkovTrackVisitor::visit_cherenkov_track(
  const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack& CT,
  bool& kill_track)
{
  visitor_->visit_cherenkov_track(CT, kill_track);
  if(kill_track)return;

  for(auto& isphere : spheres_)
  {
    Eigen::Vector3d rx = isphere.r0 - CT.x_mid;
    double rx2  = rx.squaredNorm();
    double rxu  = rx.dot(CT.dx_hat); // note bad naming of dx_hat
    double rxu2 = SQR(rxu);
    double rxv = std::sqrt(std::max(rx2 - rxu2, 0.0));

    double rx2_minus_r2 = rx2 - isphere.radius_sq;

    if(rx2_minus_r2 < 0) // Particle is inside detector sphere - always hit!
    {
      IACTDetectorSphereHit hit;
      hit.x0         = CT.x_mid;
      hit.u          = CT.dx_hat;
      hit.rxu        = rxu;
      hit.rxv        = rxv;
      hit.dmin       = std::numeric_limits<double>::quiet_NaN();
      hit.cos_phimax = -1;
      hit.phimax     = M_PI;
      isphere.processor->process_hit(hit);
      continue;
    }

    double dmin = CT.cos_thetac*rxv - CT.sin_thetac*rxu;
    if(SQR(dmin) > isphere.radius_sq)continue; // Cone misses sphere

    double cos_phimax = (sqrt(rx2_minus_r2) - CT.cos_thetac*rxu) /
      (CT.sin_thetac*rxv);

    if(cos_phimax < 1.0)
    {
      IACTDetectorSphereHit hit;
      hit.x0         = CT.x_mid;
      hit.u          = CT.dx_hat;
      hit.rxu        = rxu;
      hit.rxv        = rxv;
      hit.dmin       = dmin;
      hit.cos_phimax = std::min(cos_phimax, -1.0);
      hit.phimax     = std::acos(cos_phimax);
      isphere.processor->process_hit(hit);
      continue;
    }
  }

  visitor_->leave_cherenkov_track();
}

void IACTDetectorSphereAirCherenkovTrackVisitor::leave_event()
{
  visitor_->leave_event();
}
