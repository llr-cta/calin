/*

   calin/simulation/iact_array_tracker.cpp -- Stephen Fegan -- 2016-07-24

   Base class for all air shower track visitors

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

#include <limits>
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <util/log.hpp>
#include <simulation/iact_array_tracker.hpp>
#include <math/special.hpp>

using namespace calin::util::log;
using namespace calin::simulation::iact_array_tracker;
using calin::math::special::SQR;

IACTDetectorSphereCherenkovConeIntersectionProcessor::~IACTDetectorSphereCherenkovConeIntersectionProcessor()
{
  // nothing to see here
}

IACTDetectorSpherePotentialCherenkovConeIntersectionVisitor::~IACTDetectorSpherePotentialCherenkovConeIntersectionVisitor()
{
  // nothing to see here
}

void IACTDetectorSpherePotentialCherenkovConeIntersectionVisitor::visit_event(const calin::simulation::tracker::Event& event,
  bool& kill_event)
{
  // nothing to see here
}

void IACTDetectorSpherePotentialCherenkovConeIntersectionVisitor::visit_cherenkov_track(
  const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack& cherenkov_track,
  bool& kill_track)
{
  // nothing to see here
}

void IACTDetectorSpherePotentialCherenkovConeIntersectionVisitor::leave_cherenkov_track()
{
  // nothing to see here
}

void IACTDetectorSpherePotentialCherenkovConeIntersectionVisitor::leave_event()
{
  // nothing to see here
}

IACTDetectorSphereCherenkovConeIntersectionFinder::
IACTDetectorSphereCherenkovConeIntersectionFinder(IACTDetectorSpherePotentialCherenkovConeIntersectionVisitor* visitor,
    bool adopt_visitor):
  calin::simulation::air_cherenkov_tracker::AirCherenkovTrackVisitor(),
  visitor_(visitor), adopt_visitor_(adopt_visitor),
  spheres_(visitor->spheres())
{
  // nothing to see here
#if 0
  for(auto sphere : spheres_)
    LOG(INFO) << "Sphere: [ " << sphere.r0.transpose() << "] "
      << sphere.radius_sq;
#endif
}

IACTDetectorSphereCherenkovConeIntersectionFinder::
~IACTDetectorSphereCherenkovConeIntersectionFinder()
{
  for(auto& isphere : spheres_) {
    delete isphere.cone_processor;
  }
  if(adopt_visitor_)delete visitor_;
}

void IACTDetectorSphereCherenkovConeIntersectionFinder::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
  visitor_->visit_event(event, kill_event);
}

void IACTDetectorSphereCherenkovConeIntersectionFinder::visit_cherenkov_track(
  const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack& CT,
  bool& kill_track)
{
  visitor_->visit_cherenkov_track(CT, kill_track);
  if(kill_track)return;

  for(auto& isphere : spheres_)
  {
    // Author's note : see Calin notebook for sphere / cone intersection

    IACTDetectorSphereCherenkovConeIntersection hit;
    hit.x0                = CT.x_mid;
    hit.rx                = isphere.r0 - CT.x_mid;
    hit.u                 = CT.dx_hat;

    hit.rx2               = hit.rx.squaredNorm();
    hit.rxu               = hit.rx.dot(CT.dx_hat); // note bad naming of dx_hat
    double rxu2 = SQR(hit.rxu);
    hit.rxv               = std::sqrt(std::max(hit.rx2 - rxu2, 0.0));;

    if(hit.rxv < 0.1)
    {
      if(std::abs(hit.u(2)) < 0.1) {
        hit.v             = Eigen::Vector3d(-hit.u(1),hit.u(0),0.0)/
          sqrt(SQR(hit.u(0))+SQR(hit.u(1)));
      } else {
        hit.v             = Eigen::Vector3d(-hit.u(2),0.0,hit.u(0))/
          sqrt(SQR(hit.u(0))+SQR(hit.u(2)));
      }
    }
    else
    {
      hit.v               = (hit.rx - hit.rxu*hit.u)*(1.0/hit.rxv);
    }
    hit.w                 = hit.u.cross(hit.v);

    hit.rot              << hit.u, hit.v, hit.w;


    double rx2_minus_r2 = hit.rx2 - isphere.radius_sq;

    if(rx2_minus_r2 < 0) // Particle is inside detector sphere - always hit!
    {
      hit.dmin            = std::numeric_limits<double>::quiet_NaN();
      hit.cos_phimax      = -1.0;
      hit.phimax          = M_PI;
      hit.cherenkov_track = &CT;
  #if 0
      LOG(INFO) << "1 " << hit.rxv << ' ' << hit.rxu << " [ "
        << hit.u.transpose() << " ] [ "
        << hit.v.transpose() << "] [ "
        << hit.w.transpose() << "] " << rx2_minus_r2
        << '\n';
  #endif
      isphere.cone_processor->process_hit(hit);
      continue;
    }

    hit.dmin = CT.cos_thetac*hit.rxv - CT.sin_thetac*hit.rxu;
    if(SQR(hit.dmin) > isphere.radius_sq) {
#if 0
      LOG(INFO) << "2 " << hit.rxv << ' ' << hit.rxu << " [ "
        << hit.u.transpose() << " ] [ "
        << hit.v.transpose() << "] [ "
        << hit.w.transpose() << "] " << rx2_minus_r2 << ' '
        << hit.dmin << '\n';
#endif
      continue; // Cone misses sphere
    }

    double cos_phimax = (sqrt(rx2_minus_r2) - CT.cos_thetac*hit.rxu) /
      (CT.sin_thetac*hit.rxv);

#if 0
    LOG(INFO) << "3 " << hit.rxv << ' ' << hit.rxu << " [ "
      << hit.u.transpose() << " ] [ "
      << hit.v.transpose() << "] [ "
      << hit.w.transpose() << "] " << rx2_minus_r2 << ' '
      << hit.dmin << ' ' << cos_phimax << '\n';
#endif

    // This guards against cases where the sphere intersects the "reverse" cone
    if(cos_phimax < 1.0)
    {
      if(cos_phimax <= -1.0)
      {
        hit.cos_phimax    = -1.0;
        hit.phimax        = M_PI;
      }
      else
      {
        hit.cos_phimax    = cos_phimax;
        hit.phimax        = std::acos(cos_phimax);
      }
      hit.cherenkov_track = &CT;
      isphere.cone_processor->process_hit(hit);
      continue;
    }
  }

  visitor_->leave_cherenkov_track();
}

void IACTDetectorSphereCherenkovConeIntersectionFinder::leave_event()
{
  visitor_->leave_event();
}

IACTDetectorSphereCherenkovPhotonIntersectionFinder::
IACTDetectorSphereCherenkovPhotonIntersectionFinder(
    calin::simulation::ray_processor::RayProcessor* visitor,
    calin::simulation::atmosphere::LayeredRefractiveAtmosphere* atm,
    bool adopt_visitor, bool adopt_atm):
  CherenkovPhotonVisitor(),
  visitor_(visitor), adopt_visitor_(adopt_visitor),
  atm_(atm), adopt_atm_(adopt_atm)
{
  // nothing to see here
}

IACTDetectorSphereCherenkovPhotonIntersectionFinder::
~IACTDetectorSphereCherenkovPhotonIntersectionFinder()
{
  if(adopt_visitor_)delete visitor_;
  if(adopt_atm_)delete atm_;
}

void IACTDetectorSphereCherenkovPhotonIntersectionFinder::
set_bandpass(double epsilon0, double bandwidth, bool do_color_photons)
{

}

void IACTDetectorSphereCherenkovPhotonIntersectionFinder::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
  visitor_->start_processing();
}

void IACTDetectorSphereCherenkovPhotonIntersectionFinder::
visit_cherenkov_photon(const calin::simulation::air_cherenkov_tracker::CherenkovPhoton& cherenkov_photon)
{
  calin::math::ray::Ray ray(cherenkov_photon.x0, cherenkov_photon.u0,
    cherenkov_photon.t0, cherenkov_photon.epsilon);
  unsigned scope_id = 0;
  for(auto& isphere : visitor_->detector_spheres()) {
    Eigen::Vector3d dx = cherenkov_photon.x0 - isphere.r0;
    double u0_dot_dx = cherenkov_photon.u0.dot(dx);
    double dr2 = dx.squaredNorm() - u0_dot_dx*u0_dot_dx;
    if(dr2 < SQR(isphere.radius + refraction_safety_margin_)) {
      // we have a potential hit !
      if(do_refraction_) {
        // Here we propagate to the observation plane associated with the
        // detector sphere then refract the ray then bac-propagate to the
        // height of the photon emission so absoption can be calculated later
        atm_->propagate_ray_with_refraction(ray, isphere.iobs, /* time_reversal_ok = */ false);
        double reverse_propagate_dist = (cherenkov_photon.x0.z() - ray.z())/ray.uz();
        ray.propagate_dist(reverse_propagate_dist);
      }
      visitor_->process_ray(scope_id, ray, /* pe_weight = */ 1.0);
    }
    ++scope_id;
  }
}

void IACTDetectorSphereCherenkovPhotonIntersectionFinder::leave_event()
{
  visitor_->finish_processing();
}
