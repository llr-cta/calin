/*

   calin/simulation/vso_quadrature_iact_array_integration.cpp
                                        -- Stephen Fegan -- 2016-07-28

   IACT array tracker that does quadrature integration over the cherenkov
   cone - instance that uses VSOptics ray tracer

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
#include <simulation/quadrature_iact_array_integration.hpp>
#include <simulation/vso_quadrature_iact_array_integration.hpp>
#include <math/special.hpp>
#include <math/ray.hpp>
#include <io/log.hpp>

using namespace calin::simulation::quadrature_iact_array_integration;
using namespace calin::io::log;
using calin::math::special::SQR;
using calin::simulation::iact_array_tracker::IACTDetectorSphere;
using namespace calin::simulation::vs_optics;

#define ENABLE_RAYTRACER_STATUS

VSO_IACTDetectorSphereHitProcessor::VSO_IACTDetectorSphereHitProcessor(
    VSO_QuadratureIACTArrayIntegrationHitVisitor* quadrature,
    calin::simulation::vs_optics::VSOTelescope* scope):
  IACTDetectorSphereHitProcessor(), quadrature_(quadrature), scope_(scope)
{
  // nothing to see here
}

VSO_IACTDetectorSphereHitProcessor::~VSO_IACTDetectorSphereHitProcessor()
{
  // nothing to see here
}

void VSO_IACTDetectorSphereHitProcessor::
process_hit(
  const calin::simulation::iact_array_tracker::IACTDetectorSphereHit& hit)
{
  quadrature_->process_hit(hit, scope_);
}


VSO_QuadratureIACTArrayIntegrationHitVisitor::
VSO_QuadratureIACTArrayIntegrationHitVisitor(
    const calin::ix::simulation::tracker::QuadratureIACTArrayIntegrationConfig& config,
    calin::simulation::vs_optics::VSOArray* array,
    QuadratureIACTArrayPEProcessor* visitor,
    calin::math::rng::RNG* rng,
    bool adopt_array, bool adopt_visitor, bool adopt_rng):
  calin::simulation::iact_array_tracker::HitIACTVisitor(),
  ray_spacing_linear_(config.ray_spacing_linear()),
  ray_spacing_angular_(config.ray_spacing_angular() / 180.0 * M_PI),
  array_(array), adopt_array_(adopt_array),
  visitor_(visitor), adopt_visitor_(adopt_visitor),
  rng_(rng ? rng : new calin::math::rng::RNG),
  adopt_rng_(rng ? true : adopt_rng),
  ray_tracer_(new calin::simulation::vs_optics::VSORayTracer(array_, rng_)),
  num_hit_(array_->numTelescopes()), num_miss_(array_->numTelescopes()),
  ray_tracer_status_(array_->numTelescopes())
{
  // nothing to see here
}

VSO_QuadratureIACTArrayIntegrationHitVisitor::
~VSO_QuadratureIACTArrayIntegrationHitVisitor()
{
  if(adopt_array_)delete array_;
  if(adopt_visitor_)delete visitor_;
  delete ray_tracer_;
  if(adopt_rng_)delete rng_;
}

std::vector<IACTDetectorSphere>
VSO_QuadratureIACTArrayIntegrationHitVisitor::spheres()
{
  std::vector<IACTDetectorSphere> s;
  for(unsigned iscope=0; iscope<array_->numTelescopes(); iscope++)
  {
    auto* scope = array_->telescope(iscope);
    Eigen::Vector3d sphere_center = scope->reflectorIPCenter();
    scope->reflectorToGlobal_pos(sphere_center);
    s.emplace_back(sphere_center, SQR(0.5*scope->reflectorIP()),
      new VSO_IACTDetectorSphereHitProcessor(this, scope));
  }
  return s;
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
  std::fill(num_hit_.begin(), num_hit_.end(), 0);
  std::fill(num_miss_.begin(), num_miss_.end(), 0);
  visitor_->visit_event(event, kill_event);
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::visit_cherenkov_track(
  const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack& cherenkov_track,
  bool& kill_track)
{
  visitor_->visit_cherenkov_track(cherenkov_track, kill_track);
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::leave_cherenkov_track()
{
  visitor_->leave_cherenkov_track();
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::leave_event()
{
  visitor_->leave_event();
}

unsigned
VSO_QuadratureIACTArrayIntegrationHitVisitor::num_hit(unsigned iscope) const
{
  if(iscope >= num_hit_.size())
    throw std::out_of_range(
      "VSO_QuadratureIACTArrayIntegrationHitVisitor::num_hit: iscope out "
      "of range");
  return num_hit_[iscope];
}

unsigned
VSO_QuadratureIACTArrayIntegrationHitVisitor::num_miss(unsigned iscope) const
{
  if(iscope >= num_miss_.size())
    throw std::out_of_range(
      "VSO_QuadratureIACTArrayIntegrationHitVisitor::num_miss: iscope out "
      "of range");
  return num_miss_[iscope];
}

const std::vector<unsigned> VSO_QuadratureIACTArrayIntegrationHitVisitor::
raytracer_status(unsigned iscope) const
{
  if(iscope >= ray_tracer_status_.size())
    throw std::out_of_range("iscope out of range");

#ifdef ENABLE_RAYTRACER_STATUS
  return ray_tracer_status_[iscope];
#else
  throw std::runtime_error("Ray tracer status not enabled at compile time");
#endif
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::process_test_ray(
  const calin::simulation::iact_array_tracker::IACTDetectorSphereHit& hit,
  calin::simulation::vs_optics::VSOTelescope* scope,
  double cos_phi, double sin_phi, double weight)
{
  Eigen::Vector3d p_hat = hit.rot *
    Eigen::Vector3d(hit.cherenkov_track->cos_thetac,
      hit.cherenkov_track->sin_thetac * cos_phi,
      hit.cherenkov_track->sin_thetac * sin_phi);

  calin::math::ray::Ray ray { hit.cherenkov_track->x_mid, p_hat,
    hit.cherenkov_track->t_mid };

  VSOTraceInfo trace_info;
  const VSOPixel* pixel = ray_tracer_->trace(ray, trace_info, scope);

#if 0
  if(sin_phi==0) {
    LOG(INFO) << trace_info.status << ' ' << trace_info.mirror_hexid << ' '
      << ray.Position().r << ' ' << p_hat.norm() << " [ "
      << hit.u.transpose() << " ] [ "
      << hit.v.transpose() << "] [ "
      << hit.w.transpose() << "] " << hit.cherenkov_track->cos_thetac << ' '
      << hit.cherenkov_track->sin_thetac;
  }
#endif

#ifdef ENABLE_RAYTRACER_STATUS
  if(trace_info.status >= ray_tracer_status_[scope->id()].size())
    ray_tracer_status_[scope->id()].resize(trace_info.status+1);
  ray_tracer_status_[scope->id()][trace_info.status]++;
#endif

  if(trace_info.status==TS_PE_GENERATED and pixel!=nullptr) {
    num_hit_[scope->id()]++;
    visitor_->process_pe(scope->id(), pixel->id(),
      ray.position().x(), ray.position().z(), ray.time(), weight);
  } else {
    num_miss_[scope->id()]++;
  }
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::process_hit(
  const calin::simulation::iact_array_tracker::IACTDetectorSphereHit& hit,
  calin::simulation::vs_optics::VSOTelescope* scope)
{
  double dphi = std::numeric_limits<double>::infinity();
  if(ray_spacing_linear_ > 0)
    dphi = ray_spacing_linear_*hit.cherenkov_track->cos_thetac /
      std::abs(2.0*M_PI*hit.rxu*hit.cherenkov_track->sin_thetac);
  if(ray_spacing_angular_ > 0)
    dphi = std::min(dphi, ray_spacing_angular_);
  unsigned n = std::floor(2.0*hit.phimax / dphi);
  if(n%2 == 0)++n; // always choose odd number of integration points
  dphi = 2.0*hit.phimax/double(n);
  const double cos_dphi = std::cos(dphi);
  const double sin_dphi = std::sin(dphi);
  double cos_phi = 1.0;
  double sin_phi = 0.0;
  unsigned n_half = (n+1)/2;
  double weight = dphi/(2.0*M_PI)*hit.cherenkov_track->yield_density;
  // Adapt weight for track length / QE / number of cone reflections etc...

#if 0
  LOG(INFO) << hit.x0.transpose() << ' '
    << "[ " << hit.u.dot(hit.u) << ' ' << hit.u.dot(hit.v) << ' ' << hit.u.dot(hit.w) << " ], "
    << "[ " << hit.v.dot(hit.u) << ' ' << hit.v.dot(hit.v) << ' ' << hit.v.dot(hit.w) << " ], "
    << "[ " << hit.w.dot(hit.u) << ' ' << hit.w.dot(hit.v) << ' ' << hit.w.dot(hit.w) << " ]";
#endif

  for(unsigned i=0; i<n_half ; i++)
  {
    if(i==0) {
      process_test_ray(hit, scope, cos_phi, sin_phi, weight);
    } else {
      process_test_ray(hit, scope, cos_phi, sin_phi, weight);
      process_test_ray(hit, scope, cos_phi, -sin_phi, weight);
    }
    const double sin_phi_ = sin_phi * cos_dphi + cos_phi * sin_dphi;
    const double cos_phi_ = cos_phi * cos_dphi - sin_phi * sin_dphi;
    sin_phi = sin_phi_;
    cos_phi = cos_phi_;
  }
}
