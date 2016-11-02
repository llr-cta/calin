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

VSOQuadratureIACTArrayTraceProcessor::~VSOQuadratureIACTArrayTraceProcessor()
{
  // nothing to see here
}

void VSOQuadratureIACTArrayTraceProcessor::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
  // nothing to see here
}

void VSOQuadratureIACTArrayTraceProcessor::visit_cherenkov_track(
  const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack& cherenkov_track,
  bool& kill_track)
{
  // nothing to see here
}

void VSOQuadratureIACTArrayTraceProcessor::process_trace(unsigned scope_id,
  const calin::simulation::vs_optics::VSOTraceInfo& trace, double pe_weight)
{
  // nothing to see here
}

void VSOQuadratureIACTArrayTraceProcessor::leave_cherenkov_track()
{
  // nothing to see here
}

void VSOQuadratureIACTArrayTraceProcessor::leave_event()
{
  // nothing to see here
}

VSOTraceMultiProcessor::~VSOTraceMultiProcessor()
{
  for(auto* ivisitor : adopted_visitors_)delete ivisitor;
}

void VSOTraceMultiProcessor::
add_trace_visitor(VSOQuadratureIACTArrayTraceProcessor* visitor,
  bool adopt_visitor)
{
  visitors_.emplace_back(visitor);
  if(adopt_visitor)adopted_visitors_.emplace_back(visitor);
}

void VSOTraceMultiProcessor::
add_pe_visitor(QuadratureIACTArrayPEProcessor* pe_visitor,
  bool adopt_pe_visitor)
{
  add_trace_visitor(new VSOTraceToPEAdaptor(pe_visitor, adopt_pe_visitor), true);
}

void VSOTraceMultiProcessor::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
  for(auto* ivisitor : visitors_) {
    ivisitor->visit_event(event, kill_event);
    if(kill_event)break;
  }
}

void VSOTraceMultiProcessor::visit_cherenkov_track(
  const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack& cherenkov_track,
  bool& kill_track)
{
  for(auto* ivisitor : visitors_) {
    ivisitor->visit_cherenkov_track(cherenkov_track, kill_track);
    if(kill_track)break;
  }
}

void VSOTraceMultiProcessor::process_trace(unsigned scope_id,
  const calin::simulation::vs_optics::VSOTraceInfo& trace, double pe_weight)
{
  for(auto* ivisitor : visitors_)
    ivisitor->process_trace(scope_id, trace, pe_weight);
}

void VSOTraceMultiProcessor::leave_cherenkov_track()
{
  for(auto* ivisitor : visitors_)
    ivisitor->leave_cherenkov_track();
}

void VSOTraceMultiProcessor::leave_event()
{
  for(auto* ivisitor : visitors_)
    ivisitor->leave_event();
}

VSOTraceToPEAdaptor::~VSOTraceToPEAdaptor()
{
  if(adopt_pe_visitor_)delete pe_visitor_;
}

void VSOTraceToPEAdaptor::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
  pe_visitor_->visit_event(event, kill_event);
}

void VSOTraceToPEAdaptor::visit_cherenkov_track(
  const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack& cherenkov_track,
  bool& kill_track)
{
  pe_visitor_->visit_cherenkov_track(cherenkov_track, kill_track);
}

void VSOTraceToPEAdaptor::process_trace(unsigned scope_id,
  const calin::simulation::vs_optics::VSOTraceInfo& trace, double pe_weight)
{
  if(trace.status==TS_PE_GENERATED) {
    assert(trace.pixel != nullptr);
    pe_visitor_->process_pe(scope_id, trace.pixel->id(),
      trace.fplane_x, trace.fplane_z, trace.fplane_t, pe_weight);
  }
}

void VSOTraceToPEAdaptor::leave_cherenkov_track()
{
  pe_visitor_->leave_cherenkov_track();
}

void VSOTraceToPEAdaptor::leave_event()
{
  pe_visitor_->leave_event();
}

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

VSOTraceStatusProcessor::~VSOTraceStatusProcessor()
{
  // nothing to see here
}

void VSOTraceStatusProcessor::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
  if(reset_on_each_event_)
    for(auto& istatus : ray_tracer_status_)
      std::fill(istatus.begin(), istatus.end(), 0);
}

void VSOTraceStatusProcessor::process_trace(unsigned scope_id,
  const calin::simulation::vs_optics::VSOTraceInfo& trace, double pe_weight)
{
  if(scope_id > ray_tracer_status_.size())
    ray_tracer_status_.resize(scope_id+1);
  if(trace.status >= ray_tracer_status_[scope_id].size())
    ray_tracer_status_[scope_id].resize(trace.status+1);
  ray_tracer_status_[scope_id][trace.status]++;
}

const std::vector<unsigned>
VSOTraceStatusProcessor::raytracer_status(unsigned iscope) const
{
  if(iscope >= ray_tracer_status_.size())
    throw std::out_of_range("iscope out of range");
  return ray_tracer_status_[iscope];
}

VSOScopeDiagnosticTraceProcessor::~VSOScopeDiagnosticTraceProcessor()
{
  // nothing to see here
}
void VSOScopeDiagnosticTraceProcessor::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
  reflec_x_.clear();
  reflec_z_.clear();
}

void VSOScopeDiagnosticTraceProcessor::process_trace(unsigned scope_id,
  const calin::simulation::vs_optics::VSOTraceInfo& trace, double pe_weight)
{
  if(scope_id != iscope_)return;
  if(trace.rayWasReflected()) {
    reflec_x_.emplace_back(trace.reflec_x);
    reflec_z_.emplace_back(trace.reflec_z);
    reflec_w_.emplace_back(pe_weight);
  }
}

VSO_QuadratureIACTArrayIntegrationHitVisitor::
VSO_QuadratureIACTArrayIntegrationHitVisitor(
    const calin::ix::simulation::tracker::QuadratureIACTArrayIntegrationConfig& config,
    calin::simulation::vs_optics::VSOArray* array,
    VSOQuadratureIACTArrayTraceProcessor* visitor,
    calin::math::rng::RNG* rng,
    bool adopt_array, bool adopt_visitor, bool adopt_rng):
  calin::simulation::iact_array_tracker::HitIACTVisitor(),
  ray_spacing_linear_(config.ray_spacing_linear()),
  ray_spacing_angular_(config.ray_spacing_angular() / 180.0 * M_PI),
  array_(array), adopt_array_(adopt_array),
  visitor_(visitor), adopt_visitor_(adopt_visitor),
  rng_(rng ? rng : new calin::math::rng::RNG),
  adopt_rng_(rng ? true : adopt_rng),
  ray_tracer_(new calin::simulation::vs_optics::VSORayTracer(array_, rng_))
{
  // nothing to see here
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
  rng_(rng ? rng : new calin::math::rng::RNG),
  adopt_rng_(rng ? true : adopt_rng),
  ray_tracer_(new calin::simulation::vs_optics::VSORayTracer(array_, rng_))
{
  add_pe_visitor(visitor);
}

VSO_QuadratureIACTArrayIntegrationHitVisitor::
~VSO_QuadratureIACTArrayIntegrationHitVisitor()
{
  if(adopt_array_)delete array_;
  if(adopt_visitor_)delete visitor_;
  delete ray_tracer_;
  if(adopt_rng_)delete rng_;
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::
add_trace_visitor(VSOQuadratureIACTArrayTraceProcessor* visitor, bool adopt_visitor)
{
  if(visitor_ == nullptr) {
    visitor_ = visitor;
    adopt_visitor_ = adopt_visitor;
  } else {
    if(multi_visitor_ == nullptr) {
      multi_visitor_ = new VSOTraceMultiProcessor;
      multi_visitor_->add_trace_visitor(visitor_, adopt_visitor_);
      visitor_ = multi_visitor_;
      adopt_visitor_ = true;
    }
    multi_visitor_->add_trace_visitor(visitor, adopt_visitor);
  }
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::
add_pe_visitor(QuadratureIACTArrayPEProcessor* visitor, bool adopt_visitor)
{
  add_trace_visitor(new VSOTraceToPEAdaptor(visitor, adopt_visitor), true);
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

void VSO_QuadratureIACTArrayIntegrationHitVisitor::set_detection_efficiencies(
  const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency,
  const calin::simulation::detector_efficiency::AtmosphericAbsorption& atmospheric_absorption,
  double w0,
  const calin::math::interpolation_1d::InterpLinear1D& cone_efficiency)
{
  effective_bandwidth_.clear();
  for(unsigned iscope=0; iscope<array_->numTelescopes(); iscope++) {
    auto* scope = array_->telescope(iscope);
    effective_bandwidth_.emplace_back(atmospheric_absorption.integrateBandwidth(
      scope->position().z(), std::fabs(w0), detector_efficiency));
  }
  cone_efficiency_ = cone_efficiency;
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
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

void VSO_QuadratureIACTArrayIntegrationHitVisitor::process_test_ray(
  const calin::simulation::iact_array_tracker::IACTDetectorSphereHit& hit,
  calin::simulation::vs_optics::VSOTelescope* scope,
  double cos_phi, double sin_phi, double weight)
{
  const unsigned iscope = scope->id();

  Eigen::Vector3d p_hat = hit.rot *
    Eigen::Vector3d(hit.cherenkov_track->cos_thetac,
      hit.cherenkov_track->sin_thetac * cos_phi,
      hit.cherenkov_track->sin_thetac * sin_phi);

  calin::math::ray::Ray ray { hit.cherenkov_track->x_mid, p_hat,
    hit.cherenkov_track->t_mid };

  double z0 = ray.position().z();
  double w = p_hat.z();

  VSOTraceInfo trace_info;
  ray_tracer_->trace(ray, trace_info, scope);

  if(iscope < effective_bandwidth_.size()) {
    weight *= effective_bandwidth_[iscope].bandwidth(z0, std::fabs(w));
    if(trace_info.rayHitFocalPlane())
      weight *= cone_efficiency_.y(trace_info.fplane_uy);
  }

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

  visitor_->process_trace(iscope, trace_info, weight);
}

// This function should eventually be moved to a base class
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
