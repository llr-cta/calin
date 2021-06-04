/*

   calin/simulation/vso_ray_processor.cpp -- Stephen Fegan -- 2021-05-27

   SCT ray (weight, scope, ray) processor.

   Copyright 2021, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <math/special.hpp>
#include <simulation/sct_ray_processor.hpp>

using namespace calin::simulation::sct_ray_processor;
using calin::math::special::SQR;

#include <util/log.hpp>
using namespace calin::util::log;

SCTTracedRayVisitor::~SCTTracedRayVisitor()
{
  // nothing to see here
}

void SCTTracedRayVisitor::start_processing()
{
  // nothing to see here
}

void SCTTracedRayVisitor::process_traced_ray(unsigned scope_id,
  const calin::simulation::sct_optics::SCTRayTracerResults& trace, double pe_weight)
{
  // nothing to see here
}

void SCTTracedRayVisitor::finish_processing()
{
  // nothing to see here
}

SCTRecordingTracedRayVisitor::~SCTRecordingTracedRayVisitor()
{
  // nothing to see here
}

void SCTRecordingTracedRayVisitor::start_processing()
{
  scope_id_.clear();
  results_.clear();
  pe_weight_.clear();
}

void SCTRecordingTracedRayVisitor::process_traced_ray(unsigned scope_id,
  const calin::simulation::sct_optics::SCTRayTracerResults& trace, double pe_weight)
{
  scope_id_.push_back(scope_id);
  results_.push_back(trace);
  pe_weight_.push_back(pe_weight);
}

void SCTRecordingTracedRayVisitor::finish_processing()
{
  // nothing to see here
}

calin::simulation::sct_optics::SCTRayTracerResults
SCTRecordingTracedRayVisitor::results(unsigned ievent, unsigned& scope_id_out,
  double& pe_weight_out) const
{
  if(ievent >= results_.size()) {
    throw std::out_of_range("SCTRecordingTracedRayVisitor::results: ievent out of range");
  }
  scope_id_out = scope_id_[ievent];
  pe_weight_out = pe_weight_[ievent];
  return results_[ievent];
}

SCTTracedRayVisitor2PEProcessorAdapter::
SCTTracedRayVisitor2PEProcessorAdapter(
    calin::simulation::pe_processor::PEProcessor* visitor, bool use_fp_position,
    calin::simulation::sct_optics::SCTRayTracerStatus status_min,
    bool adopt_visitor):
  SCTTracedRayVisitor(),
  visitor_(visitor), use_fp_position_(use_fp_position), status_min_(status_min),
  adopt_visitor_(adopt_visitor)
{
  // nothing to see here
}

SCTTracedRayVisitor2PEProcessorAdapter::
~SCTTracedRayVisitor2PEProcessorAdapter()
{
  if(adopt_visitor_)delete visitor_;
}

void SCTTracedRayVisitor2PEProcessorAdapter::start_processing()
{
  visitor_->start_processing();
}

void SCTTracedRayVisitor2PEProcessorAdapter::
process_traced_ray(unsigned scope_id,
  const calin::simulation::sct_optics::SCTRayTracerResults& trace, double pe_weight)
{
  double x;
  double z;

  if(trace.status >= status_min_) {
    if(use_fp_position_) {
      x = trace.fp_position.x();
      z = trace.fp_position.z();
    } else {
      x = trace.camera_position.x();
      z = trace.camera_position.z();
    }
    visitor_->process_focal_plane_hit(scope_id, trace.camera_pixel_id,
      x, z, trace.camera_time, pe_weight);
  }
}

void SCTTracedRayVisitor2PEProcessorAdapter::finish_processing()
{
  visitor_->finish_processing();
}

SCTMultiTracedRayVisitor::SCTMultiTracedRayVisitor():
  SCTTracedRayVisitor(), visitors_()
{
  // nothing to see here
}

SCTMultiTracedRayVisitor::~SCTMultiTracedRayVisitor()
{
  for(auto* ivisitor : adopted_visitors_)delete ivisitor;
}

void SCTMultiTracedRayVisitor::start_processing()
{
  for(auto* ivisitor : visitors_)ivisitor->start_processing();
}

void SCTMultiTracedRayVisitor::process_traced_ray(unsigned scope_id,
  const calin::simulation::sct_optics::SCTRayTracerResults& trace, double pe_weight)
{
  for(auto* ivisitor : visitors_)
    ivisitor->process_traced_ray(scope_id, trace, pe_weight);
}

void SCTMultiTracedRayVisitor::finish_processing()
{
  for(auto* ivisitor : visitors_)ivisitor->finish_processing();
}

void SCTMultiTracedRayVisitor::add_visitor(SCTTracedRayVisitor* visitor,
  bool adopt_visitor)
{
  visitors_.push_back(visitor);
  if(adopt_visitor)adopted_visitors_.push_back(visitor);
}

void SCTMultiTracedRayVisitor::add_pe_processor(
  calin::simulation::pe_processor::PEProcessor* pe_processor, bool adopt_pe_processor)
{
  add_visitor(new SCTTracedRayVisitor2PEProcessorAdapter(pe_processor, adopt_pe_processor),
    /* adopt_visitor = */ true);
}

SCTRayProcessor::SCTRayProcessor(calin::ix::simulation::sct_optics::SCTArray* array,
    SCTTracedRayVisitor* visitor, calin::math::rng::RNG* rng,
    bool adopt_array, bool adopt_visitor, bool adopt_rng):
  calin::simulation::ray_processor::RayProcessor(),
  array_(array), adopt_array_(adopt_array),
  visitor_(visitor), adopt_visitor_(adopt_visitor),
  rng_(rng ? rng : new calin::math::rng::RNG(__PRETTY_FUNCTION__)),
  adopt_rng_(rng ? adopt_rng : true),
  ray_tracer_(new calin::simulation::sct_optics::SCTRayTracer(array_, rng_)) //,
  // scope_response_(array->numTelescopes())
{
  // nothing to see hee
}

SCTRayProcessor::SCTRayProcessor(calin::ix::simulation::sct_optics::SCTArray* array,
    calin::simulation::pe_processor::PEProcessor* visitor,
    calin::math::rng::RNG* rng,
    bool adopt_array, bool adopt_visitor, bool adopt_rng):
  calin::simulation::ray_processor::RayProcessor(),
  array_(array), adopt_array_(adopt_array),
  visitor_(new SCTTracedRayVisitor2PEProcessorAdapter(visitor, adopt_visitor)),
  adopt_visitor_(true),
  rng_(rng ? rng : new calin::math::rng::RNG(__PRETTY_FUNCTION__)),
  adopt_rng_(rng ? adopt_rng : true),
  ray_tracer_(new calin::simulation::sct_optics::SCTRayTracer(array_, rng_)) //,
  // scope_response_(array->numTelescopes())
{
  // nothing to see here
}

SCTRayProcessor::~SCTRayProcessor()
{
  if(adopt_array_)delete array_;
  if(adopt_visitor_)delete visitor_;
  if(adopt_rng_)delete rng_;
  delete ray_tracer_;
}

std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere>
SCTRayProcessor::detector_spheres()
{
  std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere> s;
  for(unsigned iscope=0; iscope<ray_tracer_->num_telescopes(); iscope++) {
    s.emplace_back(ray_tracer_->detector_sphere_center(iscope),
      ray_tracer_->detector_sphere_radius(iscope));
  }
  return s;
}

void SCTRayProcessor::start_processing()
{
  nhit_ = 0;
  visitor_->start_processing();
}

void SCTRayProcessor::process_ray(unsigned scope_id,
  const calin::math::ray::Ray& ray, double pe_weight)
{
  double z0 = ray.position().z();
  double w = ray.direction().z();

  calin::math::ray::Ray ray_copy(ray);
  sct_optics::SCTRayTracerResults trace_info;
  trace_info.camera_pixel_id = -1;

  ray_tracer_->trace_ray_in_global_frame(scope_id, ray_copy, trace_info);

#if 0
  const ScopeResponse& scope_response { scope_response_[scope_id] };
  if(scope_response.has_effective_bandwidth)
    pe_weight *= scope_response.effective_bandwidth.bandwidth(z0, std::fabs(w));
  else
    pe_weight *= scope_response.detector_bandwidth;

  if(trace_info.rayHitFocalPlane())
    pe_weight *= scope_response.cone_efficiency.y(trace_info.fplane_uy);
#endif

#if 0
  static unsigned counter = 0;
  if(counter++<10) {
    LOG(INFO) << trace_info.status << ' '
      << trace_info.primary_facet << ' ' << trace_info.primary_position.x() << ' ' << trace_info.primary_position.z() << ' '
      << trace_info.secondary_facet << ' ' << trace_info.secondary_position.x() << ' ' << trace_info.secondary_position.z();
  }
#endif

  if(trace_info.camera_pixel_id >= 0) {
    ++nhit_;
  }
  visitor_->process_traced_ray(scope_id, trace_info, pe_weight);
}

void SCTRayProcessor::finish_processing()
{
  visitor_->finish_processing();
}

void SCTRayProcessor::add_fp_hit_trace_visitor(SCTTracedRayVisitor* visitor,
  bool adopt_visitor)
{
  if(multi_visitor_ == nullptr) {
    multi_visitor_ = new SCTMultiTracedRayVisitor;
    multi_visitor_->add_visitor(visitor_, adopt_visitor_);
    visitor_ = multi_visitor_;
    adopt_visitor_ = true;
  }
  multi_visitor_->add_visitor(visitor, adopt_visitor);
}

void SCTRayProcessor::add_pe_visitor(
  calin::simulation::pe_processor::PEProcessor* visitor, bool adopt_visitor)
{
  add_fp_hit_trace_visitor(new SCTTracedRayVisitor2PEProcessorAdapter(visitor,
    adopt_visitor), true);
}

#if 0
void SCTRayProcessor::set_all_detection_efficiencies(
  const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency,
  const calin::simulation::detector_efficiency::AtmosphericAbsorption& atmospheric_absorption,
  double w0,
  const calin::simulation::detector_efficiency::AngularEfficiency& cone_efficiency)
{
  for(unsigned iscope=0; iscope<scope_response_.size(); iscope++)
    set_scope_detection_efficiencies(iscope, detector_efficiency,
      atmospheric_absorption, w0, cone_efficiency);
}

void SCTRayProcessor::set_all_detector_response_without_atmospheric_absorption(
  const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency)
{
  for(unsigned iscope=0; iscope<scope_response_.size(); iscope++)
    set_scope_detector_response_without_atmospheric_absorption(iscope, detector_efficiency);
}

void SCTRayProcessor::set_all_detector_and_atmosphere_response(
  const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency,
  const calin::simulation::detector_efficiency::AtmosphericAbsorption& atmospheric_absorption,
  double w0)
{
  for(unsigned iscope=0; iscope<scope_response_.size(); iscope++)
    set_scope_detector_and_atmosphere_response(iscope, detector_efficiency,
      atmospheric_absorption, w0);
}

void SCTRayProcessor::set_all_cone_angular_response(
  const calin::simulation::detector_efficiency::AngularEfficiency& cone_efficiency)
{
  for(unsigned iscope=0; iscope<scope_response_.size(); iscope++)
    set_scope_cone_angular_response(iscope, cone_efficiency);
}

void SCTRayProcessor::set_scope_detection_efficiencies(unsigned iscope,
  const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency,
  const calin::simulation::detector_efficiency::AtmosphericAbsorption& atmospheric_absorption,
  double w0,
  const calin::simulation::detector_efficiency::AngularEfficiency& cone_efficiency)
{
  if(iscope >= scope_response_.size())
    throw std::out_of_range("Scope not defined: " + std::to_string(iscope));

  auto* scope = array_->telescope(iscope);
  scope_response_[iscope].has_effective_bandwidth = true;
  scope_response_[iscope].effective_bandwidth =
    atmospheric_absorption.integrateBandwidth(
      scope->position().z(), std::fabs(w0), detector_efficiency);
  scope_response_[iscope].detector_bandwidth = detector_efficiency.integrate();
  scope_response_[iscope].cone_efficiency = cone_efficiency;
}

void SCTRayProcessor::
set_scope_detector_response_without_atmospheric_absorption(unsigned iscope,
  const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency)
{
  if(iscope >= scope_response_.size())
    throw std::out_of_range("Scope not defined: " + std::to_string(iscope));

  scope_response_[iscope].has_effective_bandwidth = false;
  scope_response_[iscope].effective_bandwidth =
    calin::simulation::detector_efficiency::ACTEffectiveBandwidth(1.0);
  scope_response_[iscope].detector_bandwidth = detector_efficiency.integrate();
}

void SCTRayProcessor::
set_scope_detector_and_atmosphere_response(unsigned iscope,
  const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency,
  const calin::simulation::detector_efficiency::AtmosphericAbsorption& atmospheric_absorption,
  double w0)
{
  if(iscope >= scope_response_.size())
    throw std::out_of_range("Scope not defined: " + std::to_string(iscope));

  auto* scope = array_->telescope(iscope);
  scope_response_[iscope].has_effective_bandwidth = true;
  scope_response_[iscope].effective_bandwidth =
    atmospheric_absorption.integrateBandwidth(
      scope->position().z(), std::fabs(w0), detector_efficiency);
  scope_response_[iscope].detector_bandwidth = detector_efficiency.integrate();
}

void SCTRayProcessor::set_scope_cone_angular_response(unsigned iscope,
  const calin::simulation::detector_efficiency::AngularEfficiency& cone_efficiency)
{
  if(iscope >= scope_response_.size())
    throw std::out_of_range("Scope not defined: " + std::to_string(iscope));

  scope_response_[iscope].cone_efficiency = cone_efficiency;
}
#endif
