/*

   calin/simulation/vso_ray_processor.cpp -- Stephen Fegan -- 2017-01-16

   VSO ray (weight, scope, ray) processor.

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <math/special.hpp>
#include <simulation/vso_ray_processor.hpp>

using namespace calin::simulation::vso_ray_processor;
using calin::math::special::SQR;

VSOFPHitTraceVisitor::~VSOFPHitTraceVisitor()
{
  // nothing to see here
}

void VSOFPHitTraceVisitor::start_processing()
{
  // nothing to see here
}

void VSOFPHitTraceVisitor::process_fp_hit_trace(unsigned scope_id,
  const calin::simulation::vs_optics::VSOTraceInfo& trace, double pe_weight)
{
  // nothing to see here
}

void VSOFPHitTraceVisitor::finish_processing()
{
  // nothing to see here
}

VSOFPHitTraceVisitor2PEProcessorAdapter::
VSOFPHitTraceVisitor2PEProcessorAdapter(
    calin::simulation::pe_processor::PEProcessor* visitor,
    bool process_pes_without_pixel, bool adopt_visitor):
  VSOFPHitTraceVisitor(),
  visitor_(visitor), adopt_visitor_(adopt_visitor),
  process_pes_without_pixel_(process_pes_without_pixel)
{
  // nothing to see here
}

VSOFPHitTraceVisitor2PEProcessorAdapter::
~VSOFPHitTraceVisitor2PEProcessorAdapter()
{
  if(adopt_visitor_)delete visitor_;
}

void VSOFPHitTraceVisitor2PEProcessorAdapter::start_processing()
{
  visitor_->start_processing();
}

void VSOFPHitTraceVisitor2PEProcessorAdapter::
process_fp_hit_trace(unsigned scope_id,
  const calin::simulation::vs_optics::VSOTraceInfo& trace, double pe_weight)
{
  if(trace.pixel != nullptr) {
    visitor_->process_pe(scope_id, trace.pixel->id(),
      trace.fplane_x, trace.fplane_z, trace.fplane_t, pe_weight);
  } else if(process_pes_without_pixel_) {
    visitor_->process_pe(scope_id, -1,
      trace.fplane_x, trace.fplane_z, trace.fplane_t, pe_weight);
  }
}

void VSOFPHitTraceVisitor2PEProcessorAdapter::finish_processing()
{
  visitor_->finish_processing();
}

VSOMultiFPHitTraceVisitor::VSOMultiFPHitTraceVisitor():
  VSOFPHitTraceVisitor(), visitors_()
{
  // nothing to see here
}

VSOMultiFPHitTraceVisitor::~VSOMultiFPHitTraceVisitor()
{
  for(auto* ivisitor : adopted_visitors_)delete ivisitor;
}

void VSOMultiFPHitTraceVisitor::start_processing()
{
  for(auto* ivisitor : visitors_)ivisitor->start_processing();
}

void VSOMultiFPHitTraceVisitor::process_fp_hit_trace(unsigned scope_id,
  const calin::simulation::vs_optics::VSOTraceInfo& trace, double pe_weight)
{
  for(auto* ivisitor : visitors_)
    ivisitor->process_fp_hit_trace(scope_id, trace, pe_weight);
}

void VSOMultiFPHitTraceVisitor::finish_processing()
{
  for(auto* ivisitor : visitors_)ivisitor->finish_processing();
}

void VSOMultiFPHitTraceVisitor::add_visitor(VSOFPHitTraceVisitor* visitor,
  bool adopt_visitor)
{
  visitors_.push_back(visitor);
  if(adopt_visitor)adopted_visitors_.push_back(visitor);
}

VSORayProcessor::VSORayProcessor(calin::simulation::vs_optics::VSOArray* array,
    VSOFPHitTraceVisitor* visitor, calin::math::rng::RNG* rng,
    bool adopt_array, bool adopt_visitor, bool adopt_rng):
  calin::simulation::ray_processor::RayProcessor(),
  array_(array), adopt_array_(adopt_array),
  visitor_(visitor), adopt_visitor_(adopt_visitor),
  rng_(rng ? rng : new calin::math::rng::RNG),
  adopt_rng_(rng ? true : adopt_rng),
  ray_tracer_(new calin::simulation::vs_optics::VSORayTracer(array_, rng_))
{
  // nothing to see hee
}

VSORayProcessor::VSORayProcessor(calin::simulation::vs_optics::VSOArray* array,
    calin::simulation::pe_processor::PEProcessor* visitor,
    calin::math::rng::RNG* rng,
    bool adopt_array, bool adopt_visitor, bool adopt_rng):
  calin::simulation::ray_processor::RayProcessor(),
  array_(array), adopt_array_(adopt_array),
  visitor_(new VSOFPHitTraceVisitor2PEProcessorAdapter(visitor, adopt_visitor)),
  adopt_visitor_(true),
  rng_(rng ? rng : new calin::math::rng::RNG),
  adopt_rng_(rng ? true : adopt_rng),
  ray_tracer_(new calin::simulation::vs_optics::VSORayTracer(array_, rng_))
{
  // nothing to see here
}

VSORayProcessor::~VSORayProcessor()
{
  if(adopt_array_)delete array_;
  if(adopt_visitor_)delete visitor_;
  if(adopt_rng_)delete rng_;
}

std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere>
VSORayProcessor::detector_spheres()
{
  std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere> s;
  for(unsigned iscope=0; iscope<array_->numTelescopes(); iscope++)
  {
    auto* scope = array_->telescope(iscope);
    Eigen::Vector3d sphere_center = scope->reflectorIPCenter();
    scope->reflectorToGlobal_pos(sphere_center);
    s.emplace_back(sphere_center, SQR(0.5*scope->reflectorIP()));
  }
  return s;
}

void VSORayProcessor::start_processing()
{
  visitor_->start_processing();
}

void VSORayProcessor::process_ray(unsigned scope_id,
  const calin::math::ray::Ray& ray, double pe_weight)
{
  double z0 = ray.position().z();
  double w = ray.direction().z();
  auto* scope = array_->telescope(scope_id);

  calin::math::ray::Ray ray_copy(ray);
  vs_optics::VSOTraceInfo trace_info;
  ray_tracer_->trace(ray_copy, trace_info, scope);

  if(scope_id < effective_bandwidth_.size()) {
    pe_weight *= effective_bandwidth_[scope_id].bandwidth(z0, std::fabs(w));
  }

  if(trace_info.rayHitFocalPlane()) {
    pe_weight *= cone_efficiency_.y(trace_info.fplane_uy);
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

  visitor_->process_fp_hit_trace(scope_id, trace_info, pe_weight);
}

void VSORayProcessor::finish_processing()
{
  visitor_->finish_processing();
}

void VSORayProcessor::add_fp_hit_trace_visitor(VSOFPHitTraceVisitor* visitor,
  bool adopt_visitor)
{
  if(multi_visitor_ == nullptr) {
    multi_visitor_ = new VSOMultiFPHitTraceVisitor;
    multi_visitor_->add_visitor(visitor_, adopt_visitor_);
    visitor_ = multi_visitor_;
    adopt_visitor_ = true;
  }
  multi_visitor_->add_visitor(visitor, adopt_visitor);
}

void VSORayProcessor::add_pe_visitor(
  calin::simulation::pe_processor::PEProcessor* visitor,
  bool adopt_visitor)
{
  add_fp_hit_trace_visitor(
    new VSOFPHitTraceVisitor2PEProcessorAdapter(visitor, adopt_visitor));
}
