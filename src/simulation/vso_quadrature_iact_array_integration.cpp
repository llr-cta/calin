/*

   calin/simulation/vso_quadrature_iact_array_integration.cpp
                                        -- Stephen Fegan -- 2016-07-28

   IACT array tracker that does quadrature integration over the cherenkov
   cone - instance that uses VSOptics ray tracer

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

#include <simulation/iact_array_tracker.hpp>
#include <simulation/quadrature_iact_array_integration.hpp>
#include <simulation/vso_quadrature_iact_array_integration.hpp>
#include <math/special.hpp>
#include <math/ray.hpp>
#include <util/log.hpp>

using namespace calin::simulation::vso_quadrature_iact_array_integration;
using namespace calin::simulation::quadrature_iact_array_integration;
using namespace calin::simulation::vso_ray_processor;

#if 0
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
#endif

VSO_QuadratureIACTArrayIntegrationHitVisitor::
VSO_QuadratureIACTArrayIntegrationHitVisitor(
    const calin::ix::simulation::tracker::QuadratureIACTArrayIntegrationConfig& config,
    calin::simulation::vs_optics::VSOArray* array,
    calin::simulation::vso_ray_processor::VSOTracedRayVisitor* visitor,
    calin::math::rng::RNG* rng,
    bool adopt_array, bool adopt_visitor, bool adopt_rng):
  QuadratureIACTArrayIntegrationHitVisitor(config,
    new VSORayProcessor(array, visitor, rng, adopt_array, adopt_visitor, adopt_rng),
    true)
{
  // nothing to see here
}

VSO_QuadratureIACTArrayIntegrationHitVisitor::
VSO_QuadratureIACTArrayIntegrationHitVisitor(
    const calin::ix::simulation::tracker::QuadratureIACTArrayIntegrationConfig& config,
    calin::simulation::vs_optics::VSOArray* array,
    calin::simulation::pe_processor::PEProcessor* visitor,
    calin::math::rng::RNG* rng,
    bool adopt_array, bool adopt_visitor, bool adopt_rng):
  QuadratureIACTArrayIntegrationHitVisitor(config,
    new VSORayProcessor(array, visitor, rng, adopt_array, adopt_visitor, adopt_rng),
    true)
{
  // nothing to see here
}

VSO_QuadratureIACTArrayIntegrationHitVisitor::
~VSO_QuadratureIACTArrayIntegrationHitVisitor()
{
  // nothing to see here
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::add_fp_hit_trace_visitor(
  calin::simulation::vso_ray_processor::VSOTracedRayVisitor* visitor,
  bool adopt_visitor)
{
  dynamic_cast<VSORayProcessor*>(visitor_)
    ->add_fp_hit_trace_visitor(visitor,adopt_visitor);
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::add_pe_visitor(
  calin::simulation::pe_processor::PEProcessor* visitor,
  bool adopt_visitor)
{
  dynamic_cast<VSORayProcessor*>(visitor_)
    ->add_pe_visitor(visitor,adopt_visitor);
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::set_all_detection_efficiencies(
  const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency,
  const calin::simulation::detector_efficiency::AtmosphericAbsorption& atmospheric_absorption,
  double w0,
  const calin::simulation::detector_efficiency::AngularEfficiency& cone_efficiency)
{
  dynamic_cast<VSORayProcessor*>(visitor_)->
    set_all_detection_efficiencies(
      detector_efficiency, atmospheric_absorption, w0, cone_efficiency);
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::set_all_detector_response_without_atmospheric_absorption(
  const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency)
{
  dynamic_cast<VSORayProcessor*>(visitor_)->
    set_all_detector_response_without_atmospheric_absorption(detector_efficiency);
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::set_all_detector_and_atmosphere_response(
  const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency,
  const calin::simulation::detector_efficiency::AtmosphericAbsorption& atmospheric_absorption,
  double w0)
{
  dynamic_cast<VSORayProcessor*>(visitor_)->
    set_all_detector_and_atmosphere_response(detector_efficiency, atmospheric_absorption, w0);
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::set_all_cone_angular_response(
  const calin::simulation::detector_efficiency::AngularEfficiency& cone_efficiency)
{
  dynamic_cast<VSORayProcessor*>(visitor_)->set_all_cone_angular_response(cone_efficiency);
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::
set_scope_detection_efficiencies(unsigned iscope,
  const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency,
  const calin::simulation::detector_efficiency::AtmosphericAbsorption& atmospheric_absorption,
  double w0,
  const calin::simulation::detector_efficiency::AngularEfficiency& cone_efficiency)
{
  dynamic_cast<VSORayProcessor*>(visitor_)->
    set_scope_detection_efficiencies(
      iscope, detector_efficiency, atmospheric_absorption, w0, cone_efficiency);
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::
set_scope_detector_response_without_atmospheric_absorption(unsigned iscope,
  const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency)
{
  dynamic_cast<VSORayProcessor*>(visitor_)->
    set_scope_detector_response_without_atmospheric_absorption(iscope, detector_efficiency);
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::
set_scope_detector_and_atmosphere_response(unsigned iscope,
  const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency,
  const calin::simulation::detector_efficiency::AtmosphericAbsorption& atmospheric_absorption,
  double w0)
{
  dynamic_cast<VSORayProcessor*>(visitor_)->
    set_scope_detector_and_atmosphere_response(iscope, detector_efficiency, atmospheric_absorption, w0);
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::
set_scope_cone_angular_response(unsigned iscope,
  const calin::simulation::detector_efficiency::AngularEfficiency& cone_efficiency)
{
  dynamic_cast<VSORayProcessor*>(visitor_)->
    set_scope_cone_angular_response(iscope, cone_efficiency);
}
