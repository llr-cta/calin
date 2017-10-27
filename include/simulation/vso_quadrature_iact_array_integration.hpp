/*

   calin/simulation/vso_quadrature_iact_array_integration.hpp
                                              -- Stephen Fegan -- 2016-07-27

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

#pragma once

#include <simulation/vso_ray_processor.hpp>
#include <simulation/quadrature_iact_array_integration.hpp>

namespace calin { namespace simulation { namespace vso_quadrature_iact_array_integration {

#if 0
class VSOTraceStatusProcessor: public VSOQuadratureIACTArrayTraceProcessor
{
public:
  VSOTraceStatusProcessor(bool reset_on_each_event = false):
    VSOQuadratureIACTArrayTraceProcessor(),
    reset_on_each_event_(reset_on_each_event) {
    // nothing to see here
  }
  virtual ~VSOTraceStatusProcessor();
  void visit_event(const calin::simulation::tracker::Event& event,
    bool& kill_event) override;
  void process_trace(unsigned scope_id,
    const calin::simulation::vs_optics::VSOTraceInfo& trace, double pe_weight) override;
  const std::vector<unsigned> raytracer_status(unsigned iscope) const;
private:
  bool reset_on_each_event_ = false;
  std::vector<std::vector<unsigned>> ray_tracer_status_;
};

class VSOScopeDiagnosticTraceProcessor:
  public VSOQuadratureIACTArrayTraceProcessor
{
public:
  VSOScopeDiagnosticTraceProcessor(unsigned iscope = 0):
      VSOQuadratureIACTArrayTraceProcessor(), iscope_(iscope) {
    // nothing to see here
  }
  virtual ~VSOScopeDiagnosticTraceProcessor();
  void visit_event(const calin::simulation::tracker::Event& event,
    bool& kill_event) override;
  void process_trace(unsigned scope_id,
    const calin::simulation::vs_optics::VSOTraceInfo& trace,
    double pe_weight) override;
  const std::vector<double>& reflec_x() const { return reflec_x_; }
  const std::vector<double>& reflec_z() const { return reflec_z_; }
  const std::vector<double>& reflec_w() const { return reflec_w_; }
private:
  unsigned iscope_;
  std::vector<double> reflec_x_;
  std::vector<double> reflec_z_;
  std::vector<double> reflec_w_;
};
#endif

class VSO_QuadratureIACTArrayIntegrationHitVisitor:
  public calin::simulation::quadrature_iact_array_integration::QuadratureIACTArrayIntegrationHitVisitor
{
public:
  VSO_QuadratureIACTArrayIntegrationHitVisitor(
    const calin::ix::simulation::tracker::QuadratureIACTArrayIntegrationConfig& config,
    calin::simulation::vs_optics::VSOArray* array,
    calin::simulation::vso_ray_processor::VSOTracedRayVisitor* visitor,
    calin::math::rng::RNG* rng = nullptr,
    bool adopt_array = false, bool adopt_visitor = false,
    bool adopt_rng = false);
  VSO_QuadratureIACTArrayIntegrationHitVisitor(
    const calin::ix::simulation::tracker::QuadratureIACTArrayIntegrationConfig& config,
    calin::simulation::vs_optics::VSOArray* array,
    calin::simulation::pe_processor::PEProcessor* visitor,
    calin::math::rng::RNG* rng = nullptr,
    bool adopt_array = false, bool adopt_visitor = false,
    bool adopt_rng = false);
  ~VSO_QuadratureIACTArrayIntegrationHitVisitor();

  void add_fp_hit_trace_visitor(
    calin::simulation::vso_ray_processor::VSOTracedRayVisitor* visitor,
    bool adopt_visitor = false);
  void add_pe_visitor(calin::simulation::pe_processor::PEProcessor* visitor,
    bool adopt_visitor = false);
  void add_visitor(
    calin::simulation::vso_ray_processor::VSOTracedRayVisitor* visitor,
    bool adopt_visitor = false) { add_fp_hit_trace_visitor(visitor,adopt_visitor); }
  void add_visitor(calin::simulation::pe_processor::PEProcessor* visitor,
    bool adopt_visitor = false) { add_pe_visitor(visitor,adopt_visitor); }

  void set_all_detection_efficiencies(
    const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency,
    const calin::simulation::detector_efficiency::AtmosphericAbsorption& atmospheric_absorption,
    double w0,
    const calin::simulation::detector_efficiency::AngularEfficiency& cone_efficiency = { 1.0 });
  void set_all_detector_response_without_atmospheric_absorption(
    const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency);
  void set_all_detector_and_atmosphere_response(
    const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency,
    const calin::simulation::detector_efficiency::AtmosphericAbsorption& atmospheric_absorption,
    double w0);
  void set_all_cone_angular_response(
    const calin::simulation::detector_efficiency::AngularEfficiency& cone_efficiency);

  void set_scope_detection_efficiencies(unsigned iscope,
    const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency,
    const calin::simulation::detector_efficiency::AtmosphericAbsorption& atmospheric_absorption,
    double w0,
    const calin::simulation::detector_efficiency::AngularEfficiency& cone_efficiency = { 1.0 });
  void set_scope_detector_response_without_atmospheric_absorption(unsigned iscope,
    const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency);
  void set_scope_detector_and_atmosphere_response(unsigned iscope,
    const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency,
    const calin::simulation::detector_efficiency::AtmosphericAbsorption& atmospheric_absorption,
    double w0);
  void set_scope_cone_angular_response(unsigned iscope,
    const calin::simulation::detector_efficiency::AngularEfficiency& cone_efficiency);
};

} } } // namespace calin::simulation::quadrature_iact_array_integration
