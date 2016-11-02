/*

   calin/simulation/vso_quadrature_iact_array_integration.hpp
                                              -- Stephen Fegan -- 2016-07-27

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

#pragma once

#include <math/rng.hpp>
#include <simulation/vso_array.hpp>
#include <simulation/vso_raytracer.hpp>
#include <simulation/detector_efficiency.hpp>
#include <simulation/quadrature_iact_array_integration.hpp>

namespace calin { namespace simulation { namespace quadrature_iact_array_integration {

class VSOQuadratureIACTArrayTraceProcessor
{
public:
  virtual ~VSOQuadratureIACTArrayTraceProcessor();
  virtual void visit_event(const calin::simulation::tracker::Event& event,
    bool& kill_event);
  virtual void visit_cherenkov_track(
    const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack& cherenkov_track,
    bool& kill_track);
  virtual void process_trace(unsigned scope_id,
    const calin::simulation::vs_optics::VSOTraceInfo& trace, double pe_weight);
  virtual void leave_cherenkov_track();
  virtual void leave_event();
};

class VSOTraceMultiProcessor: public VSOQuadratureIACTArrayTraceProcessor
{
public:
  VSOTraceMultiProcessor(): VSOQuadratureIACTArrayTraceProcessor() {
    // nothing to see here
  }
  virtual ~VSOTraceMultiProcessor();
  void add_trace_visitor(VSOQuadratureIACTArrayTraceProcessor* visitor,
    bool adopt_visitor = false);
  void add_pe_visitor(QuadratureIACTArrayPEProcessor* pe_visitor,
    bool adopt_pe_visitor = false);
  void visit_event(const calin::simulation::tracker::Event& event,
    bool& kill_event) override;
  void visit_cherenkov_track(
    const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack& cherenkov_track,
    bool& kill_track) override;
  void process_trace(unsigned scope_id,
    const calin::simulation::vs_optics::VSOTraceInfo& trace, double pe_weight) override;
  void leave_cherenkov_track() override;
  void leave_event() override;
private:
  std::vector<VSOQuadratureIACTArrayTraceProcessor*> visitors_;
  std::vector<VSOQuadratureIACTArrayTraceProcessor*> adopted_visitors_;
};

class VSOTraceToPEAdaptor: public VSOQuadratureIACTArrayTraceProcessor
{
public:
  VSOTraceToPEAdaptor(QuadratureIACTArrayPEProcessor* pe_visitor,
      bool adopt_pe_visitor = false): VSOQuadratureIACTArrayTraceProcessor(),
    pe_visitor_(pe_visitor), adopt_pe_visitor_(adopt_pe_visitor) {
    // nothing to see here
  }
  virtual ~VSOTraceToPEAdaptor();
  void visit_event(const calin::simulation::tracker::Event& event,
    bool& kill_event) override;
  void visit_cherenkov_track(
    const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack& cherenkov_track,
    bool& kill_track) override;
  void process_trace(unsigned scope_id,
    const calin::simulation::vs_optics::VSOTraceInfo& trace, double pe_weight) override;
  void leave_cherenkov_track() override;
  void leave_event() override;
private:
  QuadratureIACTArrayPEProcessor* pe_visitor_ = nullptr;
  bool adopt_pe_visitor_ = false;
};

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

class VSO_QuadratureIACTArrayIntegrationHitVisitor;

class VSO_IACTDetectorSphereHitProcessor:
  public calin::simulation::iact_array_tracker::IACTDetectorSphereHitProcessor
{
public:
  VSO_IACTDetectorSphereHitProcessor(
    VSO_QuadratureIACTArrayIntegrationHitVisitor* quadrature,
    calin::simulation::vs_optics::VSOTelescope* scope);
  virtual ~VSO_IACTDetectorSphereHitProcessor();
  virtual void process_hit(
    const calin::simulation::iact_array_tracker::IACTDetectorSphereHit& hit);
private:
  VSO_QuadratureIACTArrayIntegrationHitVisitor* quadrature_ = nullptr;
  calin::simulation::vs_optics::VSOTelescope* scope_ = nullptr;
};

class VSO_QuadratureIACTArrayIntegrationHitVisitor:
  public calin::simulation::iact_array_tracker::HitIACTVisitor
{
public:
  VSO_QuadratureIACTArrayIntegrationHitVisitor(
    const calin::ix::simulation::tracker::QuadratureIACTArrayIntegrationConfig& config,
    calin::simulation::vs_optics::VSOArray* array,
    VSOQuadratureIACTArrayTraceProcessor* visitor,
    calin::math::rng::RNG* rng = nullptr,
    bool adopt_array = false, bool adopt_visitor = false,
    bool adopt_rng = false);
  VSO_QuadratureIACTArrayIntegrationHitVisitor(
    const calin::ix::simulation::tracker::QuadratureIACTArrayIntegrationConfig& config,
    calin::simulation::vs_optics::VSOArray* array,
    QuadratureIACTArrayPEProcessor* visitor,
    calin::math::rng::RNG* rng = nullptr,
    bool adopt_array = false, bool adopt_visitor = false,
    bool adopt_rng = false);
  ~VSO_QuadratureIACTArrayIntegrationHitVisitor();
  void add_trace_visitor(VSOQuadratureIACTArrayTraceProcessor* visitor,
    bool adopt_visitor = false);
  void add_pe_visitor(QuadratureIACTArrayPEProcessor* visitor,
    bool adopt_visitor = false);
  void add_visitor(VSOQuadratureIACTArrayTraceProcessor* visitor,
    bool adopt_visitor = false) { add_trace_visitor(visitor,adopt_visitor); }
  void add_visitor(QuadratureIACTArrayPEProcessor* visitor,
    bool adopt_visitor = false) { add_pe_visitor(visitor,adopt_visitor); }

  std::vector<calin::simulation::iact_array_tracker::IACTDetectorSphere> spheres() override;
  void visit_event(const calin::simulation::tracker::Event& event,
    bool& kill_event) override;
  void visit_cherenkov_track(
    const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack& cherenkov_track,
    bool& kill_track) override;
  void leave_cherenkov_track() override;
  void leave_event() override;

  void set_detection_efficiencies(
    const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency,
    const calin::simulation::detector_efficiency::AtmosphericAbsorption& atmospheric_absorption,
    double w0,
    const calin::math::interpolation_1d::InterpLinear1D& cone_efficiency = { 1.0 });

protected:
  friend class VSO_IACTDetectorSphereHitProcessor;

  void process_test_ray(
    const calin::simulation::iact_array_tracker::IACTDetectorSphereHit& hit,
    calin::simulation::vs_optics::VSOTelescope* scope,
    double cos_phi, double sin_phi, double dphi);
  void process_hit(
    const calin::simulation::iact_array_tracker::IACTDetectorSphereHit& hit,
    calin::simulation::vs_optics::VSOTelescope* scope);

  double ray_spacing_linear_;
  double ray_spacing_angular_;
  calin::simulation::vs_optics::VSOArray* array_ = nullptr;
  bool adopt_array_ = false;
  VSOQuadratureIACTArrayTraceProcessor* visitor_ = nullptr;
  bool adopt_visitor_ = false;
  VSOTraceMultiProcessor* multi_visitor_ = nullptr;
  calin::math::rng::RNG* rng_ = nullptr;
  bool adopt_rng_ = false;
  calin::simulation::vs_optics::VSORayTracer* ray_tracer_ = nullptr;
  std::vector<calin::simulation::detector_efficiency::ACTEffectiveBandwidth>
    effective_bandwidth_;
  calin::math::interpolation_1d::InterpLinear1D cone_efficiency_;
};

} } } // namespace calin::simulation::quadrature_iact_array_integration
