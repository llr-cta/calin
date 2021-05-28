/*

   calin/simulation/sct_ray_processor.hpp -- Stephen Fegan -- 2021-05-27

   SCT ray tracking ray processor.

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

#pragma once

#include <math/ray.hpp>
#include <math/rng.hpp>
#include <simulation/pe_processor.hpp>
#include <simulation/ray_processor.hpp>
#include <simulation/detector_efficiency.hpp>
#include <simulation/sct_optics.pb.h>
#include <simulation/sct_ray_tracer.hpp>

namespace calin { namespace simulation { namespace sct_ray_processor {

class SCTTracedRayVisitor
{
public:
  virtual ~SCTTracedRayVisitor();
  virtual void start_processing();
  virtual void process_traced_ray(unsigned scope_id,
    const calin::simulation::sct_optics::SCTRayTracerResults& trace, double pe_weight);
  virtual void finish_processing();
};

class SCTMultiTracedRayVisitor: public SCTTracedRayVisitor
{
public:
  SCTMultiTracedRayVisitor();
  virtual ~SCTMultiTracedRayVisitor();
  void start_processing() override;
  void process_traced_ray(unsigned scope_id,
    const calin::simulation::sct_optics::SCTRayTracerResults& trace, double pe_weight) override;
  void finish_processing() override;
  void add_visitor(SCTTracedRayVisitor* visitor, bool adopt_visitor);
private:
  std::vector<SCTTracedRayVisitor*> visitors_;
  std::vector<SCTTracedRayVisitor*> adopted_visitors_;
};

class SCTTracedRayVisitor2PEProcessorAdapter: public SCTTracedRayVisitor
{
public:
  SCTTracedRayVisitor2PEProcessorAdapter(
    calin::simulation::pe_processor::PEProcessor* visitor,
    bool use_fp_position = false,
    calin::simulation::sct_optics::SCTRayTracerStatus status_min
      = calin::simulation::sct_optics::RTS_MISSED_CAMERA,
    bool adopt_visitor = false);
  virtual ~SCTTracedRayVisitor2PEProcessorAdapter();
  void start_processing() override;
  void process_traced_ray(unsigned scope_id,
    const calin::simulation::sct_optics::SCTRayTracerResults& trace, double pe_weight) override;
  void finish_processing() override;
private:
  calin::simulation::pe_processor::PEProcessor* visitor_ = nullptr;
  bool use_fp_position_ = false;
  calin::simulation::sct_optics::SCTRayTracerStatus status_min_;
  bool adopt_visitor_ = false;
};

class SCTRayProcessor: public calin::simulation::ray_processor::RayProcessor
{
public:
  SCTRayProcessor(calin::ix::simulation::sct_optics::SCTArray* array,
    SCTTracedRayVisitor* visitor, calin::math::rng::RNG* rng = nullptr,
    bool adopt_array = false, bool adopt_visitor = false,
    bool adopt_rng = false);
  SCTRayProcessor(calin::ix::simulation::sct_optics::SCTArray* array,
    calin::simulation::pe_processor::PEProcessor* visitor,
    calin::math::rng::RNG* rng = nullptr,
    bool adopt_array = false, bool adopt_visitor = false,
    bool adopt_rng = false);
  virtual ~SCTRayProcessor();

  std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere>
    detector_spheres() override;
  void start_processing() override;
  void process_ray(unsigned scope_id, const calin::math::ray::Ray& ray,
    double pe_weight) override;
  void finish_processing() override;

  void add_fp_hit_trace_visitor(SCTTracedRayVisitor* visitor,
    bool adopt_visitor=false);
  void add_pe_visitor(calin::simulation::pe_processor::PEProcessor* visitor,
    bool adopt_visitor=false);

#if 0
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
#endif

  uint64_t nhit() { return nhit_; }

private:
  calin::ix::simulation::sct_optics::SCTArray* array_ = nullptr;
  bool adopt_array_ = false;
  SCTTracedRayVisitor* visitor_ = nullptr;
  bool adopt_visitor_ = false;
  SCTMultiTracedRayVisitor* multi_visitor_ = nullptr;
  calin::math::rng::RNG* rng_ = nullptr;
  bool adopt_rng_ = false;
  calin::simulation::sct_optics::SCTRayTracer* ray_tracer_ = nullptr;
  uint64_t nhit_ = 0;

#if 0
  struct ScopeResponse {
    bool has_effective_bandwidth = false;
    calin::simulation::detector_efficiency::ACTEffectiveBandwidth effective_bandwidth = { 1.0 };
    double detector_bandwidth = 1.0;
    calin::math::interpolation_1d::InterpLinear1D cone_efficiency = { 1.0 };
  };

  std::vector<ScopeResponse> scope_response_;
  #endif
};

} } } // namespace calin::simulation::sct_ray_processor
