/*

   calin/simulation/vso_ray_processor.hpp -- Stephen Fegan -- 2017-01-16

   Multi-purpose ray processor.

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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
#include <simulation/vso_array.hpp>
#include <simulation/vso_raytracer.hpp>

namespace calin { namespace simulation { namespace vso_ray_processor {

class VSOTracedRayVisitor
{
public:
  virtual ~VSOTracedRayVisitor();
  virtual void start_processing();
  virtual void process_traced_ray(unsigned scope_id,
    const calin::simulation::vs_optics::VSOTraceInfo& trace, double pe_weight);
  virtual void finish_processing();
};

class VSOMultiTracedRayVisitor: public VSOTracedRayVisitor
{
public:
  VSOMultiTracedRayVisitor();
  virtual ~VSOMultiTracedRayVisitor();
  void start_processing() override;
  void process_traced_ray(unsigned scope_id,
    const calin::simulation::vs_optics::VSOTraceInfo& trace, double pe_weight) override;
  void finish_processing() override;
  void add_visitor(VSOTracedRayVisitor* visitor, bool adopt_visitor);
private:
  std::vector<VSOTracedRayVisitor*> visitors_;
  std::vector<VSOTracedRayVisitor*> adopted_visitors_;
};

class VSOTracedRayVisitor2PEProcessorAdapter: public VSOTracedRayVisitor
{
public:
  VSOTracedRayVisitor2PEProcessorAdapter(
    calin::simulation::pe_processor::PEProcessor* visitor,
    bool adopt_visitor = false);
  virtual ~VSOTracedRayVisitor2PEProcessorAdapter();
  void start_processing() override;
  void process_traced_ray(unsigned scope_id,
    const calin::simulation::vs_optics::VSOTraceInfo& trace, double pe_weight) override;
  void finish_processing() override;
private:
  calin::simulation::pe_processor::PEProcessor* visitor_ = nullptr;
  bool adopt_visitor_ = false;
};

class VSORayProcessor: public calin::simulation::ray_processor::RayProcessor
{
public:
  VSORayProcessor(calin::simulation::vs_optics::VSOArray* array,
    VSOTracedRayVisitor* visitor, calin::math::rng::RNG* rng = nullptr,
    bool adopt_array = false, bool adopt_visitor = false,
    bool adopt_rng = false);
  VSORayProcessor(calin::simulation::vs_optics::VSOArray* array,
    calin::simulation::pe_processor::PEProcessor* visitor,
    calin::math::rng::RNG* rng = nullptr,
    bool adopt_array = false, bool adopt_visitor = false,
    bool adopt_rng = false);
  virtual ~VSORayProcessor();

  std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere>
    detector_spheres() override;
  void start_processing() override;
  void process_ray(unsigned scope_id, const calin::math::ray::Ray& ray,
    double pe_weight) override;
  void finish_processing() override;

  void add_fp_hit_trace_visitor(VSOTracedRayVisitor* visitor,
    bool adopt_visitor=false);
  void add_pe_visitor(calin::simulation::pe_processor::PEProcessor* visitor,
    bool adopt_visitor=false);

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

  uint64_t nhit() { return nhit_; }

  static std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere>
    detector_spheres_for_array(const calin::simulation::vs_optics::VSOArray* array);

private:
  calin::simulation::vs_optics::VSOArray* array_ = nullptr;
  bool adopt_array_ = false;
  VSOTracedRayVisitor* visitor_ = nullptr;
  bool adopt_visitor_ = false;
  VSOMultiTracedRayVisitor* multi_visitor_ = nullptr;
  calin::math::rng::RNG* rng_ = nullptr;
  bool adopt_rng_ = false;
  calin::simulation::vs_optics::VSORayTracer* ray_tracer_ = nullptr;
  uint64_t nhit_ = 0;

  struct ScopeResponse {
    bool has_effective_bandwidth = false;
    calin::simulation::detector_efficiency::ACTEffectiveBandwidth effective_bandwidth = { 1.0 };
    double detector_bandwidth = 1.0;
    calin::math::interpolation_1d::InterpLinear1D cone_efficiency = { 1.0 };
    calin::simulation::detector_efficiency::DetectionEfficiency detector_efficiency = { 1.0 };
    calin::simulation::detector_efficiency::AtmosphericAbsorption atmospheric_absorption;
  };

  std::vector<ScopeResponse> scope_response_;
};

} } } // namespace calin::simulation::ray_processor
