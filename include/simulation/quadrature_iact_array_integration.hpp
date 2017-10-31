/*

   calin/simulation/quadrature_iact_array_integration.hpp
                                              -- Stephen Fegan -- 2016-07-27

   IACT array tracker that does quadrature integration over the cherenkov
   cone.

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

#include <simulation/iact_array_tracker.hpp>
#include <simulation/ray_processor.hpp>
#include <simulation/tracker.pb.h>

namespace calin { namespace simulation { namespace quadrature_iact_array_integration {

class QuadratureIACTArrayPEProcessor
{
public:
  virtual ~QuadratureIACTArrayPEProcessor();
  virtual void visit_event(const calin::simulation::tracker::Event& event,
    bool& kill_event);
  virtual void visit_cherenkov_track(
    const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack& cherenkov_track,
    bool& kill_track);
  virtual void process_pe(unsigned scope_id, unsigned pixel_id,
    double x, double y, double t0, double pe_weight);
  virtual void leave_cherenkov_track();
  virtual void leave_event();
};

class QuadratureIACTArrayIntegrationHitVisitor;

class QuadratureIACTDetectorSphereHitProcessor:
  public calin::simulation::iact_array_tracker::IACTDetectorSphereHitProcessor
{
public:
  QuadratureIACTDetectorSphereHitProcessor(
    QuadratureIACTArrayIntegrationHitVisitor* quadrature, unsigned scope_id);
  virtual ~QuadratureIACTDetectorSphereHitProcessor();
  virtual void process_hit(
    const calin::simulation::iact_array_tracker::IACTDetectorSphereHit& hit);
private:
  QuadratureIACTArrayIntegrationHitVisitor* quadrature_ = nullptr;
  unsigned scope_id_ = 0;
};

class QuadratureIACTArrayIntegrationHitVisitor:
  public calin::simulation::iact_array_tracker::HitIACTVisitor
{
public:
  QuadratureIACTArrayIntegrationHitVisitor(
    const calin::ix::simulation::tracker::QuadratureIACTArrayIntegrationConfig& config,
    calin::simulation::ray_processor::RayProcessor* visitor,
    bool adopt_visitor = false);
  ~QuadratureIACTArrayIntegrationHitVisitor();

  std::vector<calin::simulation::iact_array_tracker::IACTDetectorSphere> spheres() override;
  void visit_event(const calin::simulation::tracker::Event& event,
    bool& kill_event) override;
  void leave_event() override;

protected:
  friend class QuadratureIACTDetectorSphereHitProcessor;
  
  void process_test_ray(
    const calin::simulation::iact_array_tracker::IACTDetectorSphereHit& hit,
    unsigned scope_id, double cos_phi, double sin_phi, double dphi);
  void process_hit_sphere(
    const calin::simulation::iact_array_tracker::IACTDetectorSphereHit& hit,
    unsigned scope_id);

  double ray_spacing_linear_;
  double ray_spacing_angular_;
  calin::simulation::ray_processor::RayProcessor* visitor_ = nullptr;
  bool adopt_visitor_ = false;
};

} } } // namespace calin::simulation::quadrature_iact_array_integration
