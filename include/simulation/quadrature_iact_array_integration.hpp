/*

   calin/simulation/quadrature_iact_array_integration.hpp
                                              -- Stephen Fegan -- 2016-07-27

   IACT array tracker that does quadrature integration over the cherenkov
   cone.

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

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

class QuadratureIACTArrayIntegration;

class QuadratureIACTDetectorSphereCherenkovConeIntersectionProcessor:
  public calin::simulation::iact_array_tracker::IACTDetectorSphereCherenkovConeIntersectionProcessor
{
public:
  QuadratureIACTDetectorSphereCherenkovConeIntersectionProcessor(
    QuadratureIACTArrayIntegration* quadrature, unsigned scope_id);
  virtual ~QuadratureIACTDetectorSphereCherenkovConeIntersectionProcessor();
  virtual void process_hit(
    const calin::simulation::iact_array_tracker::IACTDetectorSphereCherenkovConeIntersection& hit);
private:
  QuadratureIACTArrayIntegration* quadrature_ = nullptr;
  unsigned scope_id_ = 0;
};

class QuadratureIACTArrayIntegration:
  public calin::simulation::iact_array_tracker::IACTDetectorSpherePotentialCherenkovConeIntersectionVisitor
{
public:
  QuadratureIACTArrayIntegration(
    const calin::ix::simulation::tracker::QuadratureIACTArrayIntegrationConfig& config,
    calin::simulation::ray_processor::RayProcessor* visitor,
    bool adopt_visitor = false);
  ~QuadratureIACTArrayIntegration();

  std::vector<calin::simulation::iact_array_tracker::IACTDetectorSphere> spheres() override;
  void visit_event(const calin::simulation::tracker::Event& event,
    bool& kill_event) override;
  void leave_event() override;

protected:
  friend class QuadratureIACTDetectorSphereCherenkovConeIntersectionProcessor;

  void process_test_ray(
    const calin::simulation::iact_array_tracker::IACTDetectorSphereCherenkovConeIntersection& hit,
    unsigned scope_id, double cos_phi, double sin_phi, double dphi);
  void process_hit_sphere(
    const calin::simulation::iact_array_tracker::IACTDetectorSphereCherenkovConeIntersection& hit,
    unsigned scope_id);

  double ray_spacing_linear_;
  double ray_spacing_angular_;
  calin::simulation::ray_processor::RayProcessor* visitor_ = nullptr;
  bool adopt_visitor_ = false;
};

} } } // namespace calin::simulation::quadrature_iact_array_integration
