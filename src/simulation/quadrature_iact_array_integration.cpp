/*

   calin/simulation/quadrature_iact_array_integration.cpp
                                        -- Stephen Fegan -- 2016-07-28

   IACT array tracker that does quadrature integration over the cherenkov
   cone.

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
#include <math/special.hpp>
#include <math/vs_particle.hpp>

using namespace calin::simulation::quadrature_iact_array_integration;
using calin::math::special::SQR;
using calin::simulation::iact_array_tracker::IACTDetectorSphere;
using namespace calin::simulation::vs_optics;
//using namespace calin::math::vs_physics;

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


VSO_QuadratureIACTArrayIntegrationHitVisitor::
VSO_QuadratureIACTArrayIntegrationHitVisitor(double step_size,
    calin::simulation::vs_optics::VSOArray* array, bool adopt_array):
  calin::simulation::iact_array_tracker::HitIACTVisitor(),
  step_size_(step_size), array_(array), adopt_array_(adopt_array)
{
  // nothing to see here
}

VSO_QuadratureIACTArrayIntegrationHitVisitor::
~VSO_QuadratureIACTArrayIntegrationHitVisitor()
{
  if(adopt_array_)delete array_;
}

std::vector<IACTDetectorSphere>
VSO_QuadratureIACTArrayIntegrationHitVisitor::spheres()
{
  std::vector<IACTDetectorSphere> s;
  for(unsigned iscope=0; iscope<array_->numTelescopes(); iscope++)
  {
    auto* scope = array_->telescope(iscope);
    s.emplace_back(scope->pos(), scope->curvatureRadius(),
      new VSO_IACTDetectorSphereHitProcessor(this, scope));
  }
  return s;
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
  // nothing to see here
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::visit_cherenkov_track(
  const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack& cherenkov_track,
  bool& kill_track)
{
  // nothing to see here
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::leave_cherenkov_track()
{
  // nothing to see here
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::leave_event()
{
  // nothing to see here
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::process_test_ray(
  const calin::simulation::iact_array_tracker::IACTDetectorSphereHit& hit,
  calin::simulation::vs_optics::VSOTelescope* scope,
  double cos_phi, double sin_phi, double dphi)
{
  Eigen::Vector3d p_hat = hit.u * hit.cherenkov_track->cos_thetac +
    (hit.v * cos_phi + hit.w * sin_phi)*hit.cherenkov_track->sin_thetac;

  math::vs_physics::Particle ray;
  VSOTraceInfo trace_info;
  const VSOPixel* pixel = ray_tracer_->trace(ray, trace_info, scope);


}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::process_hit(
  const calin::simulation::iact_array_tracker::IACTDetectorSphereHit& hit,
  calin::simulation::vs_optics::VSOTelescope* scope)
{
  double dphi = step_size_*hit.cherenkov_track->cos_thetac /
    (2.0*M_PI*hit.rxu*hit.cherenkov_track->sin_thetac);
  unsigned n = std::floor(2.0*hit.phimax / dphi);
  if(n%2 == 0)++n; // always choose odd number of integration Points
  dphi = 2.0*hit.phimax/double(n);
  const double cos_dphi = std::cos(dphi);
  const double sin_dphi = std::sin(dphi);
  double cos_phi = 1.0;
  double sin_phi = 0.0;
  for(unsigned i=0; i<n ; i++)
  {
    if(i==0) {
      process_test_ray(hit, scope, cos_phi, sin_phi, dphi);
    } else {
      process_test_ray(hit, scope, cos_phi, sin_phi, dphi);
      process_test_ray(hit, scope, cos_phi, -sin_phi, dphi);
    }
    const double sin_phi_ = sin_phi * cos_dphi + cos_phi * sin_dphi;
    const double cos_phi_ = cos_phi * cos_dphi - sin_phi * sin_dphi;
    sin_phi = sin_phi_;
    cos_phi = cos_phi_;
  }
}
