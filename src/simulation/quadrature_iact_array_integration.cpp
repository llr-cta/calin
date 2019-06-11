/*

   calin/simulation/quadrature_iact_array_integration.cpp
                                        -- Stephen Fegan -- 2016-07-28

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

#include <limits>

#include <simulation/iact_array_tracker.hpp>
#include <simulation/quadrature_iact_array_integration.hpp>
#include <math/special.hpp>
#include <util/log.hpp>

using namespace calin::simulation::quadrature_iact_array_integration;
using namespace calin::util::log;
using calin::math::special::SQR;
using calin::simulation::iact_array_tracker::IACTDetectorSphere;

QuadratureIACTDetectorSphereCherenkovConeIntersectionProcessor::
QuadratureIACTDetectorSphereCherenkovConeIntersectionProcessor(
    QuadratureIACTArrayIntegration* quadrature, unsigned scope_id):
  calin::simulation::iact_array_tracker::IACTDetectorSphereCherenkovConeIntersectionProcessor(),
  quadrature_(quadrature), scope_id_(scope_id)
{
  // nothing to see here
}

QuadratureIACTDetectorSphereCherenkovConeIntersectionProcessor::
~QuadratureIACTDetectorSphereCherenkovConeIntersectionProcessor()
{
  // nothing to see here
}

void QuadratureIACTDetectorSphereCherenkovConeIntersectionProcessor::process_hit(
  const calin::simulation::iact_array_tracker::IACTDetectorSphereCherenkovConeIntersection& hit)
{
  quadrature_->process_hit_sphere(hit, scope_id_);
}

QuadratureIACTArrayIntegration::
QuadratureIACTArrayIntegration(
    const calin::ix::simulation::tracker::QuadratureIACTArrayIntegrationConfig& config,
    calin::simulation::ray_processor::RayProcessor* visitor,
    bool adopt_visitor):
  calin::simulation::iact_array_tracker::IACTDetectorSpherePotentialCherenkovConeIntersectionVisitor(),
  ray_spacing_linear_(config.ray_spacing_linear()),
  ray_spacing_angular_(config.ray_spacing_angular() / 180.0 * M_PI),
  visitor_(visitor), adopt_visitor_(adopt_visitor)
{
  // nothing to see here
}

QuadratureIACTArrayIntegration::
~QuadratureIACTArrayIntegration()
{
  if(adopt_visitor_)delete visitor_;
}

std::vector<calin::simulation::iact_array_tracker::IACTDetectorSphere>
QuadratureIACTArrayIntegration::spheres()
{
  auto ray_processor_spheres = visitor_->detector_spheres();
  std::vector<calin::simulation::iact_array_tracker::IACTDetectorSphere> s;
  unsigned scope_id = 0;
  for(auto isphere : ray_processor_spheres) {
    s.emplace_back(isphere.r0, SQR(isphere.radius),
      new QuadratureIACTDetectorSphereCherenkovConeIntersectionProcessor(this, scope_id++));
  }
  return s;
}

void QuadratureIACTArrayIntegration::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
  visitor_->start_processing();
}

void QuadratureIACTArrayIntegration::leave_event()
{
  visitor_->finish_processing();
}

void QuadratureIACTArrayIntegration::process_test_ray(
  const calin::simulation::iact_array_tracker::IACTDetectorSphereCherenkovConeIntersection& hit,
  unsigned scope_id, double cos_phi, double sin_phi, double weight)
{
  Eigen::Vector3d p_hat = hit.rot *
    Eigen::Vector3d(hit.cherenkov_track->cos_thetac,
      hit.cherenkov_track->sin_thetac * cos_phi,
      hit.cherenkov_track->sin_thetac * sin_phi);

  calin::math::ray::Ray ray { hit.cherenkov_track->x_mid, p_hat,
    hit.cherenkov_track->t_mid };

  visitor_->process_ray(scope_id, ray, weight);
}

void QuadratureIACTArrayIntegration::process_hit_sphere(
  const calin::simulation::iact_array_tracker::IACTDetectorSphereCherenkovConeIntersection& hit,
  unsigned scope_id)
{
  double dphi = std::numeric_limits<double>::infinity();
  if(ray_spacing_linear_ > 0)
    dphi = ray_spacing_linear_*hit.cherenkov_track->cos_thetac /
      std::abs(2.0*M_PI*hit.rxu*hit.cherenkov_track->sin_thetac);
  if(ray_spacing_angular_ > 0)
    dphi = std::min(dphi, ray_spacing_angular_);
  unsigned n = std::floor(2.0*hit.phimax / dphi);
  if(n%2 == 0)++n; // always choose odd number of integration points
  dphi = 2.0*hit.phimax/double(n);
  const double cos_dphi = std::cos(dphi);
  const double sin_dphi = std::sin(dphi);
  double cos_phi = 1.0;
  double sin_phi = 0.0;
  unsigned n_half = (n+1)/2;
  double weight = dphi/(2.0*M_PI)*hit.cherenkov_track->yield_density;
  // Adapt weight for track length / QE / number of cone reflections etc...

#if 0
  LOG(INFO) << hit.x0.transpose() << ' '
    << "[ " << hit.u.dot(hit.u) << ' ' << hit.u.dot(hit.v) << ' ' << hit.u.dot(hit.w) << " ], "
    << "[ " << hit.v.dot(hit.u) << ' ' << hit.v.dot(hit.v) << ' ' << hit.v.dot(hit.w) << " ], "
    << "[ " << hit.w.dot(hit.u) << ' ' << hit.w.dot(hit.v) << ' ' << hit.w.dot(hit.w) << " ]";
#endif

  for(unsigned i=0; i<n_half ; i++)
  {
    if(i==0) {
      process_test_ray(hit, scope_id, cos_phi, sin_phi, weight);
    } else {
      process_test_ray(hit, scope_id, cos_phi, sin_phi, weight);
      process_test_ray(hit, scope_id, cos_phi, -sin_phi, weight);
    }
    const double sin_phi_ = sin_phi * cos_dphi + cos_phi * sin_dphi;
    const double cos_phi_ = cos_phi * cos_dphi - sin_phi * sin_dphi;
    sin_phi = sin_phi_;
    cos_phi = cos_phi_;
  }
}
