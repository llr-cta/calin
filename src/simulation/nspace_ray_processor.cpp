/*

   calin/simulation/nspace_ray_processor.cpp -- Stephen Fegan -- 2021-01-27

   Process rays into an NSpace at a given observation level.

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

#include <util/log.hpp>
#include <simulation/nspace_ray_processor.hpp>

using namespace calin::util::log;
using namespace calin::simulation::nspace_ray_processor;

NSpaceRayProcessor::
NSpaceRayProcessor(const calin::ix::simulation::ray_processor::NSpaceRayProcessorConfig& config):
  RayProcessor(),
  config_(config), space_(nspace_axes()), p_(space_.naxes()),
  x0_(config_.x_origin(), config_.y_origin(), config_.observation_altitude()),
{
  // nothing to see here
}

NSpaceRayProcessor::~NSpaceRayProcessor()
{
  // nothing to see here
}

std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere>
NSpaceRayProcessor::detector_spheres()
{
  // struct RayProcessorDetectorSphere
  // {
  //   RayProcessorDetectorSphere() { /* nothing to see here */ }
  //   RayProcessorDetectorSphere(const Eigen::Vector3d& r0_, double radius, unsigned iobs_ = 0):
  //       r0(r0_), radius(radius), iobs(iobs_) { /* nothing to see here */ }
  //   Eigen::Vector3d r0;        // Center of detector sphere [cm]
  //   double radius = 0;         // Squared radius of sphere  [cm^2]
  //   unsigned iobs = 0;         // Observation layer associated with this detector
  // };

  std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere> sphere;
  sphere.r0 = x0_;
  sphere.radius = M_SQRT2 * config_.xy_radius();
  sphere.iobs = config_.observation_level();
  return sphere;
}

void NSpaceRayProcessor::start_processing()
{
  if(config.clear_at_new_event()) {
    nspace_.clear();
  }
}

void NSpaceRayProcessor::process_ray(unsigned scope_id, const calin::math::ray::Ray& ray,
  double pe_weight)
{
  Vector3d x = ray.position() - x0_;
  Vector3d u = rot_ * ray.direction();
  switch(config.axis_variables()) {
  case calin::ix::simulation::ray_processor::XY:
  default:
    p_[0] = x[0];
    p_[1] = x[1];
    break;
  case calin::ix::simulation::ray_processor::UXUY:
    p_[0] = u[0];
    p_[1] = u[2];
    break;
  case calin::ix::simulation::ray_processor::T:
    p_[0] = ray.t();
    break;
  case calin::ix::simulation::ray_processor::XY_UXUY:
    p_[0] = x[0];
    p_[1] = x[1];
    p_[2] = u[0];
    p_[3] = u[2];
    break;
  case calin::ix::simulation::ray_processor::XY_T:
    p_[0] = x[0];
    p_[1] = x[1];
    p_[2] = ray.t();
    break;
  case calin::ix::simulation::ray_processor::XY_UXUY_T:
    p_[0] = x[0];
    p_[1] = x[1];
    p_[2] = u[0];
    p_[3] = u[2];
    p_[4] = ray.t();
    break;
  };
  space_.insert(p_, pe_weight);
}

void NSpaceRayProcessor::finish_processing()
{
  // nothing to see here
}
