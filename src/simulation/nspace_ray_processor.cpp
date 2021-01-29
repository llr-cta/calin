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
using namespace calin::simulation::ray_processor;

NSpaceRayProcessor::
NSpaceRayProcessor(const calin::ix::simulation::ray_processor::NSpaceRayProcessorConfig& config):
  RayProcessor(),
  config_(config), space_(nspace_axes()), p_(space_.naxes()),
  x0_(config_.x_origin(), config_.y_origin(), config_.observation_altitude())
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

  std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere> spheres;
  calin::simulation::ray_processor::RayProcessorDetectorSphere sphere;
  sphere.r0 = x0_;
  sphere.radius = M_SQRT2 * config_.xy_radius();
  sphere.iobs = config_.observation_level();
  spheres.push_back(sphere);
  return spheres;
}

void NSpaceRayProcessor::start_processing()
{
  if(config_.clear_at_new_event()) {
    space_.clear();
  }
}

void NSpaceRayProcessor::process_ray(unsigned scope_id, const calin::math::ray::Ray& ray,
  double pe_weight)
{
  calin::math::ray::Ray ray_copy = ray;
  if(not ray_copy.propagate_to_plane(Eigen::Vector3d::UnitZ(), -x0_[2], /*time_reversal_ok=*/ false)) {
    return;
  }
  Eigen::Vector3d x = ray_copy.position() - x0_;
  Eigen::Vector3d u = rot_ * ray_copy.direction();
  switch(config_.axis_variables()) {
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
    p_[0] = ray.time()*1e9;
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
    p_[2] = ray.time()*1e9;
    break;
  case calin::ix::simulation::ray_processor::XY_UXUY_T:
    p_[0] = x[0];
    p_[1] = x[1];
    p_[2] = u[0];
    p_[3] = u[2];
    p_[4] = ray.time()*1e9;
    break;
  };
  space_.accumulate(p_, pe_weight);
}

void NSpaceRayProcessor::finish_processing()
{
  // nothing to see here
}

std::vector<calin::math::nspace::Axis>
NSpaceRayProcessor::nspace_axes() const
{
  std::vector<calin::math::nspace::Axis> axes;

  switch(config_.axis_variables()) {
  case calin::ix::simulation::ray_processor::XY:
  case calin::ix::simulation::ray_processor::XY_UXUY:
  case calin::ix::simulation::ray_processor::XY_T:
  case calin::ix::simulation::ray_processor::XY_UXUY_T:
  default:
    axes.push_back({-config_.xy_radius(), config_.xy_radius(), config_.xy_num_bins()});
    axes.push_back({-config_.xy_radius(), config_.xy_radius(), config_.xy_num_bins()});
    break;
  case calin::ix::simulation::ray_processor::UXUY:
  case calin::ix::simulation::ray_processor::T:
    // do nothing
    break;
  };

  switch(config_.axis_variables()) {
  case calin::ix::simulation::ray_processor::UXUY:
  case calin::ix::simulation::ray_processor::XY_UXUY:
  case calin::ix::simulation::ray_processor::XY_UXUY_T:
    axes.push_back({-config_.uxuy_radius()/180.0*M_PI, config_.uxuy_radius()/180.0*M_PI, config_.uxuy_num_bins()});
    axes.push_back({-config_.uxuy_radius()/180.0*M_PI, config_.uxuy_radius()/180.0*M_PI, config_.uxuy_num_bins()});
    break;
  case calin::ix::simulation::ray_processor::XY:
  case calin::ix::simulation::ray_processor::XY_T:
  case calin::ix::simulation::ray_processor::T:
  default:
    // do nothing
    break;
  };

  switch(config_.axis_variables()) {
  case calin::ix::simulation::ray_processor::T:
  case calin::ix::simulation::ray_processor::XY_T:
  case calin::ix::simulation::ray_processor::XY_UXUY_T:
    axes.push_back({0.0, config_.t_duration(), config_.t_num_bins()});
    break;
  case calin::ix::simulation::ray_processor::XY:
  case calin::ix::simulation::ray_processor::UXUY:
  case calin::ix::simulation::ray_processor::XY_UXUY:
  default:
    // do nothing
    break;
  };

  return axes;
}

calin::ix::simulation::ray_processor::NSpaceRayProcessorConfig
NSpaceRayProcessor::default_config()
{
  calin::ix::simulation::ray_processor::NSpaceRayProcessorConfig config;
  config.set_axis_variables(calin::ix::simulation::ray_processor::XY);
  config.set_xy_radius(51200);
  config.set_xy_num_bins(1024);
  config.set_uxuy_radius(2.56);
  config.set_uxuy_num_bins(512);
  config.set_t_duration(100000);
  config.set_t_num_bins(100000);
  config.set_observation_altitude(0.0);
  config.set_x_origin(0.0);
  config.set_y_origin(0.0);
  config.set_observation_level(0);
  config.set_clear_at_new_event(false);
  return config;
}
