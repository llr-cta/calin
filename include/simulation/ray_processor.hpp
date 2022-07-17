/*

   calin/simulation/ray_processor.hpp -- Stephen Fegan -- 2017-01-16

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

#include <vector>
#include <math/ray.hpp>
#include <math/ray_generator.hpp>
#include <simulation/pe_processor.hpp>

namespace calin { namespace simulation { namespace ray_processor {

struct RayProcessorDetectorSphere
{
  RayProcessorDetectorSphere() { /* nothing to see here */ }
  RayProcessorDetectorSphere(const Eigen::Vector3d& r0_, double radius, unsigned iobs_ = 0):
    r0(r0_), radius(radius), iobs(iobs_) { /* nothing to see here */ }
  RayProcessorDetectorSphere(const Eigen::Vector3d& r0_, double radius,
      const Eigen::Vector3d& obs_dir_, double field_of_view_radius_, unsigned iobs_ = 0):
    r0(r0_), radius(radius), iobs(iobs_),
    obs_dir(obs_dir_), field_of_view_radius(field_of_view_radius_) { /* nothing to see here */ }

  // Position of sphere around detector through which incoming ray must mast
  // if it is to have any chance of being detected
  Eigen::Vector3d r0;        // Center of detector sphere [cm]
  double radius = 0;         // Radius of sphere  [cm^2]
  unsigned iobs = 0;         // Observation layer associated with this detector

  // Optional restriction on direction that incoming ray must have if it is to
  // be detected. If incoming ray has direction u, then ray must have
  // u.obs_dir <= -cos(field_of_view_radius) for it to be detected
  Eigen::Vector3d obs_dir = Eigen::Vector3d::UnitY(); // Pointing direction of detector
  double field_of_view_radius = M_PI; // Field of view of detector [radians]
};

class RayProcessor
{
public:
  virtual ~RayProcessor();
  virtual std::vector<RayProcessorDetectorSphere> detector_spheres();
  virtual void start_processing();
  virtual void process_ray(unsigned scope_id, const calin::math::ray::Ray& ray,
    double pe_weight);
  virtual void finish_processing();
  uint64_t process_all_from_ray_generator(
    calin::math::ray_generator::RayGenerator* gen, unsigned scope_id = 0);
};

class FederatedPEProcessor: public calin::simulation::pe_processor::PEProcessor
{
public:
  FederatedPEProcessor(calin::simulation::pe_processor::PEProcessor* pe_processor,
    unsigned scope_id_base);
  virtual ~FederatedPEProcessor();
  void start_processing() override;
  void process_focal_plane_hit(unsigned scope_id, int pixel_id,
    double x, double y, double ux, double uy, double t0, double pe_weight) override;
  void finish_processing() override;
protected:
  unsigned scope_id_base_ = 0;
  calin::simulation::pe_processor::PEProcessor* pe_processor_ = nullptr;
};

class FederatedRayProcessor: public RayProcessor
{
public:
  FederatedRayProcessor();
  virtual ~FederatedRayProcessor();
  std::vector<RayProcessorDetectorSphere> detector_spheres() override;
  void start_processing() override;
  void process_ray(unsigned scope_id, const calin::math::ray::Ray& ray,
    double pe_weight) override;
  void finish_processing() override;
  void add_processor(RayProcessor* processor, bool adopt_processor);
  FederatedPEProcessor* add_processor_and_pe_visitor(
    RayProcessor* ray_processor,
    calin::simulation::pe_processor::PEProcessor* pe_processor,
    bool adopt_ray_processor);
protected:
  std::vector<RayProcessor*> processors_;
  std::vector<unsigned> nsphere_;
  std::vector<RayProcessor*> adopted_processors_;
};

} } } // namespace calin::simulation::ray_processor
