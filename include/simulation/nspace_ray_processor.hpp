/*

   calin/simulation/nspace_ray_processor.hpp -- Stephen Fegan -- 2021-01-15

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

#pragma once

#include <vector>
#include <math/ray.hpp>
#include <math/nspace.hpp>
#include <simulation/ray_processor.hpp>
#include <simulation/ray_processor.pb.h>

namespace calin { namespace simulation { namespace ray_processor {

class NSpaceRayProcessor: public RayProcessor
{
public:
  NSpaceRayProcessor(const calin::ix::simulation::ray_processor::NSpaceRayProcessorConfig& config);
  virtual ~NSpaceRayProcessor();
  std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere> detector_spheres() override;
  void start_processing() override;
  void process_ray(unsigned scope_id, const calin::math::ray::Ray& ray,
    double pe_weight) override;
  void finish_processing() override;
  static calin::ix::simulation::ray_processor::NSpaceRayProcessorConfig default_config();
  const calin::math::nspace::SparseNSpace& nspace() const { return space_; }
protected:
  std::vector<calin::math::nspace::Axis> nspace_axes() const;
  calin::ix::simulation::ray_processor::NSpaceRayProcessorConfig config_;
  calin::math::nspace::SparseNSpace space_;
  Eigen::VectorXd p_;
  Eigen::Vector3d x0_;
  Eigen::Matrix3d rot_ = Eigen::Matrix3d::Identity();
};

} } } // namespace calin::simulation::ray_processor
