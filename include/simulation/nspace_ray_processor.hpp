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
#include <simulation/detector_efficiency.hpp>

namespace calin { namespace simulation { namespace ray_processor {

class NSpaceRayProcessor: public RayProcessor
{
public:
  NSpaceRayProcessor(const calin::ix::simulation::ray_processor::NSpaceRayProcessorConfig& config,
    unsigned nspace_log2_block_size = 0);
  virtual ~NSpaceRayProcessor();
  std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere> detector_spheres() override;
  void start_processing() override;
  void process_ray(unsigned scope_id, const calin::math::ray::Ray& ray,
    double pe_weight) override;
  void finish_processing() override;
  static calin::ix::simulation::ray_processor::NSpaceRayProcessorConfig default_config();
  void clear();
  unsigned nevent() const { return nevent_; }
  const calin::math::nspace::BlockSparseNSpace& nspace() const { return space_; }
  calin::math::nspace::BlockSparseNSpace& mutable_nspace() { return space_; }

  void set_detection_efficiencies(
    double epsilon0, double bandwidth,
    const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency,
    const calin::simulation::detector_efficiency::AtmosphericAbsorption& atmospheric_absorption,
    double w0);
  calin::simulation::detector_efficiency::ACTEffectiveBandwidth* effective_bandwidth() const { return effective_bandwidth_; }

protected:
  std::vector<calin::math::nspace::Axis> nspace_axes() const;
  calin::ix::simulation::ray_processor::NSpaceRayProcessorConfig config_;
  calin::math::nspace::BlockSparseNSpace space_;
  unsigned nevent_ = 0;
  Eigen::VectorXd p_;
  Eigen::Vector3d x0_;
  Eigen::Matrix3d rot_ = Eigen::Matrix3d::Identity();
  calin::simulation::detector_efficiency::ACTEffectiveBandwidth* effective_bandwidth_ = nullptr;
};

} } } // namespace calin::simulation::ray_processor
