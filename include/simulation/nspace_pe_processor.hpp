/*

   calin/simulation/nspace_pe_processor.hpp -- Stephen Fegan -- 2022-06-13

   Process PE hits into an NSpace

   Copyright 2022, Stephen Fegan <sfegan@llr.in2p3.fr>
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
#include <simulation/pe_processor.pb.h>
#include <simulation/pe_processor.hpp>

namespace calin { namespace simulation { namespace pe_processor {

class NSpacePEProcessor: public PEProcessor
{
public:
  NSpacePEProcessor(const calin::ix::simulation::pe_processor::NSpacePEProcessorConfig& config,
    unsigned nspace_log2_block_size = 0);
  virtual ~NSpacePEProcessor();

  void start_processing() override;
  void process_focal_plane_hit(unsigned scope_id, int pixel_id,
    double x, double y, double ux, double uy, double t0, double pe_weight) override;
  void finish_processing() override;

  void clear();
  unsigned nevent() const { return nevent_; }
  const calin::math::nspace::BlockSparseNSpace& nspace() const { return space_; }
  calin::math::nspace::BlockSparseNSpace& mutable_nspace() { return space_; }

  static calin::ix::simulation::pe_processor::NSpacePEProcessorConfig default_config();

protected:
  std::vector<calin::math::nspace::Axis> nspace_axes() const;
  calin::ix::simulation::pe_processor::NSpacePEProcessorConfig config_;
  calin::math::nspace::BlockSparseNSpace space_;
  Eigen::VectorXd p_;
  unsigned nevent_ = 0;
};

} } } // namespace calin::simulation::pe_processor
