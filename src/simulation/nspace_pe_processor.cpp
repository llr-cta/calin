/*

   calin/simulation/nspace_pe_processor.cpp -- Stephen Fegan -- 2022-06-13

   Process PEs into an NSpace.

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

#include <util/log.hpp>
#include <simulation/nspace_pe_processor.hpp>

using namespace calin::util::log;
using namespace calin::simulation::pe_processor;

NSpacePEProcessor::
NSpacePEProcessor(const calin::ix::simulation::pe_processor::NSpacePEProcessorConfig& config,
    unsigned nspace_log2_block_size):
  PEProcessor(),
  config_(config), space_(nspace_axes(), nspace_log2_block_size), p_(space_.naxes())
{
  // nothing to see here
}

NSpacePEProcessor::~NSpacePEProcessor()
{
  // nothing to see here
}

void NSpacePEProcessor::start_processing()
{
  if(config_.clear_at_new_event()) {
    clear();
  }
}

void NSpacePEProcessor::process_focal_plane_hit(unsigned scope_id, int pixel_id,
  double x, double y, double t0, double pe_weight)
{
  switch(config_.axis_variables()) {
  case calin::ix::simulation::pe_processor::XY:
  default:
    p_[0] = x;
    p_[1] = y;
    break;
  case calin::ix::simulation::pe_processor::T:
    p_[0] = t0*1e9;
    break;
  case calin::ix::simulation::pe_processor::XY_T:
    p_[0] = x;
    p_[1] = y;
    p_[2] = t0*1e9;
    break;
  case calin::ix::simulation::pe_processor::SCOPE_XY:
    p_[0] = scope_id;
    p_[1] = x;
    p_[2] = y;
    break;
  case calin::ix::simulation::pe_processor::SCOPE_T:
    p_[0] = scope_id;
    p_[1] = t0*1e9;
    break;
  case calin::ix::simulation::pe_processor::SCOPE_XY_T:
    p_[0] = scope_id;
    p_[1] = x;
    p_[2] = y;
    p_[3] = t0*1e9;
    break;
  };
  space_.accumulate(p_, pe_weight);
}

void NSpacePEProcessor::finish_processing()
{
  ++nevent_;
}

void NSpacePEProcessor::clear()
{
  nevent_ = 0;
  space_.clear();
}

std::vector<calin::math::nspace::Axis>
NSpacePEProcessor::nspace_axes() const
{
  std::vector<calin::math::nspace::Axis> axes;

  double xy_radius = 0.5*config_.xy_diameter();

  switch(config_.axis_variables()) {
  case calin::ix::simulation::pe_processor::SCOPE_XY:
  case calin::ix::simulation::pe_processor::SCOPE_XY_T:
  case calin::ix::simulation::pe_processor::SCOPE_T:
    axes.push_back({-0.5, config_.num_scopes()-0.5, config_.num_scopes()});
    break;
  case calin::ix::simulation::pe_processor::T:
  case calin::ix::simulation::pe_processor::XY:
  case calin::ix::simulation::pe_processor::XY_T:
  default:
    // do nothing
    break;
  };

  switch(config_.axis_variables()) {
  case calin::ix::simulation::pe_processor::XY:
  case calin::ix::simulation::pe_processor::XY_T:
  case calin::ix::simulation::pe_processor::SCOPE_XY:
  case calin::ix::simulation::pe_processor::SCOPE_XY_T:
  default:
    axes.push_back({-xy_radius, xy_radius, config_.xy_num_bins()});
    axes.push_back({-xy_radius, xy_radius, config_.xy_num_bins()});
    break;
  case calin::ix::simulation::pe_processor::T:
  case calin::ix::simulation::pe_processor::SCOPE_T:
    // do nothing
    break;
  };

  switch(config_.axis_variables()) {
  case calin::ix::simulation::pe_processor::T:
  case calin::ix::simulation::pe_processor::XY_T:
  case calin::ix::simulation::pe_processor::SCOPE_XY_T:
  case calin::ix::simulation::pe_processor::SCOPE_T:
    axes.push_back({0.0, config_.t_duration(), config_.t_num_bins()});
    break;
  case calin::ix::simulation::pe_processor::XY:
  case calin::ix::simulation::pe_processor::SCOPE_XY:
  default:
    // do nothing
    break;
  };

  return axes;
}

calin::ix::simulation::pe_processor::NSpacePEProcessorConfig
NSpacePEProcessor::default_config()
{
  calin::ix::simulation::pe_processor::NSpacePEProcessorConfig config;
  config.set_axis_variables(calin::ix::simulation::pe_processor::XY);
  config.set_xy_diameter(512);
  config.set_xy_num_bins(1024);
  config.set_t_duration(100000);
  config.set_t_num_bins(100000);
  config.set_clear_at_new_event(false);
  return config;
}
