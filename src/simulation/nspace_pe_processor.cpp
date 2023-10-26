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
  double x, double y, double ux, double uy, double t, double pe_weight)
{
  switch(config_.axis_variables()) {
  case calin::ix::simulation::pe_processor::XY:
  default:
    p_[0] = x;
    p_[1] = y;
    break;
  case calin::ix::simulation::pe_processor::UXUY:
    p_[0] = ux;
    p_[1] = uy;
    break;
  case calin::ix::simulation::pe_processor::T:
    p_[0] = t - config_.t_origin();
    break;
  case calin::ix::simulation::pe_processor::XY_UXUY:
    p_[0] = x;
    p_[1] = y;
    p_[2] = ux;
    p_[3] = uy;
    break;
  case calin::ix::simulation::pe_processor::XY_T:
    p_[0] = x;
    p_[1] = y;
    p_[2] = t - config_.t_origin();
    break;
  case calin::ix::simulation::pe_processor::UXUY_T:
    p_[0] = ux;
    p_[1] = uy;
    p_[2] = t - config_.t_origin();
    break;
  case calin::ix::simulation::pe_processor::XY_UXUY_T:
    p_[0] = x;
    p_[1] = y;
    p_[2] = ux;
    p_[3] = uy;
    p_[4] = t - config_.t_origin();
    break;
  case calin::ix::simulation::pe_processor::SCOPE_XY:
    p_[0] = scope_id;
    p_[1] = x;
    p_[2] = y;
    break;
  case calin::ix::simulation::pe_processor::SCOPE_UXUY:
    p_[0] = scope_id;
    p_[1] = ux;
    p_[2] = uy;
    break;
  case calin::ix::simulation::pe_processor::SCOPE_T:
    p_[0] = scope_id;
    p_[1] = t - config_.t_origin();
    break;
  case calin::ix::simulation::pe_processor::SCOPE_XY_UXUY:
    p_[0] = scope_id;
    p_[1] = x;
    p_[2] = y;
    p_[3] = ux;
    p_[4] = uy;
  case calin::ix::simulation::pe_processor::SCOPE_XY_T:
    p_[0] = scope_id;
    p_[1] = x;
    p_[2] = y;
    p_[3] = t - config_.t_origin();
    break;
  case calin::ix::simulation::pe_processor::SCOPE_UXUY_T:
    p_[0] = scope_id;
    p_[1] = ux;
    p_[2] = uy;
    p_[3] = t - config_.t_origin();
    break;
  case calin::ix::simulation::pe_processor::SCOPE_XY_UXUY_T:
    p_[0] = scope_id;
    p_[1] = x;
    p_[2] = y;
    p_[3] = ux;
    p_[4] = uy;
    p_[5] = t - config_.t_origin();
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
  double uxuy_radius = 0.5*config_.uxuy_diameter()*M_PI/180.0;

  switch(config_.axis_variables()) {
  case calin::ix::simulation::pe_processor::SCOPE_XY:
  case calin::ix::simulation::pe_processor::SCOPE_UXUY:
  case calin::ix::simulation::pe_processor::SCOPE_T:
  case calin::ix::simulation::pe_processor::SCOPE_XY_UXUY:
  case calin::ix::simulation::pe_processor::SCOPE_XY_T:
  case calin::ix::simulation::pe_processor::SCOPE_UXUY_T:
  case calin::ix::simulation::pe_processor::SCOPE_XY_UXUY_T:
    axes.push_back({-0.5, config_.num_scopes()-0.5, config_.num_scopes()});
    break;
  case calin::ix::simulation::pe_processor::XY:
  case calin::ix::simulation::pe_processor::UXUY:
  case calin::ix::simulation::pe_processor::T:
  case calin::ix::simulation::pe_processor::XY_UXUY:
  case calin::ix::simulation::pe_processor::XY_T:
  case calin::ix::simulation::pe_processor::UXUY_T:
  case calin::ix::simulation::pe_processor::XY_UXUY_T:
  default:
    // do nothing
    break;
  };

  switch(config_.axis_variables()) {
  case calin::ix::simulation::pe_processor::XY:
  case calin::ix::simulation::pe_processor::XY_UXUY:
  case calin::ix::simulation::pe_processor::XY_T:
  case calin::ix::simulation::pe_processor::XY_UXUY_T:
  case calin::ix::simulation::pe_processor::SCOPE_XY:
  case calin::ix::simulation::pe_processor::SCOPE_XY_UXUY:
  case calin::ix::simulation::pe_processor::SCOPE_XY_T:
  case calin::ix::simulation::pe_processor::SCOPE_XY_UXUY_T:
  default:
    axes.push_back({-xy_radius, xy_radius, config_.xy_num_bins()});
    axes.push_back({-xy_radius, xy_radius, config_.xy_num_bins()});
    break;
  case calin::ix::simulation::pe_processor::UXUY:
  case calin::ix::simulation::pe_processor::T:
  case calin::ix::simulation::pe_processor::UXUY_T:
  case calin::ix::simulation::pe_processor::SCOPE_UXUY:
  case calin::ix::simulation::pe_processor::SCOPE_T:
  case calin::ix::simulation::pe_processor::SCOPE_UXUY_T:
    // do nothing
    break;
  };

  switch(config_.axis_variables()) {
  case calin::ix::simulation::pe_processor::UXUY:
  case calin::ix::simulation::pe_processor::XY_UXUY:
  case calin::ix::simulation::pe_processor::UXUY_T:
  case calin::ix::simulation::pe_processor::XY_UXUY_T:
  case calin::ix::simulation::pe_processor::SCOPE_UXUY:
  case calin::ix::simulation::pe_processor::SCOPE_XY_UXUY:
  case calin::ix::simulation::pe_processor::SCOPE_UXUY_T:
  case calin::ix::simulation::pe_processor::SCOPE_XY_UXUY_T:
    axes.push_back({-uxuy_radius, uxuy_radius, config_.uxuy_num_bins()});
    axes.push_back({-uxuy_radius, uxuy_radius, config_.uxuy_num_bins()});
    break;
  case calin::ix::simulation::pe_processor::T:
  case calin::ix::simulation::pe_processor::XY:
  case calin::ix::simulation::pe_processor::XY_T:
  case calin::ix::simulation::pe_processor::SCOPE_XY:
  case calin::ix::simulation::pe_processor::SCOPE_T:
  case calin::ix::simulation::pe_processor::SCOPE_XY_T:
  default:
    // do nothing
    break;
  };

  switch(config_.axis_variables()) {
  case calin::ix::simulation::pe_processor::T:
  case calin::ix::simulation::pe_processor::XY_T:
  case calin::ix::simulation::pe_processor::UXUY_T:
  case calin::ix::simulation::pe_processor::XY_UXUY_T:
  case calin::ix::simulation::pe_processor::SCOPE_T:
  case calin::ix::simulation::pe_processor::SCOPE_XY_T:
  case calin::ix::simulation::pe_processor::SCOPE_UXUY_T:
  case calin::ix::simulation::pe_processor::SCOPE_XY_UXUY_T:
    axes.push_back({0.0, config_.t_duration(), config_.t_num_bins()});
    break;
  case calin::ix::simulation::pe_processor::XY:
  case calin::ix::simulation::pe_processor::UXUY:
  case calin::ix::simulation::pe_processor::XY_UXUY:
  case calin::ix::simulation::pe_processor::SCOPE_XY:
  case calin::ix::simulation::pe_processor::SCOPE_UXUY:
  case calin::ix::simulation::pe_processor::SCOPE_XY_UXUY:
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
  config.set_uxuy_diameter(102.4);
  config.set_uxuy_num_bins(1024);
  config.set_t_duration(100000);
  config.set_t_num_bins(100000);
  config.set_clear_at_new_event(false);
  return config;
}
