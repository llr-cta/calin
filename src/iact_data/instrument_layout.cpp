/*

   calin/iact_data/instrument_layout.cpp -- Stephen Fegan -- 2017-09-15

   Generic camera layout functions

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <cmath>
#include <vector>
#include <stdexcept>
#include <cassert>

#include <google/protobuf/repeated_field.h>

#include <math/regular_grid.hpp>
#include <iact_data/instrument_layout.hpp>
#include <io/log.hpp>

using namespace calin::math::regular_grid;
using namespace calin::ix::iact_data::instrument_layout;
using namespace calin::io::log;

calin::math::regular_grid::Grid* calin::iact_data::instrument_layout::
make_grid_from_instrument_layout(
  const calin::ix::iact_data::instrument_layout::CameraLayout* camera_layout,
  double grid_plate_scale)
{
  switch(camera_layout->pixel_grid_layout())
  {
  case calin::ix::iact_data::instrument_layout::CameraLayout::HEX_GRID:
    return new calin::math::regular_grid::HexGrid(
      camera_layout->pixel_grid_spacing()*grid_plate_scale,
      camera_layout->pixel_grid_rotation()/180.0*M_PI,
      camera_layout->pixel_grid_offset_x()*grid_plate_scale,
      camera_layout->pixel_grid_offset_y()*grid_plate_scale);
  default:
    return nullptr;
  }
}

namespace {

  void compute_and_store_boundary(calin::math::regular_grid::Grid* grid,
    const std::vector<unsigned>& grid_ids,
    google::protobuf::RepeatedField<google::protobuf::uint32>* indexes,
    google::protobuf::RepeatedField<double>* vx,
    google::protobuf::RepeatedField<double>* vy)
  {
    indexes->Clear();
    vx->Clear();
    vy->Clear();
    auto boundary = grid->compute_region_boundary(grid_ids);
    unsigned ncurve = grid->num_bounday_curves(boundary);
    for(unsigned icurve=0; icurve<ncurve; icurve++) {
      Eigen::VectorXd xv;
      Eigen::VectorXd yv;
      grid->extract_bounday_curve(xv, yv, boundary, icurve, false);
      for(unsigned ixv=0; ixv<xv.size(); ixv++)vx->Add(xv[ixv]);
      for(unsigned iyv=0; iyv<yv.size(); iyv++)vy->Add(yv[iyv]);
      indexes->Add(vx->size());
    }
  }

} // anonymous namespace

void calin::iact_data::instrument_layout::compute_camera_and_module_outlines(
  calin::ix::iact_data::instrument_layout::CameraLayout* camera_layout)
{
  auto* grid = make_grid_from_instrument_layout(camera_layout);
  if(grid == nullptr)
    throw std::runtime_error("compute_camera_and_module_outlines: camera must "
      "have standard grid.");

  std::vector<unsigned> camera_grid_ids;
  for(auto chan : camera_layout->channel())
    if(chan.pixel_index() != -1)
      camera_grid_ids.push_back(chan.pixel_grid_index());
  compute_and_store_boundary(grid, camera_grid_ids,
    camera_layout->mutable_outline_polygon_vertex_index(),
    camera_layout->mutable_outline_polygon_vertex_x(),
    camera_layout->mutable_outline_polygon_vertex_y());

  for(unsigned imodule=0; imodule<camera_layout->module_size(); imodule++)
  {
    auto* module_layout = camera_layout->mutable_module(imodule);
    std::vector<unsigned> module_grid_ids;
    for(auto ichan : module_layout->channels_in_module()) {
      auto chan = camera_layout->channel(ichan);
      if(chan.pixel_index() != -1)
        module_grid_ids.push_back(chan.pixel_grid_index());
    }
    compute_and_store_boundary(grid, module_grid_ids,
      module_layout->mutable_outline_polygon_vertex_index(),
      module_layout->mutable_outline_polygon_vertex_x(),
      module_layout->mutable_outline_polygon_vertex_y());
  }
}
