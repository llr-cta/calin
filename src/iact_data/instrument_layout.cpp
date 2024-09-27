/*

   calin/iact_data/instrument_layout.cpp -- Stephen Fegan -- 2017-09-15

   Generic camera layout functions

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

#include <cmath>
#include <vector>
#include <stdexcept>
#include <cassert>

#include <google/protobuf/repeated_field.h>

#include <math/regular_grid.hpp>
#include <math/special.hpp>
#include <iact_data/instrument_layout.hpp>
#include <iact_data/nectarcam_layout.hpp>
#include <iact_data/lstcam_layout.hpp>
#include <util/log.hpp>

using namespace calin::math::regular_grid;
using namespace calin::ix::iact_data::instrument_layout;
using namespace calin::util::log;

using calin::math::special::SQR;

calin::math::regular_grid::Grid* calin::iact_data::instrument_layout::
make_grid_from_instrument_layout(
  const calin::ix::iact_data::instrument_layout::CameraLayout& camera_layout,
  double grid_plate_scale)
{
  switch(camera_layout.pixel_grid_layout())
  {
  case calin::ix::iact_data::instrument_layout::CameraLayout::HEX_GRID:
    return new calin::math::regular_grid::HexGrid(
      camera_layout.pixel_grid_spacing()*grid_plate_scale,
      camera_layout.pixel_grid_rotation()/180.0*M_PI,
      camera_layout.pixel_grid_offset_x()*grid_plate_scale,
      camera_layout.pixel_grid_offset_y()*grid_plate_scale);
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
      grid->extract_boundary_curve(xv, yv, boundary, icurve, false);
      for(unsigned ixv=0; ixv<xv.size(); ixv++)vx->Add(xv[ixv]);
      for(unsigned iyv=0; iyv<yv.size(); iyv++)vy->Add(yv[iyv]);
      indexes->Add(vx->size());
    }
  }

} // anonymous namespace

void calin::iact_data::instrument_layout::compute_camera_and_module_outlines(
  calin::ix::iact_data::instrument_layout::CameraLayout* camera_layout)
{
  auto* grid = make_grid_from_instrument_layout(*camera_layout);
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

  for(int imodule=0; imodule<camera_layout->module_size(); imodule++)
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

  double x = camera_layout->outline_polygon_vertex_x(0);
  double y = camera_layout->outline_polygon_vertex_y(0);
  camera_layout->set_camera_boundary_box_left(x);
  camera_layout->set_camera_boundary_box_right(x);
  camera_layout->set_camera_boundary_box_bottom(y);
  camera_layout->set_camera_boundary_box_top(y);
  camera_layout->set_camera_boundary_maxabs_xy(std::max(std::abs(x),std::abs(y)));
  camera_layout->set_camera_boundary_max_r(x*x+y*y);
  for(int ivertex=1;ivertex<camera_layout->outline_polygon_vertex_x_size();++ivertex) {
    double x = camera_layout->outline_polygon_vertex_x(ivertex);
    double y = camera_layout->outline_polygon_vertex_y(ivertex);
    camera_layout->set_camera_boundary_box_left(
      std::min(camera_layout->camera_boundary_box_left(),x));
    camera_layout->set_camera_boundary_box_right(
      std::max(camera_layout->camera_boundary_box_right(),x));
    camera_layout->set_camera_boundary_box_bottom(
      std::min(camera_layout->camera_boundary_box_bottom(),y));
    camera_layout->set_camera_boundary_box_top(
      std::max(camera_layout->camera_boundary_box_top(),y));
    camera_layout->set_camera_boundary_maxabs_xy(
      std::max(camera_layout->camera_boundary_maxabs_xy(),std::max(std::abs(x),std::abs(y))));
    camera_layout->set_camera_boundary_max_r(
      std::max(camera_layout->camera_boundary_max_r(),x*x+y*y));
  }
  camera_layout->set_camera_boundary_max_r(std::sqrt(camera_layout->camera_boundary_max_r()));

  delete grid;
}

void calin::iact_data::instrument_layout::map_channels_using_grid(
  Eigen::VectorXi& map,
  const calin::ix::iact_data::instrument_layout::CameraLayout& from,
  const calin::ix::iact_data::instrument_layout::CameraLayout& to)
{
  std::map<unsigned, int> grid_map;
  for(int ichan=0; ichan<from.channel_size(); ichan++) {
    const auto& chan = from.channel(ichan);
    if(chan.pixel_grid_index() >= 0)
      grid_map[chan.pixel_grid_index()] = chan.channel_index();
  }
  map.resize(to.channel_size());
  for(int ichan=0; ichan<to.channel_size(); ichan++) {
    const auto& chan = to.channel(ichan);
    if(chan.pixel_grid_index() < 0)
      map[ichan] = -1;
    else if(grid_map.find(chan.pixel_grid_index()) == grid_map.end())
      map[ichan] = -1;
    else
      map[ichan] = grid_map[chan.pixel_grid_index()];
  }
}

void calin::iact_data::instrument_layout::map_channels_using_from_coordinates(
  Eigen::VectorXi& map,
  const std::vector<double>& from_x, const std::vector<double>& from_y,
  const calin::ix::iact_data::instrument_layout::CameraLayout& to,
  double tolerance)
{
  tolerance *= tolerance;
  if(from_x.size() != from_y.size())
    throw std::runtime_error("X and Y coordinate vectors must be same size.");
  map.resize(to.channel_size());
  for(int ichan=0; ichan<to.channel_size(); ichan++) {
    const auto& chan = to.channel(ichan);
    map[ichan] = -1;
    for(unsigned jchan=0; jchan<from_x.size(); jchan++) {
      double d2 = SQR(chan.x()-from_x[jchan])+SQR(chan.y()-from_y[jchan]);
      if(d2 < tolerance) {
        map[ichan] = jchan;
        break;
      }
    }
  }
}

void calin::iact_data::instrument_layout::map_channels_using_from_coordinates(
  Eigen::VectorXi& map,
  const calin::ix::iact_data::instrument_layout::CameraLayout& from,
  const calin::ix::iact_data::instrument_layout::CameraLayout& to,
  double tolerance)
{
  std::vector<double> from_x(from.channel_size());
  std::vector<double> from_y(from.channel_size());
  for(int ichan=0; ichan<from.channel_size(); ichan++) {
    const auto& chan = from.channel(ichan);
    from_x[ichan] = chan.x();
    from_y[ichan] = chan.y();
  }
  map_channels_using_from_coordinates(map, from_x, from_y, to, tolerance);
}

calin::ix::iact_data::instrument_layout::CameraLayout*
calin::iact_data::instrument_layout::reorder_camera_channels(
  const calin::ix::iact_data::instrument_layout::CameraLayout& in,
  const unsigned* channel_id, unsigned nchannel_id, bool recenter)
{
  std::vector<int> in_to_out(in.channel_size(), -1);

  double min_x = 0;
  double max_x = 0;
  double min_y = 0;
  double max_y = 0;
  for(unsigned i=0;i<nchannel_id;++i) {
    if(channel_id[i] >= unsigned(in.channel_size())) {
      throw std::runtime_error("Channel ID out of range: " +
        std::to_string(channel_id[i]));
    }
    if(in_to_out[channel_id[i]] != -1) {
      throw std::runtime_error("Duplicate channel ID: " +
        std::to_string(channel_id[i]));
    }
    in_to_out[channel_id[i]] = i;
    const auto& ic = in.channel(channel_id[i]);
    if(i==0) {
      min_x = max_x = ic.x();
      min_y = max_y = ic.y();
    } else {
      min_x = std::min(min_x, ic.x());
      max_x = std::max(max_x, ic.x());
      min_y = std::min(min_y, ic.y());
      max_y = std::max(max_y, ic.y());
    }
  }

  double shift_x = recenter ? 0.5*(min_x+max_x) : 0.0;
  double shift_y = recenter ? 0.5*(min_y+max_y) : 0.0;

  auto* out = new calin::ix::iact_data::instrument_layout::CameraLayout();

  out->set_camera_type(in.camera_type());
  out->set_camera_number(in.camera_number());

  for(unsigned ichan=0;ichan<nchannel_id;++ichan) {
    const auto& ic = in.channel(channel_id[ichan]);
    auto* oc = out->add_channel();

    oc->set_channel_index(ichan);
    oc->set_pixel_index(ic.pixel_index());
    oc->set_pixel_grid_index(ic.pixel_grid_index());
    oc->set_channel_set_index(ic.channel_set_index());
    oc->set_module_index(-1);
    oc->set_module_channel_index(-1);
    oc->set_pixel_grid_u(ic.pixel_grid_u());
    oc->set_pixel_grid_v(ic.pixel_grid_v());
    oc->set_pixel_spiral_index(ic.pixel_spiral_index());
    oc->set_x(ic.x() - shift_x);
    oc->set_y(ic.y() - shift_y);
    oc->set_diameter(ic.diameter());
    oc->set_geometric_area(ic.geometric_area());

    oc->set_is_boundary_pixel(ic.is_boundary_pixel());
    for(auto nci : ic.neighbour_channel_indexes()) {
      if(in_to_out[nci] != -1) {
        oc->add_neighbour_channel_indexes(in_to_out[nci]);
      } else {
        oc->set_is_boundary_pixel(true);
      }
    }
    if(oc->is_boundary_pixel()) {
      out->add_boundary_pixel_channel_index(ichan);
    }
    for(auto index : ic.outline_polygon_vertex_index())
      oc->add_outline_polygon_vertex_index(index);
    for(auto x : ic.outline_polygon_vertex_x())
      oc->add_outline_polygon_vertex_x(x - shift_x);
    for(auto y : ic.outline_polygon_vertex_y())
      oc->add_outline_polygon_vertex_y(y - shift_y);
  }

  out->set_pixel_grid_layout(in.pixel_grid_layout());
  out->set_pixel_grid_spacing(in.pixel_grid_spacing());
  out->set_pixel_grid_rotation(in.pixel_grid_rotation());
  out->set_pixel_grid_cos_rotation(in.pixel_grid_cos_rotation());
  out->set_pixel_grid_sin_rotation(in.pixel_grid_sin_rotation());
  out->set_pixel_grid_offset_x(in.pixel_grid_offset_x() - shift_x);
  out->set_pixel_grid_offset_y(in.pixel_grid_offset_y() - shift_y);
  out->set_pixel_grid_geometric_area(in.pixel_grid_geometric_area());

  for(auto i : in.pixel_channel_index())
    out->add_pixel_channel_index(i==-1 ? -1 : in_to_out[i]);
  for(auto i : in.pixel_spiral_channel_index())
    out->add_pixel_spiral_channel_index(i==-1 ? -1 : in_to_out[i]);

  out->set_adc_gains(in.adc_gains());
  out->set_adc_bits(in.adc_bits());
  out->set_can_read_waveforms(in.can_read_waveforms());
  out->set_can_read_charges(in.can_read_charges());
  out->set_can_read_peak_sample(in.can_read_peak_sample());

  compute_camera_and_module_outlines(out);
  return out;
}

calin::ix::iact_data::instrument_layout::CameraLayout*
calin::iact_data::instrument_layout::reorder_camera_modules(
  const calin::ix::iact_data::instrument_layout::CameraLayout& in,
  const unsigned* module_id, unsigned nmodule_id, bool recenter)
{
  std::vector<int> in_to_out(in.module_size(), -1);
  std::vector<unsigned> channel_id;

  for(unsigned i=0;i<nmodule_id;++i) {
    if(module_id[i] >= unsigned(in.module_size())) {
      throw std::runtime_error("Module ID out of range: " +
        std::to_string(module_id[i]));
    }
    if(in_to_out[module_id[i]] != -1) {
      throw std::runtime_error("Duplicate module ID: " +
        std::to_string(module_id[i]));
    }
    in_to_out[module_id[i]] = i;
    for(auto ic : in.module(module_id[i]).channels_in_module()) {
      channel_id.push_back(ic);
    }
  }

  Eigen::VectorXi eigen_channel_id = calin::std_to_eigenvec_unsigned(channel_id);
  auto* out = reorder_camera_channels(in, eigen_channel_id, recenter);

  double shift_x = in.pixel_grid_offset_x() - out->pixel_grid_offset_x();
  double shift_y = in.pixel_grid_offset_y() - out->pixel_grid_offset_y();

  for(unsigned imod=0, ichan=0;imod<nmodule_id;++imod) {
    const auto& im = in.module(module_id[imod]);
    auto* om = out->add_module();

    om->set_module_index(imod);
    om->set_module_grid_index(im.module_grid_index());
    om->set_module_grid_u(im.module_grid_u());
    om->set_module_grid_v(im.module_grid_v());
    om->set_module_grid_i(im.module_grid_i());
    om->set_module_grid_j(im.module_grid_j());
    om->set_module_spiral_index(im.module_spiral_index());
    for(int imodchan=0; imodchan<im.channels_in_module_size(); imodchan++) {
      out->mutable_channel(ichan)->set_module_index(imod);
      out->mutable_channel(ichan)->set_module_channel_index(imodchan);
      om->add_channels_in_module(ichan++);
    }
    om->set_x(im.x() - shift_x);
    om->set_y(im.y() - shift_y);

    for(auto index : im.outline_polygon_vertex_index())
      om->add_outline_polygon_vertex_index(index);
    for(auto x : im.outline_polygon_vertex_x())
      om->add_outline_polygon_vertex_x(x - shift_x);
    for(auto y : im.outline_polygon_vertex_y())
      om->add_outline_polygon_vertex_y(y - shift_y);
  }

  for(auto i : in.module_spiral_module_index())
    out->module_spiral_module_index(i==-1 ? -1 : in_to_out[i]);

  compute_camera_and_module_outlines(out);
  return out;
}

calin::ix::iact_data::instrument_layout::OutlinePolygon*
calin::iact_data::instrument_layout::channel_outline(
  const calin::ix::iact_data::instrument_layout::CameraLayout& camera_layout,
  const unsigned* channel_id, unsigned nchannel_id)
{
  auto* grid = make_grid_from_instrument_layout(camera_layout);
  if(grid == nullptr)
    throw std::runtime_error("channel_outline: camera must "
      "have standard grid.");

  std::vector<unsigned> camera_grid_ids;
  std::vector<bool> channel_ids_selected(camera_layout.channel_size(), false);
  for(unsigned ichannel_id=0; ichannel_id<nchannel_id; ichannel_id++) {
    unsigned ichan = channel_id[ichannel_id];
    if(ichan >= unsigned(camera_layout.channel_size())) {
      delete grid;
      throw std::runtime_error("Channel ID out of range: "
        + std::to_string(ichan));
    }
    if(camera_layout.channel(ichan).pixel_grid_index() == -1) {
      delete grid;
      throw std::runtime_error("Channel has no pixel grid ID: "
        + std::to_string(ichan));
    }
    if(channel_ids_selected[ichan] == false) {
      channel_ids_selected[ichan] = true;
      camera_grid_ids.push_back(camera_layout.channel(ichan).pixel_grid_index());
    }
  }

  auto* outline = new calin::ix::iact_data::instrument_layout::OutlinePolygon();

  compute_and_store_boundary(grid, camera_grid_ids,
    outline->mutable_outline_polygon_vertex_index(),
    outline->mutable_outline_polygon_vertex_x(),
    outline->mutable_outline_polygon_vertex_y());

  delete grid;
  return outline;
}

calin::ix::iact_data::instrument_layout::CameraLayout*
calin::iact_data::instrument_layout::camera_layout(
  calin::ix::iact_data::instrument_layout::CameraLayout::CameraType camera_type,
  calin::ix::iact_data::instrument_layout::CameraLayout* layout)
{
  switch(camera_type) {
  case calin::ix::iact_data::instrument_layout::CameraLayout::NECTARCAM:
    return calin::iact_data::nectarcam_layout::nectarcam_layout(layout);
  case calin::ix::iact_data::instrument_layout::CameraLayout::LSTCAM:
    return calin::iact_data::lstcam_layout::lstcam_layout(layout);
  case calin::ix::iact_data::instrument_layout::CameraLayout::NECTARCAM_TESTBENCH_19CHANNEL:
    return calin::iact_data::nectarcam_layout::nectarcam_19module_layout(layout);
  case calin::ix::iact_data::instrument_layout::CameraLayout::NECTARCAM_TESTBENCH_61CHANNEL:
    return calin::iact_data::nectarcam_layout::nectarcam_61module_layout(layout);
  case calin::ix::iact_data::instrument_layout::CameraLayout::NO_CAMERA:
  default:
    return nullptr;
  }
}
