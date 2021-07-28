/*

   calin/iact_data/instrument_layout.hpp -- Stephen Fegan -- 2017-09-15

   Generic instrument layout functions

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

#include <Eigen/Core>

#include <calin_global_definitions.hpp>
#include <calin_global_config.hpp>
#include <math/regular_grid.hpp>
#include <iact_data/instrument_layout.pb.h>

namespace calin { namespace iact_data { namespace instrument_layout {

calin::math::regular_grid::Grid* make_grid_from_instrument_layout(
  const calin::ix::iact_data::instrument_layout::CameraLayout& camera_layout,
  double grid_plate_scale = 1.0);

void compute_camera_and_module_outlines(
  calin::ix::iact_data::instrument_layout::CameraLayout* camera_layout);

void map_channels_using_grid(
  Eigen::VectorXi& map,
  const calin::ix::iact_data::instrument_layout::CameraLayout& from,
  const calin::ix::iact_data::instrument_layout::CameraLayout& to);

void map_channels_using_from_coordinates(
  Eigen::VectorXi& map,
  const std::vector<double>& from_x, const std::vector<double>& from_y,
  const calin::ix::iact_data::instrument_layout::CameraLayout& to,
  double tolerance = 0.1);

void map_channels_using_from_coordinates(
  Eigen::VectorXi& map,
  const calin::ix::iact_data::instrument_layout::CameraLayout& from,
  const calin::ix::iact_data::instrument_layout::CameraLayout& to,
  double tolerance = 0.1);

#ifndef SWIG
calin::ix::iact_data::instrument_layout::CameraLayout* reduce_camera_channels(
  const calin::ix::iact_data::instrument_layout::CameraLayout& layout_in,
  const unsigned* channel_id, unsigned nchannel_id, bool recenter = false);

calin::ix::iact_data::instrument_layout::CameraLayout* reduce_camera_modules(
  const calin::ix::iact_data::instrument_layout::CameraLayout& layout_in,
  const unsigned* module_id, unsigned nmodule_id, bool recenter = false);

calin::ix::iact_data::instrument_layout::OutlinePolygon* channel_outline(
  const calin::ix::iact_data::instrument_layout::CameraLayout& camera_layout,
  const unsigned* channel_id, unsigned nchannel_id);
#endif

inline calin::ix::iact_data::instrument_layout::CameraLayout* reduce_camera_channels(
  const calin::ix::iact_data::instrument_layout::CameraLayout& layout_in,
  const Eigen::VectorXi& channel_id, bool recenter = false)
{
  return reduce_camera_channels(layout_in, reinterpret_cast<const unsigned*>(channel_id.data()), channel_id.size(), recenter);
}

inline calin::ix::iact_data::instrument_layout::CameraLayout* reduce_camera_modules(
  const calin::ix::iact_data::instrument_layout::CameraLayout& layout_in,
  const Eigen::VectorXi& module_id, bool recenter = false)
{
  return reduce_camera_modules(layout_in, reinterpret_cast<const unsigned*>(module_id.data()), module_id.size(), recenter);
}

inline calin::ix::iact_data::instrument_layout::OutlinePolygon* channel_outline(
  const calin::ix::iact_data::instrument_layout::CameraLayout& camera_layout,
  const Eigen::VectorXi& channel_id)
{
  return channel_outline(camera_layout, reinterpret_cast<const unsigned*>(channel_id.data()), channel_id.size());
}

#ifndef SWIG
calin::ix::iact_data::instrument_layout::CameraLayout*
camera_layout(calin::ix::iact_data::instrument_layout::CameraLayout::CameraType camera_type,
  calin::ix::iact_data::instrument_layout::CameraLayout* layout = nullptr);
#else
calin::ix::iact_data::instrument_layout::CameraLayout*
camera_layout(calin::ix::iact_data::instrument_layout::CameraLayout::CameraType camera_type);
#endif

} } } // namespace calin::iact_data::instrument_layout
