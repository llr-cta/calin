/*

   calin/iact_data/instrument_layout.hpp -- Stephen Fegan -- 2017-09-15

   Generic instrument layout functions

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

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

#include <calin_global_definitions.hpp>
#include <calin_global_config.hpp>
#include <math/regular_grid.hpp>
#include <iact_data/instrument_layout.pb.h>

namespace calin { namespace iact_data { namespace instrument_layout {

calin::math::regular_grid::Grid* make_grid_from_instrument_layout(
  const calin::ix::iact_data::instrument_layout::CameraLayout* camera_layout,
  double grid_plate_scale = 1.0);

void compute_camera_and_module_outlines(
  calin::ix::iact_data::instrument_layout::CameraLayout* camera_layout);

} } } // namespace calin::iact_data::instrument_layout
