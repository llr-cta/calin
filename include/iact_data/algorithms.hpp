/*

  calin/iact_data/algorithms.hpp -- Stephen Fegan -- 2021-07-23

  IACT specific algorithms

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

#include <string>
#include <Eigen/Dense>

#include <calin_global_definitions.hpp>
#include <calin_global_config.hpp>
#include <iact_data/instrument_layout.hpp>
#include <iact_data/instrument_layout.pb.h>

namespace calin { namespace iact_data { namespace algorithms {

/*

  Find "islands" of contigous channels and return the number of channels in
  each island.

  INPUT :
  camera : the IACT camera layout specifying the graph of channel neighbours
  channel_id : array of channel IDs that make up the "land mass" of the islands.
      The array *must* contain channels IDs in monotonic ascending order.
  nchannel_id : number of channels in the channel_id array

  RETURNS :
  nisland : the number of islands found.

  OUTPUT :
  channel_island_id : island number to which each channel_id belongs. Must be
      pre-allocated to size of at least nchannel_id
  island_count : count of number of channels in each island. Must be
      pre-allocated to size of at least nchannel_id. Only the first "nisland"
      entries of this array will contain useful values (although all entries may
      be modified by the algorithm).

  */
#ifndef SWIG
unsigned find_channel_islands(
  const calin::ix::iact_data::instrument_layout::CameraLayout& camera,
  const int* channel_id, unsigned nchannel_id, int* channel_island_id, int* island_count);
#endif

unsigned find_channel_islands(
  const calin::ix::iact_data::instrument_layout::CameraLayout& camera,
  const Eigen::VectorXi& channel_id,
  Eigen::VectorXi& channel_island_id, Eigen::VectorXi& island_count);

} } } // namespace calin::iact_data::algorithms
