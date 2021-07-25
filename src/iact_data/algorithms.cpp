/*

   calin/iact_data/algorithms.cpp -- Stephen Fegan -- 2021-07-23

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

#include <algorithm>

#include <util/log.hpp>
#include <iact_data/algorithms.hpp>

using namespace calin::util::log;

unsigned calin::iact_data::algorithms::find_channel_islands(
  const calin::ix::iact_data::instrument_layout::CameraLayout& camera,
  const int* channel_id, unsigned nchannel_id, int* channel_island_id, int* island_count)
{
  const int* last_channel_id = channel_id+nchannel_id;
  std::fill(channel_island_id, channel_island_id+nchannel_id, -1);

  unsigned nisland = 0;
  unsigned nstack = 0;
  unsigned ichannel_id = 0;
  while(ichannel_id<nchannel_id or nstack>0) {
    int cid;
    if(nstack) {
      --nstack;
      cid = island_count[nisland + nstack];
    } else {
      cid = channel_id[ichannel_id];
      if(channel_island_id[ichannel_id] >= 0) {
        ++ichannel_id;
        continue;
      } else {
        channel_island_id[ichannel_id] = nisland;
        island_count[nisland] = 1;
        ++nisland;
        ++ichannel_id;
      }
    }
    for(auto ncid : camera.channel(cid).neighbour_channel_indexes()) {
      if(ncid > channel_id[ichannel_id-1]) {
        const auto* found_channel_id =
          std::lower_bound(channel_id+ichannel_id, last_channel_id, ncid);
        if(found_channel_id<last_channel_id and *found_channel_id==ncid) {
          unsigned ifound_channel_id = found_channel_id-channel_id;
          if(channel_island_id[ifound_channel_id] == -1) {
            int ciid = nisland-1;
            channel_island_id[ifound_channel_id] = ciid;
            island_count[ciid] += 1;
            island_count[nisland + nstack] = ncid;
            ++nstack;
          }
        }
      }
    }
  }
  return nisland;
}

unsigned calin::iact_data::algorithms::find_channel_islands(
  const calin::ix::iact_data::instrument_layout::CameraLayout& camera,
  const Eigen::VectorXi& channel_id,
  Eigen::VectorXi& channel_island_id, Eigen::VectorXi& island_count)
{
  channel_island_id.resize(channel_id.size());
  island_count.resize(channel_id.size());
  unsigned nisland = find_channel_islands(camera, channel_id.data(), channel_id.size(),
    channel_island_id.data(), island_count.data());
  island_count = island_count.head(nisland);
  return nisland;
}
