/*

   calin/iact_data/nectarcam_layout->cpp -- Stephen Fegan -- 2016-06-06

   Camera layout for NectarCam

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <math/special.hpp>
#include <math/hex_array.hpp>
#include <iact_data/nectarcam_layout.hpp>
#include <io/log.hpp>

using calin::math::special::SQR;
using namespace calin::math::hex_array;
using namespace calin::ix::iact_data::instrument_layout;
using namespace calin::io::log;

namespace {

struct mod_info {
  mod_info(int imod_, int uc_, int vc_, double xc_, double yc_,
      int ix_, double iy_): imod(imod_), uc(uc_), vc(vc_), xc(xc_), yc(yc_),
    ix(ix_), iy(iy_) { /* nothing to see here */ }
  int imod;
  int uc;
  int vc;
  double xc;
  double yc;
  int ix;
  double iy;
};

std::vector<mod_info>
nectarcam_mod_map(unsigned nring, double radius, double spacing, double& rot)
{
  if(nring == 0 and radius <= 0)
    throw std::runtime_error("nectarcam_mod_map: At least one of nring or "
      "radius must be positive.");

  int u1 = 0;
  int v1 = 0;
  cluster_hexid_to_center_uv(1,1,u1,v1,false);
  double x1 = 0;
  double y1 = 0;
  uv_to_xy(u1,v1,x1,y1);
  rot = std::atan2(-y1,x1) + 30.0/180.0*M_PI;
  //LOG(INFO) << u1 << ' ' << v1 << ' ' << rot << ' ' << x1 << ' ' << y1;

  double crot = std::cos(rot);
  double srot = std::sin(rot);

  std::vector<mod_info> modvec;
  unsigned imod = 0;
  bool ring_has_module = true;
  for(unsigned iring=0; ring_has_module and (nring==0 or iring<=nring); iring++)
  {
    ring_has_module = false;
    for(unsigned iringmod=0; iringmod<std::max(1U,iring*6); iringmod++, imod++)
    {
      int uc;
      int vc;
      cluster_hexid_to_center_uv(imod,1,uc,vc,false);
      double xc;
      double yc;
      uv_to_xy_trans(uc,vc,xc,yc,crot,srot,spacing);
      if(radius>0 and xc*xc+yc*yc > radius*radius)continue;
      ring_has_module = true;
      int ix = (vc*2+uc*3)/7;
      double iy = -(vc+uc)+2.5*ix;
      modvec.emplace_back(imod, uc, vc, xc, yc, ix, iy);
    }
  }

  std::sort(modvec.begin(), modvec.end(),
    [](const mod_info& a, const mod_info& b) { if(a.ix<b.ix)return true;
      if(a.ix==b.ix)return a.iy<b.iy;
      return false; });

  return modvec;
}

CameraLayout* nectarcam_general_layout(CameraLayout* layout,
  const std::vector<mod_info>& modvec, CameraLayout::CameraType camera_type,
  double spacing, double rot)
{
  const double crot = std::cos(rot);
  const double srot = std::sin(rot);

  if(layout == nullptr)
    layout = new CameraLayout;
  else
    layout->Clear();

  layout->set_camera_type(camera_type);
  //  repeated ChannelLayout channel                           = 10 [
  layout->set_pixel_grid_layout(CameraLayout::HEX_GRID);
  layout->set_pixel_grid_spacing(spacing);
  layout->set_pixel_grid_rotation(rot/M_PI*180.0);
  layout->set_pixel_grid_cos_rotation(crot);
  layout->set_pixel_grid_sin_rotation(srot);
  layout->set_pixel_grid_offset_x(0);
  layout->set_pixel_grid_offset_y(0);
  layout->set_pixel_grid_geometric_area(
    SQR(layout->pixel_grid_spacing())*sqrt(3)/2);
  layout->set_channels_per_modules(7);
  //  ModuleLayout module                                      = 21 [
  //    (CFO).desc = "Layout of each module." ];
  layout->set_adc_gains(CameraLayout::PARALLEL_DUAL_GAIN);
  layout->set_adc_bits(12);
  layout->set_can_read_waveforms(true);
  layout->set_can_read_charges(true);
  layout->set_can_read_peak_sample(false);

  auto& mod_ctr_id_to_name = *layout->mutable_module_counter_id_to_name();
  mod_ctr_id_to_name[0] = "global_event_counter";
  mod_ctr_id_to_name[1] = "bunch_counter";
  mod_ctr_id_to_name[2] = "event_counter";
  mod_ctr_id_to_name[3] = "ts1";
  mod_ctr_id_to_name[4] = "ts2_bunch";
  mod_ctr_id_to_name[5] = "ts2_event";
  mod_ctr_id_to_name[6] = "ts2_empty";
  for(auto&& i : layout->module_counter_id_to_name()) {
    (*layout->mutable_module_counter_name_to_id())[i.second] = i.first;
  }

#if 0 // Obsolete as CDTS has its own structure type
  auto& cam_ctr_id_to_name = *layout->mutable_camera_counter_id_to_name();
  cam_ctr_id_to_name[0] = "cdts_event_counter";
  cam_ctr_id_to_name[1] = "cdts_pps_counter";
  cam_ctr_id_to_name[2] = "cdts_clock_counter";
  cam_ctr_id_to_name[3] = "cdts_ucts_timestamp_ns";
  cam_ctr_id_to_name[4] = "cdts_camera_timestamp_ns";
  cam_ctr_id_to_name[5] = "cdts_trigger_type";
  cam_ctr_id_to_name[6] = "cdts_white_rabbit_status";
  cam_ctr_id_to_name[7] = "cdts_arbitrary_information";
  for(auto&& i : layout->camera_counter_id_to_name()) {
    (*layout->mutable_camera_counter_name_to_id())[i.second] = i.first;
  }
#endif

  //const unsigned modchanmap[] = { 3, 0, 6, 5, 4, 2, 1 };
  const unsigned gridchanmap[] = { 1, 6, 5, 0, 4, 3, 2 };

  std::map<unsigned, unsigned> grid_chan_index;
  for(unsigned imod = 0; imod<modvec.size(); imod++)
  {
    const auto& mod = modvec[imod];
    auto mhid = cluster_hexid_to_member_hexid(mod.imod, 1, false);
    for(unsigned imodchan = 0; imodchan<mhid.size(); imodchan++)
    {
      unsigned ichan = imod*7+imodchan;
      unsigned igridchan = mhid[gridchanmap[imodchan]];
      grid_chan_index[igridchan] = ichan;
    }
  }

  for(unsigned imod = 0; imod<modvec.size(); imod++)
  {
    const auto& mod = modvec[imod];
    auto* m = layout->add_module();
    auto mhid = cluster_hexid_to_member_hexid(mod.imod, 1, false);
    for(unsigned imodchan = 0; imodchan<mhid.size(); imodchan++)
    {
      unsigned ichan = imod*7+imodchan;
      unsigned igridchan = mhid[gridchanmap[imodchan]];
      int u;
      int v;
      hexid_to_uv(igridchan, u, v);
      double x;
      double y;
      uv_to_xy_trans(u,v,x,y,crot,srot,spacing);
      m->add_channels_in_module(ichan);
      auto* c = layout->add_channel();
      c->set_channel_index(ichan);
      c->set_pixel_index(ichan);
      c->set_pixel_grid_index(igridchan);
      c->set_channel_set_index(0);
      c->set_module_index(imod);
      c->set_module_channel_index(imodchan);
      c->set_x(x);
      c->set_y(y);
      c->set_diameter(layout->pixel_grid_spacing());
      c->set_geometric_area(layout->pixel_grid_geometric_area());

      for(auto nigrid : hexid_to_neighbor_hexids(igridchan))
      {
        auto nid = grid_chan_index.find(nigrid);
        if(nid != grid_chan_index.end())
          c->add_neighbour_channel_indexes(nid->second);
      }

      std::vector<double> vertex_x;
      std::vector<double> vertex_y;
      uv_to_vertexes_xy_trans(u, v, vertex_x, vertex_y, crot, srot, spacing);
      for(auto vx : vertex_x)c->add_pixel_polygon_vertex_x(vx);
      for(auto vy : vertex_y)c->add_pixel_polygon_vertex_y(vy);
    }
  }

  return layout;
}

} // anonymous namespace

CameraLayout*
calin::iact_data::nectarcam_layout::nectarcam_19module_layout(
  CameraLayout* layout)
{
  const double spacing = 5;
  double rot = 0;
  auto modvec = nectarcam_mod_map(2, 0.0, spacing, rot);
  return nectarcam_general_layout(layout, modvec,
    CameraLayout::NECTARCAM_TESTBENCH_19CHANNEL, spacing, rot);
}

CameraLayout*
calin::iact_data::nectarcam_layout::nectarcam_layout(
  CameraLayout* layout)
{
  const double spacing = 5;
  double rot = 0;
  auto modvec = nectarcam_mod_map(9, 23.5*spacing, spacing, rot);
  return nectarcam_general_layout(layout, modvec,
    CameraLayout::NECTARCAM, spacing, rot);
}
