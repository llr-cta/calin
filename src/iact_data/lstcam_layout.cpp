/*

   calin/iact_data/lstcam_layout.cpp -- Stephen Fegan -- 2018-19-05

   Camera layout for LSTCam

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <cmath>
#include <vector>
#include <stdexcept>

#include <math/special.hpp>
#include <math/hex_array.hpp>
#include <iact_data/instrument_layout.hpp>
#include <iact_data/lstcam_layout.hpp>
#include <util/log.hpp>

using calin::math::special::SQR;
using namespace calin::math::hex_array;
using namespace calin::ix::iact_data::instrument_layout;
using namespace calin::util::log;

namespace {

struct mod_info {
  mod_info(int imod_, int uc_, int vc_, double xc_, double yc_,
      int ix_, int iy_): imod(imod_), uc(uc_), vc(vc_), xc(xc_), yc(yc_),
    ix(ix_), iy(iy_) { /* nothing to see here */ }
  int imod;
  int uc;
  int vc;
  double xc;
  double yc;
  int ix;
  int iy;
};

std::vector<mod_info>
lstcam_mod_map(unsigned nring, double radius, double spacing, double& rot)
{
  if(nring == 0 and radius <= 0)
    throw std::runtime_error("lstcam_mod_map: At least one of nring or "
      "radius must be positive.");

  int u1 = 0;
  int v1 = 0;
  cluster_hexid_to_center_uv(1,1,u1,v1,false);
  double x1 = 0;
  double y1 = 0;
  uv_to_xy(u1,v1,x1,y1);
  rot = std::atan2(-y1,x1) + /*30.0*/90/180.0*M_PI;
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
      // int ix = (vc*2+uc*3)/7;
      // double iy = -(vc+uc)+2.5*ix;
      int ix = 2*uc - vc;
      int iy = -(5*vc + 4*uc);
      modvec.emplace_back(imod, uc, vc, xc, yc, ix, iy);
    }
  }

  std::sort(modvec.begin(), modvec.end(),
    [](const mod_info& a, const mod_info& b) { if(a.ix>b.ix)return true;
      if(a.ix==b.ix)return a.iy<b.iy;
      return false; });

  return modvec;
}

CameraLayout* lstcam_general_layout(CameraLayout* layout,
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
  layout->set_can_read_charges(false);
  layout->set_can_read_peak_sample(false);

  layout->add_module_counter_name("backplane_10MHz_counter");
  layout->add_module_counter_name("pps_counter");
  layout->add_module_counter_name("event_counter");
  layout->add_module_counter_name("local_133MHz_counter");
  layout->add_module_counter_name("trigger_counter");

  // https://forge.in2p3.fr/attachments/download/59969/Pixel_Numbering_explain.png
  const unsigned gridchanmap[] = { 6, 5, 4, 0, 1, 3, 2 };

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

  std::map<int, int> pixel_hexid_to_channel;
  std::map<int, int> module_hexid_to_module;

  for(unsigned imod = 0; imod<modvec.size(); imod++)
  {
    const auto& mod = modvec[imod];
    auto* m = layout->add_module();
    m->set_module_index(imod);
    m->set_module_grid_index(mod.imod);
    module_hexid_to_module[/* rotate_ccw1_hexid_ccw(mod.imod) */ mod.imod] = imod;
    int um;
    int vm;
    hexid_to_uv(mod.imod, um, vm);
    m->set_module_grid_u(um);
    m->set_module_grid_v(vm);
    m->set_module_grid_i(vm);
    m->set_module_grid_j(2*um+vm);

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
      if(gridchanmap[imodchan] == 0) {
        m->set_x(x);
        m->set_y(y);
      }
      m->add_channels_in_module(ichan);
      auto* c = layout->add_channel();
      c->set_channel_index(ichan);

      c->set_pixel_grid_index(igridchan);
      pixel_hexid_to_channel[igridchan] = ichan;
      c->set_channel_set_index(0);
      c->set_module_index(imod);
      c->set_module_channel_index(imodchan);
      c->set_pixel_grid_u(u);
      c->set_pixel_grid_v(v);
      c->set_x(x);
      c->set_y(y);
      c->set_diameter(layout->pixel_grid_spacing());
      c->set_geometric_area(layout->pixel_grid_geometric_area());

      c->set_is_boundary_pixel(false);
      for(auto nigrid : hexid_to_neighbor_hexids(igridchan))
      {
        auto nid = grid_chan_index.find(nigrid);
        if(nid != grid_chan_index.end()) {
          c->add_neighbour_channel_indexes(nid->second);
        } else {
          c->set_is_boundary_pixel(true);
        }
      }
      if(c->is_boundary_pixel()) {
        layout->add_boundary_pixel_channel_index(ichan);
      }

      std::vector<double> vertex_x;
      std::vector<double> vertex_y;
      uv_to_vertexes_xy_trans(u, v, vertex_x, vertex_y, crot, srot, spacing);
      for(auto vx : vertex_x)c->add_outline_polygon_vertex_x(vx);
      for(auto vy : vertex_y)c->add_outline_polygon_vertex_y(vy);
      c->add_outline_polygon_vertex_index(vertex_x.size());
    }
  }

  unsigned ispiral = 0;
  for(auto igrid : module_hexid_to_module) {
    auto imod = igrid.second;
    layout->add_module_spiral_module_index(imod);
    layout->mutable_module(imod)->set_module_spiral_index(ispiral++);
  }
  ispiral = 0;
  for(auto igrid : pixel_hexid_to_channel) {
    auto ichan = igrid.second;

    layout->add_pixel_spiral_channel_index(ichan);
    layout->mutable_channel(ichan)->set_pixel_spiral_index(ispiral);

    // For LSTCAM the pixel and spiral indexes are equivalent
    layout->add_pixel_channel_index(ichan);
    layout->mutable_channel(ichan)->set_pixel_index(ispiral);

    ++ispiral;
  }

  calin::iact_data::instrument_layout::compute_camera_and_module_outlines(layout);
  return layout;
}

} // anonymous namespace

CameraLayout*
calin::iact_data::lstcam_layout::lstcam_layout(
  CameraLayout* layout)
{
  const double spacing = 5;
  double rot = 0;
  auto modvec = lstcam_mod_map(9, 23.5*spacing, spacing, rot);
  return lstcam_general_layout(layout, modvec,
    CameraLayout::LSTCAM, spacing, rot);
}
