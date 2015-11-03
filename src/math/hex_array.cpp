/* 

   calin/io/hex_array.cpp -- Stephen Fegan -- 2015-10-21

   Collection of functions which translate between hexagonal and Cartesian
   geometries, and provide other useful calculations for hex grids.

   Copyright 2015, Stephen Fegan <sfegan@gmail.com>

   This file is part of "calin"
   
   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.
    
   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <cassert>

#include "math/hex_array.hpp"

void calin::math::hex_array::
hexid_to_uv(unsigned hexid, int& u, int& v)
{
  if(hexid==0) { u = v = 0; return; }
  unsigned ringid;
  unsigned segid;
  unsigned runid;
  positive_hexid_to_ringid_segid_runid(hexid, ringid, segid, runid);
  switch(segid)
  {
    case 0: u = ringid-runid; v = runid;        break;
    case 1: u = -runid;       v = ringid;       break;
    case 2: u = -ringid;      v = ringid-runid; break;
    case 3: u = runid-ringid; v = -runid;       break;
    case 4: u = runid;        v = -ringid;      break;
    case 5: u = ringid;       v = runid-ringid; break;
    default: assert(0);
  }
}

unsigned calin::math::hex_array::
uv_to_hexid(int u, int v)
{
  if(u==0 and v==0)return 0;
  int ringid = uv_to_ringid(u,v);
  unsigned segid;
  unsigned runid;
  if(u+v==ringid)        { segid=0; runid=v; }
  else if(v==ringid)     { segid=1; runid=-u; }
  else if(u==-ringid)    { segid=2; runid=ringid-v; }
  else if(u+v==-ringid)  { segid=3; runid=-v; }
  else if(v==-ringid)    { segid=4; runid=u; }
  else /*if(u==ringid)*/ { segid=5; runid=ringid+v; }
  return positive_ringid_segid_runid_to_hexid(ringid, segid, runid);
}

void calin::math::hex_array::
uv_to_neighbor_uv(int u, int v,
                  std::vector<int>& u_neighbors, std::vector<int>& v_neighbors)
{
  u_neighbors = { u+1, u, u-1, u-1, u, u+1 };
  v_neighbors = { v, v+1, v+1, v, v-1, v-1 };
}

std::vector<unsigned> calin::math::hex_array::
hexid_to_neighbor_hexids(unsigned hexid)
{
  int u;
  int v;
  hexid_to_uv(hexid,u,v);
  return { uv_to_hexid(u+1,v), uv_to_hexid(u, v+1), uv_to_hexid(u-1, v+1),
        uv_to_hexid(u-1,v), uv_to_hexid(u, v-1), uv_to_hexid(u+1, v-1) };
}

// *****************************************************************************
//
// Cluster
//
// *****************************************************************************

void calin::math::hex_array::
cluster_hexid_to_member_uv(unsigned cluster_hexid, unsigned cluster_nring,
                           std::vector<int>& u, std::vector<int>& v)
{
  int cluster_u;
  int cluster_v;
  cluster_hexid_to_center_uv(cluster_hexid, cluster_nring,
                             cluster_u, cluster_v);
  unsigned nsites = ringid_to_nsites(cluster_nring);
  u.resize(nsites);
  v.resize(nsites);
  for(unsigned i=0;i<nsites;i++)
  {
    int ui;
    int vi;
    hexid_to_uv(i, ui, vi);
    u[i] = ui+cluster_u;
    v[i] = vi+cluster_v;
  }
}

std::vector<unsigned> calin::math::hex_array::
cluster_hexid_to_member_hexid(unsigned cluster_hexid, unsigned cluster_nring)
{
  int cluster_u;
  int cluster_v;
  cluster_hexid_to_center_uv(cluster_hexid, cluster_nring,
                             cluster_u, cluster_v);
  unsigned nsites = ringid_to_nsites(cluster_nring);
  std::vector<unsigned> hexids(nsites,0);
  for(unsigned i=0;i<nsites;i++)
  {
    int ui;
    int vi;
    hexid_to_uv(i, ui, vi);
    hexids[i] = uv_to_hexid(ui+cluster_u, vi+cluster_v);
  }
  return hexids;
}

// *****************************************************************************
//
// XY <-> UV
//
// *****************************************************************************

void calin::math::hex_array::xy_to_uv(double x, double y, int& u, int& v)
{
  // Convert X,Y first into U,V space then round to nearest
  // integer. That gets us close to correct answer, mapping XY to a
  // lozenge-shaped space rather than hexagonal. We then correct the
  // four regions that lie outside the hexagonal cell assigning them
  // to their correct neighboring cell.
  // Writer's note: see Code/Projects/CTA/Calib/Scribbles/Hex\ Test.ipynb
  double dv = y*c_vy_inv;
  double du = x-dv*c_vx;
  u = std::lround(du);
  v = std::lround(dv);
  du -= u;
  dv -= v;
  double c1 = 2*du+dv;
  double c2 = 2*dv+du;
  double c3 = dv-du;
  if(c3<0) {
    if(c1>=1) u++;
    else if(c2<-1) v--;
  } else {
    if(c2>=1) v++;
    else if(c1<-1) u--;
  }
}
