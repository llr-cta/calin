/* 

   calin/io/hex_array.hpp -- Stephen Fegan -- 2015-10-21

   Collection of functions which translate between hexagonal and Cartesian
   geometries, and provide other useful calculations for hex grids.

   Some portions of this file are Copyright 2000, Vladimir Vassiliev

   The remaining portions are:

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

namespace calin { namespace math { namespace hex_array {

// *****************************************************************************
//
// Reference xy <-> hex functions from Vladimir Vassiliev (March 24, 2000)
//
// *****************************************************************************

#ifndef SWIG
namespace vvv {
void xy_to_nh(double *, double *, int *);
void nh_to_xy(int *, double *, double *);
} // namespace vvv
#endif

// *****************************************************************************
//
// General functions for the hex array building to UV <-> HexID and neighbors
//
// *****************************************************************************

inline unsigned positive_hexid_to_ringid_root(unsigned hexid)
{
  // See: http://www.codecodex.com/wiki/Calculate_an_integer_square_root
  const unsigned iarg = 1+4*(hexid-1)/3;
  // return static_cast<unsigned>(floor((std::sqrt(iarg)-1)/2)+1);
  // return (static_cast<unsigned>(floor(std::sqrt(iarg)))-1)/2+1;
  // return (unsigned(std::sqrt(double(iarg)))-1)/2+1;
  // sqrt(float) is quicker here but loses precision by iarg ~= 2^24 -
  // you have been warned :-)
  return (unsigned(std::sqrt(float(iarg)))-1)/2+1;
}

inline unsigned positive_hexid_to_ringid_loop(unsigned hexid)
{
  unsigned ringid = 0;
  unsigned nsites = 1;
  while(nsites<=hexid) { ++ringid; nsites+=6*ringid; }
  return ringid;
}

inline unsigned positive_hexid_to_ringid(unsigned hexid)
{
  return positive_hexid_to_ringid_root(hexid);
}

inline unsigned hexid_to_ringid(unsigned hexid)
{
  if(hexid==0)return 0;
  return positive_hexid_to_ringid(hexid);
}

constexpr inline unsigned ringid_to_nsites_contained(unsigned ringid)
{
  return 3*ringid*(ringid+1)+1;
}

inline void
positive_hexid_to_ringid_segid_runid(unsigned hexid, unsigned& ringid,
                                     unsigned& segid, unsigned& runid)
{
  ringid = positive_hexid_to_ringid(hexid);
  unsigned iring = hexid - ringid_to_nsites_contained(ringid-1);
  segid = int(iring/ringid);
  runid = iring - segid*ringid;
}

inline void
hexid_to_ringid_segid_runid(unsigned hexid, unsigned& ringid,
                            unsigned& segid, unsigned& runid)
{
  if(hexid==0) { ringid = segid = runid = 0; return; }
  return positive_hexid_to_ringid_segid_runid(hexid, ringid, segid, runid);
}

constexpr inline unsigned
positive_ringid_segid_runid_to_hexid(unsigned ringid, unsigned segid,
                                     unsigned runid)
{
  return ringid_to_nsites_contained(ringid-1)+segid*ringid+runid;
}

constexpr inline unsigned
ringid_segid_runid_to_hexid(unsigned ringid, unsigned segid,
                            unsigned runid)
{
  return (ringid==0) ? 0 :
      positive_ringid_segid_runid_to_hexid(ringid, segid, runid);
}

inline unsigned uv_to_ringid(int u, int v)
{
  return static_cast<unsigned>(std::max({std::abs(u), std::abs(v),
            std::abs(u+v)}));
}

void hexid_to_uv(unsigned hexid, int& u, int& v);
unsigned uv_to_hexid(int u, int v);

void uv_to_neighbor_uv(int u, int v, std::vector<int>& u_neighbors,
                       std::vector<int>& v_neighbors);
std::vector<unsigned> hexid_to_neighbor_hexids(unsigned hexid);

// *****************************************************************************
//
// Cluster
//
// *****************************************************************************

inline void cluster_uv_to_center_uv(int cluster_u, int cluster_v,
                                    unsigned cluster_nring, int& u, int& v)
{
  u = (1+cluster_nring)*cluster_u - cluster_nring*cluster_v;
  v = (1+2*cluster_nring)*cluster_v + cluster_nring*cluster_u;
}

inline void cluster_hexid_to_center_uv(unsigned cluster_hexid,
                                       unsigned cluster_nring, int& u, int& v)
{
  int cluster_u;
  int cluster_v;
  hexid_to_uv(cluster_hexid, cluster_u, cluster_v);
  cluster_uv_to_center_uv(cluster_u, cluster_v, cluster_nring, u, v);
}

inline unsigned cluster_hexid_to_center_hexid(unsigned cluster_hexid,
                                              unsigned cluster_nring)
{
  int u;
  int v;
  cluster_hexid_to_center_uv(cluster_hexid, cluster_nring, u, v);
  return uv_to_hexid(u, v);
}

void cluster_hexid_to_member_uv(unsigned cluster_hexid, unsigned cluster_nring,
                                std::vector<int>& u, std::vector<int>& v);

std::vector<unsigned> cluster_hexid_to_member_hexid(unsigned cluster_hexid,
                                                    unsigned cluster_nring);

// *****************************************************************************
//
// XY <-> UV
//
// *****************************************************************************

#define CALIN_HEX_ARRAY_SQRT3 1.73205080756887729352744634150587

constexpr double c_vx = 0.5;
constexpr double c_vy = 0.5*CALIN_HEX_ARRAY_SQRT3;
constexpr double c_vy_inv = 1.0/c_vy;

inline void uv_to_xy(int u, int v, double& x, double& y) 
{
  x = u + v*c_vx;
  y = v*c_vy;
}

void xy_to_uv(double x, double y, int& u, int& v);

inline void xy_to_uv_with_remainder(double& x_in_dx_out,
                                    double& y_in_dy_out, int& u, int& v)
{
  xy_to_uv(x_in_dx_out,  y_in_dy_out, u, v);
  x_in_dx_out -= u + v*c_vx;
  y_in_dy_out -= v*c_vy;
}

inline void uv_to_xy(int u, int v, double& x, double& y, bool clockwise)
{
  uv_to_xy(u,v,x,y);
  if(clockwise)y=-y;
}

inline void xy_to_uv(double x, double y, int& u, int& v, bool clockwise)
{
  if(clockwise)xy_to_uv(x, -y, u, v);
  else xy_to_uv(x, y, u, v);
}

inline void xy_to_uv_with_remainder(double& x_in_dx_out, double& y_in_dy_out,
                                    int& u, int& v, bool clockwise)
{
  if(clockwise)
  {
    double yy = -y_in_dy_out;
    xy_to_uv_with_remainder(x_in_dx_out, yy, u, v);
    y_in_dy_out = -yy;
  }
  else xy_to_uv_with_remainder(x_in_dx_out, y_in_dy_out, u, v);
}

inline unsigned xy_to_hexid(double x, double y)
{
  int u;
  int v;
  xy_to_uv(x, y, u, v);
  return uv_to_hexid(u,v);
}

inline unsigned xy_to_hexid(double x, double y, bool clockwise)
{
  int u;
  int v;
  xy_to_uv(x, y, u, v, clockwise);
  return uv_to_hexid(u,v);
}

inline unsigned xy_to_hexid_with_remainder(double& x_in_dx_out,
                                           double& y_in_dy_out)
{
  int u;
  int v;
  xy_to_uv_with_remainder(x_in_dx_out, y_in_dy_out, u, v);
  return uv_to_hexid(u,v);
}

inline unsigned xy_to_hexid_with_remainder(double& x_in_dx_out,
                                           double& y_in_dy_out, bool clockwise)
{
  int u;
  int v;
  xy_to_uv_with_remainder(x_in_dx_out, y_in_dy_out, u, v, clockwise);
  return uv_to_hexid(u,v);
}

inline void hexid_to_xy(unsigned hexid, double& x, double& y)
{
  int u;
  int v;
  hexid_to_uv(hexid, u, v);
  uv_to_xy(u,v,x,y);
}

inline void hexid_to_xy(unsigned hexid, double& x, double& y, bool clockwise)
{
  int u;
  int v;
  hexid_to_uv(hexid, u, v);
  uv_to_xy(u,v,x,y,clockwise);
}

inline void uv_to_vertexes_xy(int u, int v,
                              std::vector<double>& x, std::vector<double>& y,
                              bool clockwise = false)
{
  static constexpr double dx = 0.5;
  static constexpr double dy = 0.5/CALIN_HEX_ARRAY_SQRT3;
  double xc;
  double yc;
  uv_to_xy(u,v,xc,yc,clockwise);
  x = { xc+dx, xc, xc-dx, xc-dx, xc, xc+dx };
  y = { yc+dy, yc+2*dy, yc+dy, yc-dy, yc-2*dy, yc-dy };
}

inline void hexid_to_vertexes_xy(unsigned hexid,
                                 std::vector<double>& x, std::vector<double>& y,
                                 bool clockwise = false)
{
  int u;
  int v;
  hexid_to_uv(hexid, u, v);
  uv_to_vertexes_xy(u, v, x, y, clockwise);
}

inline void cluster_hexid_to_center_xy(unsigned cluster_hexid,
                                       unsigned cluster_nring,
                                       double& x, double& y)
{
  int u;
  int v;
  cluster_hexid_to_center_uv(cluster_hexid, cluster_nring, u, v);
  uv_to_xy(u, v, x, y);
}

} } } // namespace calin::math::hex_array
