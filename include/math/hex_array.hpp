/* 

   calin/io/hex_array.hpp -- Stephen Fegan -- 2015-10-21

   Collection of functions which translate between hexagonal and Cartesian
   geometries, and provide other useful calculations for hex grids.

*/

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace calin { namespace math { namespace hex_array {

namespace vvv {
// Reference xy <-> hex functions (Vladimir Vassiliev, March 24, 2000)
void xy_to_nh(double *, double *, int *);
void nh_to_xy(int *, double *, double *);
} // namespace vvv

// *****************************************************************************
//
// General functions for the hex array 
//
// *****************************************************************************

inline unsigned positive_hexid_to_ringid_root(unsigned hexid)
{
  const unsigned iarg = 1+4*(hexid-1)/3;
  //return static_cast<unsigned>(floor((std::sqrt(iarg)-1)/2)+1);
  //return (static_cast<unsigned>(floor(std::sqrt(iarg)))-1)/2+1;
  return (unsigned(std::sqrt(iarg))-1)/2+1;
}

inline unsigned positive_hexid_to_ringid_loop(unsigned hexid)
{
  unsigned ringid = 0;
  unsigned nsites = 1;
  unsigned dn = 6;
  while(nsites<=hexid) { ++ringid; nsites+=dn; dn+=6; }
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

constexpr inline unsigned ringid_to_nsites(unsigned ringid)
{
  return 3*ringid*(ringid+1)+1;
}

inline void
positive_hexid_to_ringid_segid_runid(unsigned hexid, unsigned& ringid,
                                     unsigned& segid, unsigned& runid)
{
  ringid = positive_hexid_to_ringid(hexid);
  unsigned iring = hexid - ringid_to_nsites(ringid-1);
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
  return ringid_to_nsites(ringid-1)+segid*ringid+runid;
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

// *****************************************************************************
//
// XY <-> UV
//
// *****************************************************************************

constexpr double c_ux = 1.0;
constexpr double c_uy = 0.0;
constexpr double c_vx = 0.5;
constexpr double c_vy = 0.5*1.73205080756887729352744634150587;

inline void uv_to_xy(int u, int v, double& x, double& y) 
{
  x = u*c_ux + v*c_vx;
  y = u*c_uy + v*c_vy;
}

inline void xv_to_uv(double x, double y, int& u, int& v)
{
  
}

void uv_to_xy(int u, int v, double& x, double& y,
              bool clockwise);
void xv_to_uv(double x, double y, int& u, int& v,
              bool clockwise=false);




unsigned uv_to_hexid(int u, int v);
void hexid_to_uv(unsigned hexid, int& u, int& v);



} } } // namespace calin::math::hex_array
