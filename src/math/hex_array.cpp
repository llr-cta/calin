/* 

   calin/io/hex_array.cpp -- Stephen Fegan -- 2015-10-21

   Collection of functions which translate between hexagonal and Cartesian
   geometries, and provide other useful calculations for hex grids.

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
  unsigned ringid = uv_to_ringid(u,v);
  unsigned segid;
  unsigned runid;
  if(u+v==ringid)       { segid=0; runid=v; }
  else if(v==ringid)    { segid=1; runid=-u; }
  else if(u==-ringid)   { segid=2; runid=ringid-v; }
  else if(u+v==-ringid) { segid=3; runid=-v; }
  else if(v==-ringid)   { segid=4; runid=u; }
  else if(u==ringid)    { segid=5; runid=ringid+v; }
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

