/*

   calin/math/hex_array_vcl.hpp -- Stephen Fegan -- 2018-08-20

   Collection of functions which translate between hexagonal and Cartesian
   geometries, and provide other useful calculations for hex grids.

   VCL version.

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

#pragma once

#include <math/simd.hpp>
#include <math/hex_array.hpp>

namespace calin { namespace math { namespace hex_array {

template<typename VCLArchitecture> class VCLHexArray: public VCLArchitecture
{
public:
  using typename VCLArchitecture::uint32_vt;
  using typename VCLArchitecture::int32_vt;
  using typename VCLArchitecture::uint64_vt;
  using typename VCLArchitecture::int64_vt;
  using typename VCLArchitecture::float_vt;
  using typename VCLArchitecture::double_vt;

  using typename VCLArchitecture::uint32_bvt;
  using typename VCLArchitecture::int32_bvt;
  using typename VCLArchitecture::uint64_bvt;
  using typename VCLArchitecture::int64_bvt;
  using typename VCLArchitecture::float_bvt;
  using typename VCLArchitecture::double_bvt;

  // ***************************************************************************
  //
  // General functions for the hex array building to UV <-> HexID and neighbors
  //
  // ***************************************************************************

  static inline int32_vt positive_hexid_to_ringid_loop(const int32_vt hexid)
  {
    // This algorithm is relatively slow in comparisson to the scalar version
    // but still faster overall conidering we compute 8 rigids in one go
    int32_vt ringid = 0;
    int32_vt nsites = 1;
    int32_vt nring = 0;
    int32_bvt mask = nsites > hexid;
    while(!horizontal_and(mask)) {
      ringid = select(mask, ringid, ringid+1);
      nring += 6;
      nsites += nring;
      mask = nsites > hexid;
    }
    return ringid;
  }

  static inline int32_vt positive_hexid_to_ringid_root(const int32_vt hexid)
  {
    // The following algorithm works until hexid=12,589,056
    // const unsigned iarg = 1+4*(hexid-1)/3;
    // return (unsigned(std::sqrt(float(iarg)))+1)/2;
    float_vt arg = to_float(hexid);
    arg = mul_sub(arg, 4.0/3.0, 1.0/3.0);
    arg = sqrt(arg);
    arg = mul_add(arg, 0.5, 0.5);
    return truncate_to_int(arg);
  }

  static inline int32_vt positive_hexid_to_ringid(const int32_vt hexid)
  {
    return positive_hexid_to_ringid_root(hexid);
  }

  static inline int32_vt hexid_to_ringid(const int32_vt hexid)
  {
    return select(hexid>0, positive_hexid_to_ringid(hexid), 0);
  }

  static inline int32_vt ringid_to_nsites_contained(const int32_vt ringid)
  {
    return 3*ringid*(ringid+1)+1;
  }

  static inline void positive_hexid_to_ringid_segid_runid(
    const int32_vt hexid, int32_vt& ringid, int32_vt& segid, int32_vt& runid)
  {
    // ringid = positive_hexid_to_ringid(hexid);
    // unsigned iring = hexid - ringid_to_nsites_contained(ringid-1);
    // segid = int(iring/ringid);
    // runid = iring - segid*ringid;

    ringid = positive_hexid_to_ringid(hexid);
    const int32_vt ringid_minus_one = ringid - 1;
    runid = hexid - ringid_to_nsites_contained(ringid_minus_one);
    segid = 0;

    int32_bvt mask = runid > ringid_minus_one;
    runid = if_add(mask, runid, -ringid);
    segid = if_add(mask, segid, 1);

    mask = runid > ringid_minus_one;
    runid = if_add(mask, runid, -ringid);
    segid = if_add(mask, segid, 1);

    mask = runid > ringid_minus_one;
    runid = if_add(mask, runid, -ringid);
    segid = if_add(mask, segid, 1);

    mask = runid > ringid_minus_one;
    runid = if_add(mask, runid, -ringid);
    segid = if_add(mask, segid, 1);

    mask = runid > ringid_minus_one;
    runid = if_add(mask, runid, -ringid);
    segid = if_add(mask, segid, 1);
  }

  static inline void hexid_to_ringid_segid_runid(
    const int32_vt hexid, int32_vt& ringid, int32_vt& segid, int32_vt& runid)
  {
    // if(hexid==0) { ringid = segid = runid = 0; return; }
    // return positive_hexid_to_ringid_segid_runid(hexid, ringid, segid, runid);
    positive_hexid_to_ringid_segid_runid(hexid, ringid, segid, runid);
    int32_bvt mask = hexid > 0;
    ringid = select(mask, ringid, 0);
    segid = select(mask, segid, 0);
    runid = select(mask, runid, 0);
  }

  static inline int32_vt positive_ringid_segid_runid_to_hexid(
    const int32_vt ringid, const int32_vt segid, const int32_vt runid)
  {
    return ringid_to_nsites_contained(ringid-1)+segid*ringid+runid;
  }

  static inline int32_vt ringid_segid_runid_to_hexid(
    const int32_vt ringid, const int32_vt segid, const int32_vt runid)
  {
    // return (ringid==0) ? 0 :
    //     positive_ringid_segid_runid_to_hexid(ringid, segid, runid);
    return select(ringid > 0, positive_ringid_segid_runid_to_hexid(ringid, segid, runid), 0);
  }

  static inline int32_vt uv_to_ringid(const int32_vt u, const int32_vt v)
  {
    // return static_cast<unsigned>(std::max({std::abs(u), std::abs(v),
    //           std::abs(u+v)}));
    int32_vt ringid = abs(u);
    ringid = max(ringid, abs(v));
    ringid = max(ringid, abs(u+v));
    return ringid;
  }

  static inline void hexid_to_uv_ccw(const int32_vt hexid, int32_vt& u, int32_vt& v)
  {
    // if(hexid==0) { u = v = 0; return; }
    // unsigned ringid;
    // unsigned segid;
    // unsigned runid;
    // positive_hexid_to_ringid_segid_runid(hexid, ringid, segid, runid);
    // switch(segid)
    // {
    //   case 0: u = ringid-runid; v = runid;        break;
    //   case 1: u = -runid;       v = ringid;       break;
    //   case 2: u = -ringid;      v = ringid-runid; break;
    //   case 3: u = runid-ringid; v = -runid;       break;
    //   case 4: u = runid;        v = -ringid;      break;
    //   case 5: u = ringid;       v = runid-ringid; break;
    //   default: assert(0);
    // }
    int32_vt ringid = positive_hexid_to_ringid(hexid);
    int32_vt iring = hexid - ringid_to_nsites_contained(ringid - 1);

    u = ringid;
    v = 0;

    int32_vt irun = min(iring, ringid);
    u -= irun;
    v += irun;
    iring -= irun;

    irun = min(iring, ringid);
    u -= irun;
    iring -= irun;

    irun = min(iring, ringid);
    v -= irun;
    iring -= irun;

    irun = min(iring, ringid);
    u += irun;
    v -= irun;
    iring -= irun;

    irun = min(iring, ringid);
    u += irun;
    iring -= irun;

    v += iring;

    int32_bvt mask = hexid > 0;
    u = select(mask, u, 0);
    v = select(mask, v, 0);
  }

  static inline void hexid_to_uv_cw(const int32_vt hexid, int32_vt& u, int32_vt& v)
  {
    hexid_to_uv_ccw(hexid, u, v);
    u += v;
    v = -v;
  }

  static inline int32_vt uv_to_hexid_ccw(const int32_vt u, const int32_vt v)
  {
    // if(u==0 and v==0)return 0;
    // int ringid = uv_to_ringid(u,v);
    // unsigned segid;
    // int runid;
    // int upv = u+v;
    // if(upv==ringid and v!=ringid)         { segid=0; runid=v; }
    // else if(v==ringid and u!=-ringid)     { segid=1; runid=-u; }
    // else if(u==-ringid and upv!=-ringid)  { segid=2; runid=ringid-v; }
    // else if(u+v==-ringid and v!=-ringid)  { segid=3; runid=-v; }
    // else if(v==-ringid and u!=ringid)     { segid=4; runid=u; }
    // else /*if(u==ringid and upv!=ringid)*/{ segid=5; runid=ringid+v; }
    // return positive_ringid_segid_runid_to_hexid(ringid, segid, runid);
    const int32_vt ringid = uv_to_ringid(u,v);
    const int32_vt minus_ringid = -ringid;
    const int32_vt upv = u + v;

    int32_bvt not_found_mask = true;

    int32_vt hexid = ringid_to_nsites_contained(ringid - 1);

    // Seg ID = 0
    // if(upv==ringid and v!=ringid)         { segid=0; runid=v; }
    int32_bvt here_mask = upv==ringid;
    hexid = if_add(not_found_mask, hexid, select(here_mask, v, ringid));
    not_found_mask = andnot(not_found_mask, here_mask);

    // Seg ID = 1
    // else if(v==ringid and u!=-ringid)     { segid=1; runid=-u; }
    here_mask = v==ringid;
    hexid = if_add(not_found_mask, hexid, select(here_mask, -u, ringid));
    not_found_mask = andnot(not_found_mask, here_mask);

    // Seg ID = 2
    // else if(u==-ringid and upv!=-ringid)  { segid=2; runid=ringid-v; }
    here_mask = u==minus_ringid;
    hexid = if_add(not_found_mask, hexid, select(here_mask, -upv, ringid));
    not_found_mask = andnot(not_found_mask, here_mask);

    // Seg ID = 3
    // else if(u+v==-ringid and v!=-ringid)  { segid=3; runid=-v; }
    here_mask = upv==minus_ringid;
    hexid = if_add(not_found_mask, hexid, select(here_mask, -v, ringid));
    not_found_mask = andnot(not_found_mask, here_mask);

    // Seg ID = 4
    // else if(v==-ringid and u!=ringid)     { segid=4; runid=u; }
    here_mask = v==minus_ringid;
    hexid = if_add(not_found_mask, hexid, select(here_mask, u, ringid));
    not_found_mask = andnot(not_found_mask, here_mask);

    // Seg ID = 5
    // else /*if(u==ringid and upv!=ringid)*/{ segid=5; runid=ringid+v; }
    hexid = if_add(not_found_mask, hexid, upv);

    return select(ringid>0, hexid, 0);
  }

  static inline int32_vt uv_to_hexid_cw(int32_vt u, int32_vt v)
  {
    u += v;
    v = -v;
    return uv_to_hexid_ccw(u, v);
  }

  // static inline void hexid_to_uv(int32_vt hexid, int32_vt& u, int32_vt& v)
  // {
  //   hexid_to_uv_ccw(hexid, u, v);
  // }
  //
  // static inline int32_vt uv_to_hexid(int32_vt u, int32_vt v)
  // {
  //   return uv_to_hexid_ccw(u,v);
  // }

  static inline void hexid_to_uv(int32_vt hexid, int32_vt& u, int32_vt& v, bool clockwise=false)
  {
    if(clockwise)hexid_to_uv_cw(hexid, u, v);
    else hexid_to_uv_ccw(hexid, u, v);
  }

  static inline int32_vt uv_to_hexid(int32_vt u, int32_vt v, bool clockwise=false)
  {
    if(clockwise)return uv_to_hexid_cw(u, v);
    else return uv_to_hexid_ccw(u, v);
  }

  // ***************************************************************************
  //
  // XY <-> UV
  //
  // ***************************************************************************

  // XY <-> UV without rotation matrix

  static inline void uv_to_xy_f(const int32_vt u, const int32_vt v, float_vt& x, float_vt& y)
  {
    x = to_float(u);
    y = to_float(v);
    x = mul_add(y, c_vx_f, x);
    y *= c_vy_f;
  }

  static inline void xy_to_uv_f(float_vt x, float_vt y, int32_vt& u, int32_vt& v)
  {
    // Convert X,Y first into U,V space then round to nearest
    // integer. That gets us close to correct answer, mapping XY to a
    // lozenge-shaped space rather than hexagonal. We then correct the
    // four regions that lie outside the hexagonal cell assigning them
    // to their correct neighboring cell.
    // Writer's note: see ~/Google Drive/Work/calin

    // double dv = y*c_vy_inv;
    // double du = x-dv*c_vx;
    // u = std::lround(du);
    // v = std::lround(dv);
    // du -= u;
    // dv -= v;

    y *= c_vy_inv_f;
    x = nmul_add(y, c_vx_f, x);
    u = round_to_int(x);
    v = round_to_int(y);
    x -= to_float(u);
    y -= to_float(v);

    // double c3 = dv-du;
    int32_vt c3 { reinterpret_i(y-x) }; // only sign bit needed

    int32_vt uvshift;
    int32_bvt mask;

    // double c1 = du+0.5*dv;
    // double c2 = dv+0.5*du;
    // if(c3<0) {
    //   if(c1>=1) u++;
    //   else if(c2<-1) v--;
    // } else {
    //   if(c2>=1) v++;
    //   else if(c1<-1) u--;
    // }

    uvshift = round_to_int(mul_add(y, 0.5f, x));
    mask = (uvshift ^ c3) < 0;
    u = select(mask, u+uvshift, u);

    // uvshift = _mm256_cvtps_epi32(_mm256_fmadd_ps(y, calin::math::simd::c_m256(_c_m256_one_half), x));
    // mask = _mm256_srai_epi32(_mm256_xor_si256(uvshift, c3), 31);
    // u = _mm256_blendv_epi8(u, _mm256_add_epi32(u, uvshift), mask);

    uvshift = round_to_int(mul_add(x, 0.5f, y));
    mask = (uvshift ^ c3) < 0;
    v = select(mask, v, v+uvshift);

    // uvshift = _mm256_cvtps_epi32(_mm256_fmadd_ps(x, calin::math::simd::c_m256(_c_m256_one_half), y));
    // mask = _mm256_srai_epi32(_mm256_xor_si256(uvshift, c3), 31);
    // v = _mm256_blendv_epi8(_mm256_add_epi32(v, uvshift), v, mask);
  }

  static inline void xy_to_uv_with_remainder_f(
    float_vt& x_in_dx_out, float_vt& y_in_dy_out, int32_vt& u, int32_vt& v)
  {
    xy_to_uv_f(x_in_dx_out, y_in_dy_out, u, v);
    x_in_dx_out = x_in_dx_out - to_float(u);
    float_vt vf = to_float(v);
    x_in_dx_out = nmul_add(vf, c_vx_f, x_in_dx_out);
    y_in_dy_out = nmul_add(vf, c_vy_f, y_in_dy_out);
  }

  // XY <-> UV with (shear-free) affine transformation

  static inline void uv_to_xy_trans_f(const int32_vt u, const int32_vt v, float_vt& x, float_vt& y,
    const float_vt& crot, const float_vt& srot, const float_vt& scale, const float_vt& dx = 0, const float_vt& dy = 0)
  {
    // uv_to_xy(u,v,x,y);
    // double xx = x*crot - y*srot;
    // y = scale * (y*crot + x*srot) + dy;
    // x = scale * xx + dx;
    uv_to_xy_f(u,v,x,y);
    float_vt xx = y * srot;
    xx = mul_sub(x, crot, xx);
    y = y * crot;
    y = mul_add(x, srot, y);
    x = xx;
    x = mul_add(x, scale, dx);
    y = mul_add(y, scale, dy);
  }

  static void xy_trans_to_uv_f(float_vt x, float_vt y, int32_vt& u, int32_vt& v,
    const float_vt& crot, const float_vt& srot, const float_vt& scale, const float_vt& dx = 0, const float_vt& dy = 0)
  {
    // x = (x - dx)/scale;
    // y = (y - dy)/scale;
    // double xx = x*crot + y*srot;
    // y = y*crot - x*srot;
    // xy_to_uv(xx,y,u,v);
    float_vt scale_inv = 1.0f/scale; // approx_recipr(scale);
    x = (x - dx) * scale_inv;
    y = (y - dy) * scale_inv;
    float_vt yy = x * srot;
    yy = mul_sub(y, crot, yy);
    x = x * crot;
    x = mul_add(y, srot, x);
    xy_to_uv_f(x, yy, u, v);
  }

  static void xy_trans_to_uv_scaleinv_f(float_vt x, float_vt y, int32_vt& u, int32_vt& v,
    const float_vt& crot, const float_vt& srot, const float_vt& scale_inv, const float_vt& dx = 0, const float_vt& dy = 0)
  {
    // x = (x - dx)/scale;
    // y = (y - dy)/scale;
    // double xx = x*crot + y*srot;
    // y = y*crot - x*srot;
    // xy_to_uv(xx,y,u,v);
    x = (x - dx) * scale_inv;
    y = (y - dy) * scale_inv;
    float_vt yy = x * srot;
    yy = mul_sub(y, crot, yy);
    x = x * crot;
    x = mul_add(y, srot, x);
    xy_to_uv_f(x, yy, u, v);
  }

  static inline void xy_trans_to_uv_with_remainder_f(
    float_vt& x_in_dx_out, float_vt& y_in_dy_out, int32_vt& u, int32_vt& v,
    const float_vt& crot, const float_vt& srot, const float_vt& scale, const float_vt& dx = 0, const float_vt& dy = 0)
  {
    float_vt scale_inv = 1.0f/scale; // approx_recipr(scale);
    x_in_dx_out = (x_in_dx_out - dx) * scale_inv;
    y_in_dy_out = (y_in_dy_out - dy) * scale_inv;
    float_vt yy = x_in_dx_out * srot;
    yy = mul_sub(y_in_dy_out, crot, yy);
    x_in_dx_out = x_in_dx_out * crot;
    x_in_dx_out = mul_add(y_in_dy_out, srot, x_in_dx_out);
    xy_to_uv_with_remainder_f(x_in_dx_out, yy, u, v);

    y_in_dy_out = yy * crot;
    y_in_dy_out = mul_add(x_in_dx_out, srot, y_in_dy_out);
    y_in_dy_out = y_in_dy_out * scale;

    x_in_dx_out = x_in_dx_out * crot;
    x_in_dx_out = nmul_add(yy, srot, x_in_dx_out);
    x_in_dx_out = x_in_dx_out * scale;
  }

  static inline void xy_trans_to_uv_scaleinv_with_remainder_f(
    float_vt& x_in_dx_out, float_vt& y_in_dy_out, int32_vt& u, int32_vt& v,
    const float_vt& crot, const float_vt& srot, const float_vt& scale, const float_vt& scale_inv, const float_vt& dx = 0, const float_vt& dy = 0)
  {
    x_in_dx_out = (x_in_dx_out - dx) * scale_inv;
    y_in_dy_out = (y_in_dy_out - dy) * scale_inv;
    float_vt yy = x_in_dx_out * srot;
    yy = mul_sub(y_in_dy_out, crot, yy);
    x_in_dx_out = x_in_dx_out * crot;
    x_in_dx_out = mul_add(y_in_dy_out, srot, x_in_dx_out);
    xy_to_uv_with_remainder_f(x_in_dx_out, yy, u, v);

    y_in_dy_out = yy * crot;
    y_in_dy_out = mul_add(x_in_dx_out, srot, y_in_dy_out);
    y_in_dy_out = y_in_dy_out * scale;

    x_in_dx_out = x_in_dx_out * crot;
    x_in_dx_out = nmul_add(yy, srot, x_in_dx_out);
    x_in_dx_out = x_in_dx_out * scale;
  }

  // XY <-> HEXID without rotation matrix

  static inline void hexid_to_xy_ccw_f(int32_vt hexid, float_vt& x, float_vt& y)
  {
    int32_vt u;
    int32_vt v;
    hexid_to_uv_ccw(hexid, u, v);
    uv_to_xy_f(u,v,x,y);
  }

  static inline void hexid_to_xy_cw_f(int32_vt hexid, float_vt& x, float_vt& y)
  {
    int32_vt u;
    int32_vt v;
    hexid_to_uv_cw(hexid, u, v);
    uv_to_xy_f(u,v,x,y);
  }

  static inline void hexid_to_xy_f(int32_vt hexid, float_vt& x, float_vt& y, bool clockwise=false)
  {
    int32_vt u;
    int32_vt v;
    hexid_to_uv(hexid, u, v, clockwise);
    uv_to_xy_f(u,v,x,y);
  }

  static inline int32_vt xy_to_hexid_ccw_f(const float_vt x, const float_vt y)
  {
    int32_vt u;
    int32_vt v;
    xy_to_uv_f(x, y, u, v);
    return uv_to_hexid_ccw(u, v);
  }

  static inline int32_vt xy_to_hexid_cw_f(const float_vt x, const float_vt y)
  {
    int32_vt u;
    int32_vt v;
    xy_to_uv_f(x, y, u, v);
    return uv_to_hexid_cw(u, v);
  }

  static inline int32_vt xy_to_hexid_f(const float_vt x, const float_vt y, bool clockwise=false)
  {
    int32_vt u;
    int32_vt v;
    xy_to_uv_f(x, y, u, v);
    return uv_to_hexid(u, v, clockwise);
  }

  static inline int32_vt xy_to_hexid_with_remainder_ccw_f(float_vt& x_in_dx_out, float_vt& y_in_dy_out)
  {
    int32_vt u;
    int32_vt v;
    xy_to_uv_with_remainder_f(x_in_dx_out, y_in_dy_out, u, v);
    return uv_to_hexid_ccw(u, v);
  }

  static inline int32_vt xy_to_hexid_with_remainder_cw_f(float_vt& x_in_dx_out, float_vt& y_in_dy_out)
  {
    int32_vt u;
    int32_vt v;
    xy_to_uv_with_remainder_f(x_in_dx_out, y_in_dy_out, u, v);
    return uv_to_hexid_cw(u, v);
  }

  static inline int32_vt xy_to_hexid_with_remainder_f(float_vt& x_in_dx_out, float_vt& y_in_dy_out, bool clockwise=false)
  {
    int32_vt u;
    int32_vt v;
    xy_to_uv_with_remainder_f(x_in_dx_out, y_in_dy_out, u, v);
    return uv_to_hexid(u, v, clockwise);
  }

  // XY <-> HEXID with (shear-free) affine transformation

  static inline int32_vt xy_trans_to_hexid_ccw_f(float_vt x, float_vt y,
    const float_vt& crot, const float_vt& srot, const float_vt& scale,
    const float_vt& dx = 0, const float_vt& dy = 0)
  {
    int32_vt u;
    int32_vt v;
    xy_trans_to_uv_f(x, y, u, v, crot, srot, scale, dx, dy);
    return uv_to_hexid_ccw(u, v);
  }

  static inline int32_vt xy_trans_to_hexid_ccw_scaleinv_f(float_vt x, float_vt y,
    const float_vt& crot, const float_vt& srot, const float_vt& scale_inv,
    const float_vt& dx = 0, const float_vt& dy = 0)
  {
    int32_vt u;
    int32_vt v;
    xy_trans_to_uv_scaleinv_f(x, y, u, v, crot, srot, scale_inv, dx, dy);
    return uv_to_hexid_ccw(u, v);
  }

  static inline int32_vt xy_trans_to_hexid_cw_f(float_vt x, float_vt y,
    const float_vt& crot, const float_vt& srot, const float_vt& scale,
    const float_vt& dx = 0, const float_vt& dy = 0)
  {
    int32_vt u;
    int32_vt v;
    xy_trans_to_uv_f(x, y, u, v, crot, srot, scale, dx, dy);
    return uv_to_hexid_cw(u, v);
  }

  static inline int32_vt xy_trans_to_hexid_cw_scaleinv_f(float_vt x, float_vt y,
    const float_vt& crot, const float_vt& srot, const float_vt& scale_inv,
    const float_vt& dx = 0, const float_vt& dy = 0)
  {
    int32_vt u;
    int32_vt v;
    xy_trans_to_uv_scaleinv_f(x, y, u, v, crot, srot, scale_inv, dx, dy);
    return uv_to_hexid_cw(u, v);
  }

  static inline int32_vt xy_trans_to_hexid_f(float_vt x, float_vt y,
    const float_vt& crot, const float_vt& srot, const float_vt& scale,
    const float_vt& dx = 0, const float_vt& dy = 0,
    bool clockwise=false)
  {
    int32_vt u;
    int32_vt v;
    xy_trans_to_uv_f(x, y, u, v, crot, srot, scale, dx, dy);
    return uv_to_hexid(u, v, clockwise);
  }

  static inline int32_vt xy_trans_to_hexid_scaleinv_f(float_vt x, float_vt y,
    const float_vt& crot, const float_vt& srot, const float_vt& scale_inv,
    const float_vt& dx = 0, const float_vt& dy = 0, bool clockwise=false)
  {
    int32_vt u;
    int32_vt v;
    xy_trans_to_uv_scaleinv_f(x, y, u, v, crot, srot, scale_inv, dx, dy);
    return uv_to_hexid(u, v, clockwise);
  }

  static inline int32_vt xy_trans_to_hexid_with_remainder_ccw_f(
    float_vt& x_in_dx_out, float_vt& y_in_dy_out,
    const float_vt& crot, const float_vt& srot, const float_vt& scale,
    const float_vt& dx = 0, const float_vt& dy = 0)
  {
    int32_vt u;
    int32_vt v;
    xy_trans_to_uv_with_remainder_f(x_in_dx_out, y_in_dy_out, u, v,
      crot, srot, scale, dx, dy);
    return uv_to_hexid_ccw(u, v);
  }

  static inline int32_vt xy_trans_to_hexid_with_remainder_ccw_scaleinv_f(
    float_vt& x_in_dx_out, float_vt& y_in_dy_out,
    const float_vt& crot, const float_vt& srot, const float_vt& scale, const float_vt& scale_inv,
    const float_vt& dx = 0, const float_vt& dy = 0)
  {
    int32_vt u;
    int32_vt v;
    xy_trans_to_uv_scaleinv_with_remainder_f(x_in_dx_out, y_in_dy_out, u, v,
      crot, srot, scale, scale_inv, dx, dy);
    return uv_to_hexid_ccw(u, v);
  }

  static inline int32_vt xy_trans_to_hexid_with_remainder_cw_f(
    float_vt& x_in_dx_out, float_vt& y_in_dy_out,
    const float_vt& crot, const float_vt& srot, const float_vt& scale,
    const float_vt& dx = 0, const float_vt& dy = 0)
  {
    int32_vt u;
    int32_vt v;
    xy_trans_to_uv_with_remainder_f(x_in_dx_out, y_in_dy_out, u, v,
      crot, srot, scale, dx, dy);
    return uv_to_hexid_cw(u, v);
  }

  static inline int32_vt xy_trans_to_hexid_with_remainder_cw_scaleinv_f(
    float_vt& x_in_dx_out, float_vt& y_in_dy_out,
    const float_vt& crot, const float_vt& srot, const float_vt& scale, const float_vt& scale_inv,
    const float_vt& dx = 0, const float_vt& dy = 0)
  {
    int32_vt u;
    int32_vt v;
    xy_trans_to_uv_scaleinv_with_remainder_f(x_in_dx_out, y_in_dy_out, u, v,
      crot, srot, scale, scale_inv, dx, dy);
    return uv_to_hexid_cw(u, v);
  }

  static inline int32_vt xy_trans_to_hexid_with_remainder_f(
    float_vt& x_in_dx_out, float_vt& y_in_dy_out,
    const float_vt& crot, const float_vt& srot, const float_vt& scale,
    const float_vt& dx = 0, const float_vt& dy = 0, bool clockwise = false)
  {
    int32_vt u;
    int32_vt v;
    xy_trans_to_uv_with_remainder_f(x_in_dx_out, y_in_dy_out, u, v,
      crot, srot, scale, dx, dy);
    return uv_to_hexid(u, v, clockwise);
  }

  static inline int32_vt xy_trans_to_hexid_with_remainder_scaleinv_f(
    float_vt& x_in_dx_out, float_vt& y_in_dy_out,
    const float_vt& crot, const float_vt& srot, const float_vt& scale, const float_vt& scale_inv,
    const float_vt& dx = 0, const float_vt& dy = 0, bool clockwise = false)
  {
    int32_vt u;
    int32_vt v;
    xy_trans_to_uv_scaleinv_with_remainder_f(x_in_dx_out, y_in_dy_out, u, v,
      crot, srot, scale, scale_inv, dx, dy);
    return uv_to_hexid(u, v, clockwise);
  }

  static inline void hexid_to_xy_trans_ccw_f(int32_vt hexid, float_vt& x, float_vt& y,
    const float_vt& crot, const float_vt& srot, const float_vt& scale,
    const float_vt& dx = 0, const float_vt& dy = 0)
  {
    int32_vt u;
    int32_vt v;
    hexid_to_uv_ccw(hexid, u, v);
    uv_to_xy_trans_f(u, v, x, y, crot, srot, scale, dx, dy);
  }

  static inline void hexid_to_xy_trans_cw_f(int32_vt hexid, float_vt& x, float_vt& y,
    const float_vt& crot, const float_vt& srot, const float_vt& scale,
    const float_vt& dx = 0, const float_vt& dy = 0)
  {
    int32_vt u;
    int32_vt v;
    hexid_to_uv_cw(hexid, u, v);
    uv_to_xy_trans_f(u, v, x, y, crot, srot, scale, dx, dy);
  }

  static inline void hexid_to_xy_trans_f(int32_vt hexid, float_vt& x, float_vt& y,
    const float_vt& crot, const float_vt& srot, const float_vt& scale,
    const float_vt& dx = 0, const float_vt& dy = 0, bool clockwise = false)
  {
    int32_vt u;
    int32_vt v;
    hexid_to_uv(hexid, u, v, clockwise);
    uv_to_xy_trans_f(u, v, x, y, crot, srot, scale, dx, dy);
  }

};

#if 0

unsigned test_avx2_positive_hexid_to_ringid_loop(unsigned hexid);
unsigned test_avx2_positive_hexid_to_ringid_root(unsigned hexid);
unsigned test_avx2_hexid_to_ringid(unsigned hexid);
unsigned test_avx2_ringid_to_nsites_contained(unsigned ringid);
void test_avx2_positive_hexid_to_ringid_segid_runid(unsigned hexid,
  unsigned& ringid, unsigned& segid, unsigned& runid);
void test_avx2_hexid_to_ringid_segid_runid(unsigned hexid,
  unsigned& ringid, unsigned& segid, unsigned& runid);
unsigned test_avx2_positive_ringid_segid_runid_to_hexid(
  unsigned ringid, unsigned segid, unsigned runid);
unsigned test_avx2_ringid_segid_runid_to_hexid(
  unsigned ringid, unsigned segid, unsigned runid);
unsigned test_avx2_uv_to_ringid(int u, int v);
void test_avx2_hexid_to_uv_ccw(unsigned hexid, int& u, int& v);
void test_avx2_hexid_to_uv_cw(unsigned hexid, int& u, int& v);
unsigned test_avx2_uv_to_hexid_ccw(int u, int v);
unsigned test_avx2_uv_to_hexid_cw(int u, int v);


void test_avx2_uv_to_xy_f(int u, int v, float& x, float& y);
void test_avx2_xy_to_uv_f(float x, float y, int& u, int& v);
void test_avx2_xy_to_uv_with_remainder_f(float& x_in_dx_out, float& y_in_dy_out,
  int& u, int& v);
void test_avx2_uv_to_xy_trans_f(int u, int v, float& x, float& y,
  float crot, float srot, float scale, float dx = 0, float dy = 0);
void test_avx2_xy_trans_to_uv_f(float x, float y, int& u, int& v,
  float crot, float srot, float scale, float dx = 0, float dy = 0);
void test_avx2_xy_trans_to_uv_with_remainder_f(float& x_in_dx_out, float& y_in_dy_out,
  int& u, int& v,
  float crot, float srot, float scale, float dx = 0, float dy = 0);

unsigned test_avx2_xy_to_hexid_f(float x, float y);
unsigned test_avx2_xy_to_hexid_with_remainder_f(float& x_in_dx_out, float& y_in_dy_out);
unsigned test_avx2_xy_to_hexid_with_remainder_f(float& x_in_dx_out, float& y_in_dy_out, bool clockwise);
void test_avx2_hexid_to_xy_f(unsigned hexid, float& x, float& y);
void test_avx2_hexid_to_xy_f(unsigned hexid, float& x, float& y, bool clockwise);

unsigned test_avx2_xy_trans_to_hexid_f(float x, float y,
  float crot, float srot, float scale, float dx = 0, float dy = 0);
unsigned test_avx2_xy_trans_to_hexid_f(float x, float y, bool clockwise,
  float crot, float srot, float scale, float dx = 0, float dy = 0);
unsigned test_avx2_xy_trans_to_hexid_with_remainder_f(float& x_in_dx_out,
  float& y_in_dy_out,
  float crot, float srot, float scale, float dx = 0, float dy = 0);
unsigned test_avx2_xy_trans_to_hexid_with_remainder_f(float& x_in_dx_out,
  float& y_in_dy_out, bool clockwise,
  float crot, float srot, float scale, float dx = 0, float dy = 0);
void test_avx2_hexid_to_xy_trans_f(unsigned hexid, float& x, float& y,
  float crot, float srot, float scale, float dx = 0, float dy = 0);
void test_avx2_hexid_to_xy_trans_f(unsigned hexid, float& x, float& y, bool clockwise,
  float crot, float srot, float scale, float dx = 0, float dy = 0);

#endif

} } } // namespace calin::math::hex_array
