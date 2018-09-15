/*

   calin/math/geometry_vcl.hpp -- Stephen Fegan -- 2018-08-26

   Misc geometrical functions for VCL :
   - Box intersections
   - Rotation matrix for angles

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

#include <limits>

#include <util/vcl.hpp>
#include <math/rng_vcl.hpp>

namespace calin { namespace math { namespace geometry {

template<typename VCLReal> struct VCL: public VCLReal
{
public:
  using typename VCLReal::real_t;
  using typename VCLReal::real_vt;
  using typename VCLReal::bool_vt;
  using typename VCLReal::vec3_vt;
  using typename VCLReal::mat3_vt;

  static inline bool_vt box_has_future_intersection_dirinv(real_vt& tmin, real_vt& tmax,
    const vec3_vt& min_corner, const vec3_vt& max_corner,
    const vec3_vt& pos,
    const real_vt& ux_inv, const real_vt& uy_inv, const real_vt& uz_inv)
  {
    const vec3_vt min_rel = min_corner - pos;
    const vec3_vt max_rel = max_corner - pos;

    const real_vt tx1 = min_rel.x() * ux_inv;
    const real_vt tx2 = max_rel.x() * ux_inv;
    tmin = min(tx1, tx2);
    tmax = max(tx1, tx2);

    const real_vt ty1 = min_rel.y() * uy_inv;
    const real_vt ty2 = max_rel.y() * uy_inv;
    tmin = max(tmin, min(min(ty1, ty2), tmax));
    tmax = min(tmax, max(max(ty1, ty2), tmin));

    const real_vt tz1 = min_rel.z() * uz_inv;
    const real_vt tz2 = max_rel.z() * uz_inv;
    tmin = max(tmin, min(min(tz1, tz2), tmax));
    tmax = min(tmax, max(max(tz1, tz2), tmin));

    return tmax > max(tmin, 0.0);
  }

  static inline bool_vt box_has_future_intersection(real_vt& tmin, real_vt& tmax,
    const vec3_vt& min_corner, const vec3_vt& max_corner,
    const vec3_vt& pos, const vec3_vt& dir)
  {
    return box_has_future_intersection_dirinv(tmin, tmax,
      min_corner, max_corner, pos, 1.0f/dir.x(), 1.0f/dir.y(), 1.0f/dir.z());
  }

  static inline void rotation_z_to_xyz_Rzy(mat3_vt& m,
    const real_vt x, const real_vt y, const real_vt z)
  {
    real_vt st = sqrt(x*x+y*y);
    real_vt st_inv = 1.0/st;

    // If we are aong the z-axis (z=+/-1) then we arrange for matrix to be
    // diagonal with (+1,1,+1) or (-1,1,-1)
    bool_vt is_uz = st < std::numeric_limits<real_t>::epsilon();

    real_vt sp = select(is_uz, 0.0, y * st_inv);
    real_vt cp = select(is_uz, 1.0, x * st_inv);

    m << z*cp, -sp, x,
         z*sp,  cp, y,
          -st,   0, z;
  }

  static inline void rotation_z_to_xyz_Rzyz(mat3_vt& m,
    const real_vt x, const real_vt y, const real_vt z)
  {
    real_vt st = sqrt(x*x+y*y);
    real_vt st_inv = 1.0/st;

    // If we are aong the z-axis (z=+/-1) then we arrange for matrix to be
    // diagonal with (+1,1,+1) or (-1,1,-1)
    bool_vt is_uz = st < std::numeric_limits<real_t>::epsilon();

    real_vt sp = select(is_uz, 0.0, y * st_inv);
    real_vt cp = select(is_uz, 1.0, x * st_inv);
    real_vt sp2 = sp*sp;
    real_vt cp2 = cp*cp;
    real_vt zmo_cpsp = (z-1)*cp*sp;

    m << z*cp2+sp2,  zmo_cpsp, x,
          zmo_cpsp, z*sp2+cp2, y,
                -x,        -y, z;
  }

  static inline void rotation_y_to_xyz_Ryxy(mat3_vt& m,
    const real_vt x, const real_vt y, const real_vt z)
  {
    real_vt st = sqrt(z*z+x*x);
    real_vt st_inv = 1.0/st;

    // If we are aong the z-axis (z=+/-1) then we arrange for matrix to be
    // diagonal with (+1,1,+1) or (-1,1,-1)
    bool_vt is_uy = st < std::numeric_limits<real_t>::epsilon();

    real_vt sp = select(is_uy, 0.0, x * st_inv);
    real_vt cp = select(is_uy, 1.0, z * st_inv);
    real_vt sp2 = sp*sp;
    real_vt cp2 = cp*cp;
    real_vt ymo_cpsp = (y-1)*cp*sp;

    m << y*sp2+cp2, x,  ymo_cpsp,
                -x, y,        -z,
          ymo_cpsp, z, y*cp2+sp2;
  }

  static inline void
  scatter_direction(vec3_vt& v, real_vt dispersion,
    calin::math::rng::VCLRealRNG<VCLReal>& rng)
  {
    vec3_vt x;
    rng.normal_two_bm(x.x(), x.y());
    x.x() *= dispersion;
    x.y() *= dispersion;
    // x.z() = sqrt(1.0-x.x()*x.x()-x.y()*x.y());
    x.z() = sqrt(nmul_add(x.y(),x.y(),nmul_add(x.x(),x.x(),1.0)));

    mat3_vt m;
    rotation_z_to_xyz(m, v.x(), v.y(), v.z());

    return m*x;
  }

};

#if 0
inline bool box_has_future_intersection(
  const Eigen::Vector3d& min_corner, const Eigen::Vector3d& max_corner,
  const Eigen::Vector3d& pos, const Eigen::Vector3d& dir)
{
  double tmin;
  double tmax;
  return box_has_future_intersection(tmin,tmax,min_corner,max_corner,pos,dir);
}

inline void rotation_theta_phi(Eigen::Matrix3d& m,
  const double ct, const double st, const double cp, const double sp)
{
  m << ct*cp, -sp, st*cp,
       ct*sp,  cp, st*sp,
         -sp,   0,    ct;
}

inline void rotation_theta_phi(Eigen::Matrix3d& m, double theta, double phi)
{
  const double ct = std::cos(theta);
  const double st = std::sin(theta);
  const double cp = std::cos(phi);
  const double sp = std::sin(phi);
  return rotation_theta_phi(m, ct, st, cp, sp);
}

inline void rotation_z_to_xyz(Eigen::Matrix3d& m,
  const double x, const double y, const double z)
{
  double st = std::sqrt(x*x+y*y);
  if(st == 0.0) {
    if(z>=0) {
      m.setIdentity();
    } else {
      m << -1,  0,  0,
            0,  1,  0,
            0,  0, -1;
    }
  } else {
    double sp = y/st;
    double cp = x/st;
    m << z*cp, -sp, x,
         z*sp,  cp, y,
          -st,   0, z;
  }
}

inline void rotation_z_to_vec(Eigen::Matrix3d& m, const Eigen::Vector3d& v)
{
  return rotation_z_to_xyz(m, v.x(), v.y(), v.z());
}

#ifndef SWIG

// -----------------------------------------------------------------------------
// Skip in SWIG as templates map above functions in equivalents of these
// -----------------------------------------------------------------------------

inline Eigen::Matrix3d rotation_theta_phi(
  const double ct, const double st, const double cp, const double sp)
{
  Eigen::Matrix3d m;
  rotation_theta_phi(m, ct, st, cp, sp);
  return m;
}

inline Eigen::Matrix3d rotation_theta_phi(double theta, double phi)
{
  Eigen::Matrix3d m;
  rotation_theta_phi(m, theta, phi);
  return m;
}

inline Eigen::Matrix3d rotation_z_to_xyz(
  const double x, const double y, const double z)
{
  Eigen::Matrix3d m;
  rotation_z_to_xyz(m, x, y, z);
  return m;
}

inline Eigen::Matrix3d rotation_z_to_vec(const Eigen::Vector3d v)
{
  Eigen::Matrix3d m;
  rotation_z_to_vec(m, v);
  return m;
}
#endif

#endif

} } } // namespace calin::math::geometry
