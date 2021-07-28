/*

   calin/math/geometry.hpp -- Stephen Fegan -- 2016-11-12

   Misc geometrical functions :
   - Box intersections
   - Rotation matrix for angles

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <algorithm>
#include <Eigen/Dense>

#include <common_types.pb.h>
#include <math/rng.hpp>
#include <math/least_squares.hpp>
#include <math/special.hpp>

namespace calin { namespace math { namespace geometry {

inline bool box_has_future_intersection(double& tmin, double& tmax,
  const Eigen::Vector3d& min_corner, const Eigen::Vector3d& max_corner,
  const Eigen::Vector3d& pos, const Eigen::Vector3d& dir)
{
  // See: https://tavianator.com/fast-branchless-raybounding-box-intersections/
  // and: https://tavianator.com/fast-branchless-raybounding-box-intersections-part-2-nans/

  // Normalized direction vector
  const double vx = 1.0 / dir.x();
  const double vy = 1.0 / dir.y();
  const double vz = 1.0 / dir.z();

  Eigen::Vector3d min_rel = min_corner - pos;
  Eigen::Vector3d max_rel = max_corner - pos;

  const double tx1 = min_rel.x() * vx;
  const double tx2 = max_rel.x() * vx;
  tmin = std::min(tx1, tx2);
  tmax = std::max(tx1, tx2);

  const double ty1 = min_rel.y() * vy;
  const double ty2 = max_rel.y() * vy;
  tmin = std::max(tmin, std::min(std::min(ty1, ty2), tmax));
  tmax = std::min(tmax, std::max(std::max(ty1, ty2), tmin));

  const double tz1 = min_rel.z() * vz;
  const double tz2 = max_rel.z() * vz;
  tmin = std::max(tmin, std::min(std::min(tz1, tz2), tmax));
  tmax = std::min(tmax, std::max(std::max(tz1, tz2), tmin));

  return tmax > std::max(tmin, 0.0);
}

inline bool box_has_future_intersection(
  const Eigen::Vector3d& min_corner, const Eigen::Vector3d& max_corner,
  const Eigen::Vector3d& pos, const Eigen::Vector3d& dir)
{
  double tmin;
  double tmax;
  return box_has_future_intersection(tmin,tmax,min_corner,max_corner,pos,dir);
}

inline bool oct_box_has_future_intersection(double& tmin, double& tmax,
  const Eigen::Vector3d& center, double flat_to_flat, double height,
  const Eigen::Vector3d& pos, const Eigen::Vector3d& dir)
{
  const double half_flat_to_flat = 0.5*flat_to_flat;
  const double half_height = 0.5*height;

  const double vx = 1.0 / dir.x();
  const double min_rel_x = center.x() - half_flat_to_flat - pos.x();
  const double max_rel_x = center.x() + half_flat_to_flat - pos.x();
  const double tx1 = min_rel_x * vx;
  const double tx2 = max_rel_x * vx;
  tmin = std::min(tx1, tx2);
  tmax = std::max(tx1, tx2);

  const double vz = 1.0 / dir.z();
  const double min_rel_z = center.z() - half_flat_to_flat - pos.z();
  const double max_rel_z = center.z() + half_flat_to_flat - pos.z();
  const double tz1 = min_rel_z * vz;
  const double tz2 = max_rel_z * vz;
  tmin = std::max(tmin, std::min(std::min(tz1, tz2), tmax));
  tmax = std::min(tmax, std::max(std::max(tz1, tz2), tmin));

  const double x45 = pos.x() * M_SQRT1_2 + pos.z() * M_SQRT1_2;
  const double xc45 = center.x() * M_SQRT1_2 + center.z() * M_SQRT1_2;
  const double vx45 = 1.0 / (dir.x() * M_SQRT1_2 + dir.z() * M_SQRT1_2);
  const double min_rel_x45 = xc45 - half_flat_to_flat - x45;
  const double max_rel_x45 = xc45 + half_flat_to_flat - x45;
  const double tx45_1 = min_rel_x45 * vx45;
  const double tx45_2 = max_rel_x45 * vx45;
  tmin = std::max(tmin, std::min(std::min(tx45_1, tx45_2), tmax));
  tmax = std::min(tmax, std::max(std::max(tx45_1, tx45_2), tmin));

  const double z45 = pos.z() * M_SQRT1_2 - pos.x() * M_SQRT1_2;
  const double zc45 = center.z() * M_SQRT1_2 - center.x() * M_SQRT1_2;
  const double vz45 = 1.0 / (dir.z() * M_SQRT1_2 - dir.x() * M_SQRT1_2);
  const double min_rel_z45 = zc45 - half_flat_to_flat - z45;
  const double max_rel_z45 = zc45 + half_flat_to_flat - z45;
  const double tz45_1 = min_rel_z45 * vz45;
  const double tz45_2 = max_rel_z45 * vz45;
  tmin = std::max(tmin, std::min(std::min(tz45_1, tz45_2), tmax));
  tmax = std::min(tmax, std::max(std::max(tz45_1, tz45_2), tmin));

  const double vy = 1.0 / dir.y();
  const double min_rel_y = center.y() - half_height - pos.y();
  const double max_rel_y = center.y() + half_height - pos.y();
  const double ty1 = min_rel_y * vy;
  const double ty2 = max_rel_y * vy;
  tmin = std::max(tmin, std::min(std::min(ty1, ty2), tmax));
  tmax = std::min(tmax, std::max(std::max(ty1, ty2), tmin));

  return tmax > std::max(tmin, 0.0);
}

inline bool oct_box_has_future_intersection(
  const Eigen::Vector3d& center, double flat_to_flat, double height,
  const Eigen::Vector3d& pos, const Eigen::Vector3d& dir)
{
  double tmin;
  double tmax;
  return oct_box_has_future_intersection(tmin,tmax,center,flat_to_flat,height,pos,dir);
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

inline void rotation_z_to_xyz_Rzy(Eigen::Matrix3d& m,
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

inline void rotation_y_to_xyz_Ryx(Eigen::Matrix3d& m,
  const double x, const double y, const double z)
{
  double st = std::sqrt(z*z+x*x);
  if(st == 0.0) {
    if(z>=0) {
      m.setIdentity();
    } else {
      m <<  1,  0,  0,
            0, -1,  0,
            0,  0, -1;
    }
  } else {
    double sp = x/st;
    double cp = z/st;
    m <<  cp, x, y*sp,
           0, y,  -st,
         -sp, z, y*cp;
  }
}

inline void rotation_x_to_xyz_Rxz(Eigen::Matrix3d& m,
  const double x, const double y, const double z)
{
  double st = std::sqrt(y*y+z*z);
  if(st == 0.0) {
    if(z>=0) {
      m.setIdentity();
    } else {
      m << -1,  0,  0,
            0, -1,  0,
            0,  0,  1;
    }
  } else {
    double sp = z/st;
    double cp = y/st;
    m <<  x,  -st,   0,
          y, x*cp, -sp,
          z, x*sp,  cp;
  }
}

inline void rotation_z_to_vec_Rzy(Eigen::Matrix3d& m, const Eigen::Vector3d& v)
{
  return rotation_z_to_xyz_Rzy(m, v.x(), v.y(), v.z());
}

inline void rotation_y_to_vec_Ryx(Eigen::Matrix3d& m, const Eigen::Vector3d& v)
{
  return rotation_y_to_xyz_Ryx(m, v.x(), v.y(), v.z());
}

inline void rotation_x_to_vec_Rxz(Eigen::Matrix3d& m, const Eigen::Vector3d& v)
{
  return rotation_x_to_xyz_Rxz(m, v.x(), v.y(), v.z());
}

inline void rotation_z_to_xyz_Rzyz(Eigen::Matrix3d& m,
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
    double st_inv = 1.0/st;
    double sp = y * st_inv;
    double cp = x * st_inv;
    double sp2 = sp*sp;
    double cp2 = cp*cp;
    double zmo_cpsp = (z-1)*cp*sp;
    m << z*cp2+sp2,  zmo_cpsp, x,
          zmo_cpsp, z*sp2+cp2, y,
                -x,        -y, z;
  }
}

inline void rotation_y_to_xyz_Ryxy(Eigen::Matrix3d& m,
  const double x, const double y, const double z)
{
  double st = std::sqrt(z*z+x*x);
  if(st == 0.0) {
    if(y>=0) {
      m.setIdentity();
    } else {
      m <<  1,  0,  0,
            0, -1,  0,
            0,  0, -1;
    }
  } else {
    double st_inv = 1.0/st;
    double sp = x * st_inv;
    double cp = z * st_inv;
    double sp2 = sp*sp;
    double cp2 = cp*cp;
    double ymo_cpsp = (y-1)*cp*sp;

    m << y*sp2+cp2, x,  ymo_cpsp,
                -x, y,        -z,
          ymo_cpsp, z, y*cp2+sp2;
  }
}

inline void rotation_x_to_xyz_Rxzx(Eigen::Matrix3d& m,
  const double x, const double y, const double z)
{
  double st = sqrt(y*y+z*z);
  if(st == 0.0) {
    if(y>=0) {
      m.setIdentity();
    } else {
      m << -1,  0,  0,
            0, -1,  0,
            0,  0,  1;
    }
  } else {
    double st_inv = 1.0/st;
    double sp = z * st_inv;
    double cp = y * st_inv;
    double sp2 = sp*sp;
    double cp2 = cp*cp;
    double xmo_cpsp = (x-1)*cp*sp;

    m << x,        -y,        -z,
         y, x*cp2+sp2,  xmo_cpsp,
         z,  xmo_cpsp, x*sp2+cp2;
  }
}

inline void rotation_z_to_vec_Rzyz(Eigen::Matrix3d& m, const Eigen::Vector3d& v)
{
  return rotation_z_to_xyz_Rzyz(m, v.x(), v.y(), v.z());
}

inline void rotation_y_to_vec_Ryxy(Eigen::Matrix3d& m, const Eigen::Vector3d& v)
{
  return rotation_y_to_xyz_Ryxy(m, v.x(), v.y(), v.z());
}

inline void rotation_x_to_vec_Rxzx(Eigen::Matrix3d& m, const Eigen::Vector3d& v)
{
  return rotation_x_to_xyz_Rxzx(m, v.x(), v.y(), v.z());
}

inline void rotate_in_place_2d(double& x, double& y,
  const double& cos_theta, const double& sin_theta)
{
  double xnew = cos_theta*x - sin_theta*y;
  y = cos_theta*y + sin_theta*x;
  x = xnew;
}

inline void derotate_in_place_2d(double& x, double& y,
  const double& cos_theta, const double& sin_theta)
{
  double xnew = cos_theta*x + sin_theta*y;
  y = cos_theta*y - sin_theta*x;
  x = xnew;
}

inline void rotate_in_place_Rz(Eigen::Vector3d& v,
  const double& cos_theta, const double& sin_theta)
{
  rotate_in_place_2d(v.x(), v.y(), cos_theta, sin_theta);
}

inline void derotate_in_place_Rz(Eigen::Vector3d& v,
  const double& cos_theta, const double& sin_theta)
{
  derotate_in_place_2d(v.x(), v.y(), cos_theta, sin_theta);
}

inline void rotate_in_place_Ry(Eigen::Vector3d& v,
  const double& cos_theta, const double& sin_theta)
{
  rotate_in_place_2d(v.z(), v.x(), cos_theta, sin_theta);
}

inline void derotate_in_place_Ry(Eigen::Vector3d& v,
  const double& cos_theta, const double& sin_theta)
{
  derotate_in_place_2d(v.z(), v.x(), cos_theta, sin_theta);
}

inline void rotate_in_place_Rx(Eigen::Vector3d& v,
  const double& cos_theta, const double& sin_theta)
{
  rotate_in_place_2d(v.y(), v.z(), cos_theta, sin_theta);
}

inline void derotate_in_place_Rx(Eigen::Vector3d& v,
  const double& cos_theta, const double& sin_theta)
{
  derotate_in_place_2d(v.y(), v.z(), cos_theta, sin_theta);
}

inline void rotate_in_place_z_to_u_Rzy(Eigen::Vector3d& v, const Eigen::Vector3d& u)
{
  const double st = sqrt(u.x()*u.x() + u.y()*u.y());
  if(st < std::numeric_limits<double>::epsilon()) {
    return;
  }
  const double st_inv = 1.0/st;
  const double sp = u.y() * st_inv;
  const double cp = u.x() * st_inv;
  rotate_in_place_Ry(v, u.z(), st);
  rotate_in_place_Rz(v, cp, sp);
}

inline void derotate_in_place_z_to_u_Rzy(Eigen::Vector3d& v, const Eigen::Vector3d& u)
{
  const double st = sqrt(u.x()*u.x() + u.y()*u.y());
  if(st < std::numeric_limits<double>::epsilon()) {
    return;
  }
  const double st_inv = 1.0/st;
  const double sp = u.y() * st_inv;
  const double cp = u.x() * st_inv;
  derotate_in_place_Rz(v, cp, sp);
  derotate_in_place_Ry(v, u.z(), st);
}

// Note this function generates scattered vectors with RMS of
// dispersion_per_axis on each of the axes perpendicular to v
inline void
scatter_direction_in_place(Eigen::Vector3d& v, double dispersion_per_axis,
  calin::math::rng::RNG& rng)
{
  Eigen::Vector3d x;
  rng.normal_two_bm(x.x(), x.y());
  x.x() *= dispersion_per_axis;
  x.y() *= dispersion_per_axis;
  x.z() = sqrt(1.0-x.x()*x.x()-x.y()*x.y());

  // Use the simpler rotate fuction (Rzy) as X and Y directions are arbitrary
  rotate_in_place_z_to_u_Rzy(x, v);
  v = x;
}

#ifndef SWIG

// -----------------------------------------------------------------------------
// Skip in SWIG as SWIG output templates map above functions in equivalents of these
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

inline Eigen::Matrix3d rotation_z_to_xyz_Rzy(
  const double x, const double y, const double z)
{
  Eigen::Matrix3d m;
  rotation_z_to_xyz_Rzy(m, x, y, z);
  return m;
}

inline Eigen::Matrix3d rotation_z_to_vec_Rzy(const Eigen::Vector3d v)
{
  Eigen::Matrix3d m;
  rotation_z_to_vec_Rzy(m, v);
  return m;
}

inline Eigen::Matrix3d rotation_y_to_xyz_Ryx(
  const double x, const double y, const double z)
{
  Eigen::Matrix3d m;
  rotation_y_to_xyz_Ryx(m, x, y, z);
  return m;
}

inline Eigen::Matrix3d rotation_y_to_vec_Ryx(const Eigen::Vector3d v)
{
  Eigen::Matrix3d m;
  rotation_y_to_vec_Ryx(m, v);
  return m;
}

inline Eigen::Matrix3d rotation_x_to_xyz_Rxz(
  const double x, const double y, const double z)
{
  Eigen::Matrix3d m;
  rotation_x_to_xyz_Rxz(m, x, y, z);
  return m;
}

inline Eigen::Matrix3d rotation_x_to_vec_Rxz(const Eigen::Vector3d v)
{
  Eigen::Matrix3d m;
  rotation_x_to_vec_Rxz(m, v);
  return m;
}

inline Eigen::Matrix3d rotation_z_to_xyz_Rzyz(
  const double x, const double y, const double z)
{
  Eigen::Matrix3d m;
  rotation_z_to_xyz_Rzyz(m, x, y, z);
  return m;
}

inline Eigen::Matrix3d rotation_z_to_vec_Rzyz(const Eigen::Vector3d v)
{
  Eigen::Matrix3d m;
  rotation_z_to_vec_Rzyz(m, v);
  return m;
}

inline Eigen::Matrix3d rotation_y_to_vec_Ryxy(const Eigen::Vector3d& v)
{
  Eigen::Matrix3d m;
  rotation_y_to_vec_Ryxy(m, v);
  return m;
}

inline Eigen::Matrix3d rotation_x_to_vec_Rxzx(const Eigen::Vector3d& v)
{
  Eigen::Matrix3d m;
  rotation_x_to_vec_Rxzx(m, v);
  return m;
}
#endif

Eigen::Quaterniond euler_to_quaternion(const calin::ix::common_types::EulerAngles3D& euler);
Eigen::Matrix3d euler_to_matrix(const calin::ix::common_types::EulerAngles3D& euler);

void quaternion_to_euler(calin::ix::common_types::EulerAngles3D* euler,
  const Eigen::Quaterniond& q);
void matrix_to_euler(calin::ix::common_types::EulerAngles3D* euler,
  const Eigen::Matrix3d& m);

void scattering_euler(calin::ix::common_types::EulerAngles3D* euler, double dispersion,
  calin::math::rng::RNG& rng, double twist_dispersion = 0);

calin::ix::common_types::EulerAngles3D scattering_euler(double dispersion, calin::math::rng::RNG& rng,
  double twist_dispersion = 0,
  calin::ix::common_types::EulerAngles3D::RotationOrder rotation_order =
    calin::ix::common_types::EulerAngles3D::YXY);

bool euler_is_zero(const calin::ix::common_types::EulerAngles3D& euler);

#ifndef SWIG
inline Eigen::Vector3d norm_and_y_of_polynomial_surface(double& yps, double x, double z,
  const double* p, unsigned np)
{
  double dyps_drho2;
  double rho2 = x*x + z*z;
  calin::math::least_squares::polyval_and_derivative(yps, dyps_drho2, p, np, rho2);
  Eigen::Vector3d norm(2*x*dyps_drho2, -1, 2*z*dyps_drho2);
  norm.normalize();
  return norm;
}

inline Eigen::Vector3d norm_of_polynomial_surface(double x, double z, const double* p, unsigned np)
{
  double yps;
  return norm_and_y_of_polynomial_surface(yps, x, z, p, np);
}
#endif

inline Eigen::Vector3d norm_of_polynomial_surface(double x, double z, const Eigen::VectorXd& p)
{
  return norm_of_polynomial_surface(x, z, p.data(), p.size());
}

inline int find_square_grid_site(double x, double y, double pitch_inv, unsigned nside,
  double xc = 0, double yc = 0, double dead_space_fraction = 0)
{
  const double half_side = 0.5*nside;
  x = (x - xc)*pitch_inv + half_side;
  y = (y - yc)*pitch_inv + half_side;
  const int ux = int(std::floor(x));
  const int uy = int(std::floor(y));
  if(std::min(ux,uy)<0 or std::max(ux,uy)>=int(nside)) {
    return -1;
  }
  if(dead_space_fraction>0) {
    x -= ux;
    y -= uy;
    if(std::min(std::min(x,y), 1.0-std::max(x,y))<dead_space_fraction) {
      return -1;
    }
  }
  return uy*nside+ux;
}

inline bool square_grid_site_center(double& x_out, double& y_out,
  int isite, double pitch, unsigned nside, double xc = 0, double yc = 0)
{
  if((isite<0)or(isite>int(calin::math::special::SQR(nside))))return false;
  const double half_side = 0.5*nside - 0.5;
  div_t div_res = std::div(isite, nside);
  const int ux = div_res.rem;
  const int uy = div_res.quot;
  x_out = (double(ux) - half_side)*pitch + xc;
  y_out = (double(uy) - half_side)*pitch + yc;
  return true;
}

} } } // namespace calin::math::geometry
