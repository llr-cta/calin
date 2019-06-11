/*

   calin/math/ray.cpp -- Stephen Fegan -- 2016-10-24

   Class describing light ray : with position, time, direction and energy.
   A simplification of the old UCLA "Particle" class which hopefully will
   be faster.

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <math/ray.hpp>
#include <math/special.hpp>
#include <util/log.hpp>

using namespace calin::math::ray;
using namespace calin::util::log;
using calin::math::special::SQR;

//! Propagates free particle to the given plane
bool Ray::propagate_to_plane(const Eigen::Vector3d& normal, double d,
  bool time_reversal_ok, double n)
{
  assert(std::fabs(normal.squaredNorm() - 1.0) < 1e-8);

  const double n_hat_dot_v_hat = dir_.dot(normal);

  if(n_hat_dot_v_hat == 0)     //the particle never hits the plane
    return false;

  // Compute the distance between the plane and one parallel to it
  // which goes through the particle's current position

  const double plane_sep = -d - pos_.dot(normal);

  // Distance the particle must travel to reach the plane i.e.
  // n_hat * vhat = cos( theta ) = plane_dist / propagation_dist
  const double propagation_dist = plane_sep / n_hat_dot_v_hat;

  // note: propagation_dist < 0 means particle has passed the plane already
  if((propagation_dist<0)&&(!time_reversal_ok))return false;

  propagate_dist(propagation_dist, n);
  return true;
}

//! Propagates free particle to the given plane
bool Ray::propagate_to_y_plane(double d, bool time_reversal_ok, double n)
{
  const double n_hat_dot_v_hat = dir_.y();

  if(n_hat_dot_v_hat == 0)     //the particle never hits the plane
    return false;

  // Compute the distance between the plane and one parallel to it
  // which goes through the particle's current position

  const double plane_sep = -d - pos_.y();

  // Distance the particle must travel to reach the plane i.e.
  // n_hat * vhat = cos( theta ) = plane_dist / propagation_dist
  const double propagation_dist = plane_sep / n_hat_dot_v_hat;

  // note: propagation_dist < 0 means particle has passed the plane already
  if((propagation_dist<0)&&(!time_reversal_ok))return false;

  propagate_dist(propagation_dist, n);
  return true;
}

//! Propagates free particle to the closest approach with line
bool Ray::propagate_to_point_closest_approach(const Eigen::Vector3d& r0,
  bool time_reversal_ok, double n)
{
  // Distance the particle must travel to reach the closest approach
  double propagation_dist = (pos_-r0).dot(dir_);

  // note: propagation_dist < 0 means particle has passed the plane already
  if((propagation_dist<0)&&(!time_reversal_ok))return false;

  propagate_dist(propagation_dist, n);
  return true;
}

//! Propagates free particle to the closest approach with line
bool Ray::propagate_to_line_closest_approach(const Eigen::Vector3d& normal,
  const Eigen::Vector3d& r0, bool time_reversal_ok, double n)
{
  assert(std::fabs(normal.squaredNorm() - 1.0) < 1e-8);

  const double n_hat_dot_v_hat = dir_.dot(normal);

  const double denom = 1.0-SQR(n_hat_dot_v_hat);
  if(denom<=0)return true; // Ray parallel to line; always at closest approach

  // Distance the particle must travel to reach the closest approach
  double propagation_dist = (r0-pos_).dot(normal*n_hat_dot_v_hat-dir_)/denom;

  // note: propagation_dist < 0 means particle has passed the plane already
  if((propagation_dist<0)&&(!time_reversal_ok))return false;

  propagate_dist(propagation_dist, n);
  return true;
}

bool Ray::propagate_to_standard_sphere_1st_interaction_fwd_bwd(double radius,
  double n)
{
  Eigen::Vector3d pos_rel(pos_.x(), pos_.y()-radius, pos_.z());
  const double b_2 = pos_rel.dot(dir_);
  const double c = pos_rel.squaredNorm() - SQR(radius);

  const double disc_4 = SQR(b_2)-c;
  if(disc_4 < 0)return false;

  const double time = -std::sqrt(disc_4) - b_2;
  propagate_dist(time, n);
  return true;
}


bool Ray::propagate_to_standard_sphere_1st_interaction_fwd_only(double radius,
  double n)
{
  Eigen::Vector3d pos_rel(pos_.x(), pos_.y()-radius, pos_.z());
  const double b_2 = pos_rel.dot(dir_);
  const double c = pos_rel.squaredNorm() - SQR(radius);

  const double disc_4 = SQR(b_2)-c;
  if(disc_4 < 0)return false;

  const double time = -std::sqrt(disc_4) - b_2;

  if(time>0) {
    propagate_dist(time, n);
    return true;
  }

  return false;
}

bool Ray::propagate_to_y_sphere_1st_interaction_fwd_bwd(double radius,
  double surface_y_min, double n)
{
  Eigen::Vector3d pos_rel(pos_.x(), pos_.y()-(radius+surface_y_min), pos_.z());
  const double b_2 = pos_rel.dot(dir_);
  const double c = pos_rel.squaredNorm() - SQR(radius);

  const double disc_4 = SQR(b_2)-c;
  if(disc_4 < 0)return false;

  const double time = -std::sqrt(disc_4) - b_2;
  propagate_dist(time, n);
  return true;
}

bool Ray::propagate_to_y_sphere_1st_interaction_fwd_only(double radius,
  double surface_y_min, double n)
{
  Eigen::Vector3d pos_rel(pos_.x(), pos_.y()-(radius+surface_y_min), pos_.z());
  const double b_2 = pos_rel.dot(dir_);
  const double c = pos_rel.squaredNorm() - SQR(radius);

  const double disc_4 = SQR(b_2)-c;
  if(disc_4 < 0)return false;

  const double time = -std::sqrt(disc_4) - b_2;

  if(time>0) {
    propagate_dist(time, n);
    return true;
  }

  return false;
}

bool Ray::propagate_to_standard_sphere_2nd_interaction_fwd_only(double radius,
  double n)
{
  Eigen::Vector3d pos_rel(pos_.x(), pos_.y()-radius, pos_.z());
  // double a = 1.0
  const double b_2 = pos_rel.dot(dir_);
  const double c = pos_rel.squaredNorm() - SQR(radius);

  const double disc_4 = SQR(b_2)-c;
  if(disc_4 < 0)return false;

#if 0
  // This is the "correct" way to solve the roots, as described by NR etc.
  // But the branch is slow, so it's faster to simply use the conventional
  // formula and accept that there may be some inaccuacy in the solution
  double time;
  if(b_2<0.0) {
    const double q = -(b_2-std::sqrt(disc_4));
    time = q;
  } else {
    const double q = -(b_2+std::sqrt(disc_4));
    time = c/q;
  }
#else
  const double time = std::sqrt(disc_4) - b_2;
#endif

  if(time>0) {
    propagate_dist(time, n);
    return true;
  }

  return false;
}

bool Ray::propagate_to_standard_sphere_2nd_interaction_fwd_bwd(double radius,
  double n)
{
  Eigen::Vector3d pos_rel(pos_.x(), pos_.y()-radius, pos_.z());
  const double b_2 = pos_rel.dot(dir_);
  const double c = pos_rel.squaredNorm() - SQR(radius);

  const double disc_4 = SQR(b_2)-c;
  if(disc_4 < 0)return false;

  const double time = std::sqrt(disc_4) - b_2;
  propagate_dist(time, n);
  return true;
}

bool Ray::propagate_to_y_sphere_2nd_interaction_fwd_bwd(double radius,
  double surface_y_min, double n)
{
  Eigen::Vector3d pos_rel(pos_.x(), pos_.y()-(radius+surface_y_min), pos_.z());
  const double b_2 = pos_rel.dot(dir_);
  const double c = pos_rel.squaredNorm() - SQR(radius);

  const double disc_4 = SQR(b_2)-c;
  if(disc_4 < 0)return false;

  const double time = std::sqrt(disc_4) - b_2;
  propagate_dist(time, n);
  return true;
}

bool Ray::propagate_to_y_sphere_2nd_interaction_fwd_only(double radius,
  double surface_y_min, double n)
{
  Eigen::Vector3d pos_rel(pos_.x(), pos_.y()-(radius+surface_y_min), pos_.z());
  const double b_2 = pos_rel.dot(dir_);
  const double c = pos_rel.squaredNorm() - SQR(radius);

  const double disc_4 = SQR(b_2)-c;
  if(disc_4 < 0)return false;

  const double time = std::sqrt(disc_4) - b_2;

  if(time>0) {
    propagate_dist(time, n);
    return true;
  }

  return false;
}

Ray::IPOut Ray::propagate_to_sphere(const Eigen::Vector3d& center, double radius,
  IntersectionPoint ip, bool time_reversal_ok, double n)
{
  // Vector from here to center of circle
  Eigen::Vector3d to_center = center - pos_;

  // Distance^2 along particle trajectory to closest approach with center
  double prop_dist = dir_.dot(to_center);
  double prop_dist_2 = SQR(prop_dist);

  // Distance^2 to closest approach with center along tangent (negative inside)
  double tangent_dist_2 = to_center.squaredNorm() - SQR(radius);

  if(prop_dist_2 < tangent_dist_2)   //no intersection
    return IPO_NONE;

  //time of first intersection
  double t1 = prop_dist - sqrt(prop_dist_2 - tangent_dist_2);
  //time of second intersection (t2>t1)
  double t2 = prop_dist + sqrt(prop_dist_2 - tangent_dist_2);

  double time = 0;

  // Select which intersection point to go to depending on what we were
  // asked to do and whether we can move back in time... otherwise do not
  // propagate

  IPOut ipo = IPO_NONE;
  switch(ip)
  {
  case IP_CLOSEST:
    if(time_reversal_ok) {
      if(fabs(t1)<fabs(t2))ipo=IPO_FIRST;
      else ipo=IPO_SECOND;
    } else {
      if(t1>=0)ipo=IPO_FIRST;
      else if(t2>=0)ipo=IPO_SECOND;
      else return IPO_NONE;
    }
    break;

  case IP_FARTHEST:
    if(time_reversal_ok) {
      if(fabs(t1)<fabs(t2))ipo=IPO_SECOND;
      else ipo=IPO_FIRST;
    } else {
      if(t2>=0)ipo=IPO_SECOND;
      else return IPO_NONE;
    }
    break;

  case IP_NEXT:
    if(t1>0)ipo=IPO_FIRST;
    else if(t2>0)ipo=IPO_SECOND;
    else return IPO_NONE;
    break;

  case IP_PREVIOUS:
    if(!time_reversal_ok)return IPO_NONE;
    else if(t2<0)ipo=IPO_SECOND;
    else if(t1<0)ipo=IPO_FIRST;
    else return IPO_NONE;
    break;

  case IP_EARLIEST:
    if((time_reversal_ok)||(t1>=0))ipo=IPO_FIRST;
    else return IPO_NONE;

  case IP_LATEST:
    if((time_reversal_ok)||(t2>=0))ipo=IPO_SECOND;
    else return IPO_NONE;
  };

  switch(ipo)
  {
  case IPO_FIRST:
    time = t1;
    break;
  case IPO_SECOND:
    time = t2;
    break;
  case IPO_NONE:
    assert(0);
  }

  propagate_dist(time, n);
  return ipo;
}

Ray::IPOut Ray::propagate_to_cylinder(const Eigen::Vector3d& center,
  const Eigen::Vector3d& normal, double radius,
  IntersectionPoint ip, bool time_reversal_ok, double n)
{
  Eigen::Vector3d dr = pos_ - center;  // Current (initial) position of particle

  double normal_dot_dir = normal.dot(dir_);
  double dr_dot_normal = dr.dot(normal);

  double a = 1.0-SQR(normal_dot_dir);
  double b = 2.0*(dr.dot(dir_)-dr_dot_normal*normal_dot_dir);
  double c = dr.squaredNorm()-SQR(dr_dot_normal)-SQR(radius);

  double disc = SQR(b)-4.0*a*c;
  if(disc < 0)   //no intersection
    return IPO_NONE;

  double q;
  if(b<0.0) q = -0.5*(b-std::sqrt(disc));
  else q = -0.5*(b+std::sqrt(disc));

  //time of first intersection
  double t1 = c/q;
  //time of second intersection (t2>t1)
  double t2 = q/a;
  if(t2<t1)std::swap(t1,t2);

  double time = 0;

  // Select which intersection point to go to depending on what we were
  // asked to do and whether we can move back in time... otherwise do not
  // propagate

  IPOut ipo = IPO_NONE;
  switch(ip)
  {
  case IP_CLOSEST:
    if(time_reversal_ok) {
      if(fabs(t1)<fabs(t2))ipo=IPO_FIRST;
    	else ipo=IPO_SECOND;
    } else {
      if(t1>=0)ipo=IPO_FIRST;
	    else if(t2>=0)ipo=IPO_SECOND;
	    else return IPO_NONE;
    }
    break;

  case IP_FARTHEST:
    if(time_reversal_ok) {
      if(fabs(t1)<fabs(t2))ipo=IPO_SECOND;
      else ipo=IPO_FIRST;
    } else {
      if(t2>=0)ipo=IPO_SECOND;
      else return IPO_NONE;
    }
    break;

  case IP_NEXT:
    if(t1>0)ipo=IPO_FIRST;
    else if(t2>0)ipo=IPO_SECOND;
    else return IPO_NONE;
    break;

  case IP_PREVIOUS:
    if(!time_reversal_ok)return IPO_NONE;
    else if(t2<0)ipo=IPO_SECOND;
    else if(t1<0)ipo=IPO_FIRST;
    else return IPO_NONE;
    break;

  case IP_EARLIEST:
    if((time_reversal_ok)||(t1>=0))ipo=IPO_FIRST;
    else return IPO_NONE;

  case IP_LATEST:
    if((time_reversal_ok)||(t2>=0))ipo=IPO_SECOND;
    else return IPO_NONE;
  };

  switch(ipo)
  {
  case IPO_FIRST:
    time = t1;
    break;
  case IPO_SECOND:
    time = t2;
    break;
  case IPO_NONE:
    assert(0);
  }

  propagate_dist(time, n);
  return ipo;
}
