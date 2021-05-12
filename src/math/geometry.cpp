/*

   calin/math/geometry.cpp -- Stephen Fegan -- 2021-04-30

   Functions for geometrical manuipulations

   Copyright 2021, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <numeric>
#include <cmath>

#include <math/geometry.hpp>

using namespace calin::math::geometry;

Eigen::Quaterniond calin::math::geometry::
euler_to_quaternion(const calin::ix::common_types::EulerAngles3D& euler)
{
  using Eigen::Quaterniond;
  using Eigen::AngleAxisd;
  using Eigen::Vector3d;
  Eigen::Quaterniond q;
  switch(euler.rotation_order()) {
  case calin::ix::common_types::EulerAngles3D::ZXZ:
    return AngleAxisd(euler.alpha()*(M_PI/180.0), Vector3d::UnitZ())
          * AngleAxisd(euler.beta()*(M_PI/180.0), Vector3d::UnitX())
          * AngleAxisd(euler.gamma()*(M_PI/180.0), Vector3d::UnitZ());
  case calin::ix::common_types::EulerAngles3D::XYX:
    return AngleAxisd(euler.alpha()*(M_PI/180.0), Vector3d::UnitX())
          * AngleAxisd(euler.beta()*(M_PI/180.0), Vector3d::UnitY())
          * AngleAxisd(euler.gamma()*(M_PI/180.0), Vector3d::UnitX());
  case calin::ix::common_types::EulerAngles3D::YZY:
    return AngleAxisd(euler.alpha()*(M_PI/180.0), Vector3d::UnitY())
          * AngleAxisd(euler.beta()*(M_PI/180.0), Vector3d::UnitZ())
          * AngleAxisd(euler.gamma()*(M_PI/180.0), Vector3d::UnitY());
  case calin::ix::common_types::EulerAngles3D::ZYZ:
    return AngleAxisd(euler.alpha()*(M_PI/180.0), Vector3d::UnitZ())
          * AngleAxisd(euler.beta()*(M_PI/180.0), Vector3d::UnitY())
          * AngleAxisd(euler.gamma()*(M_PI/180.0), Vector3d::UnitZ());
  case calin::ix::common_types::EulerAngles3D::XZX:
    return AngleAxisd(euler.alpha()*(M_PI/180.0), Vector3d::UnitX())
          * AngleAxisd(euler.beta()*(M_PI/180.0), Vector3d::UnitZ())
          * AngleAxisd(euler.gamma()*(M_PI/180.0), Vector3d::UnitX());
  case calin::ix::common_types::EulerAngles3D::YXY:
    return AngleAxisd(euler.alpha()*(M_PI/180.0), Vector3d::UnitY())
          * AngleAxisd(euler.beta()*(M_PI/180.0), Vector3d::UnitX())
          * AngleAxisd(euler.gamma()*(M_PI/180.0), Vector3d::UnitY());
  default:
    throw std::runtime_error("Unsupported rotation order");
  }
}

Eigen::Matrix3d calin::math::geometry::
euler_to_matrix(const calin::ix::common_types::EulerAngles3D& euler)
{
  return euler_to_quaternion(euler).toRotationMatrix();
}

void calin::math::geometry::quaternion_to_euler(
  calin::ix::common_types::EulerAngles3D* euler, const Eigen::Quaterniond& q)
{
  return matrix_to_euler(euler, q.toRotationMatrix());
}

void calin::math::geometry::matrix_to_euler(
  calin::ix::common_types::EulerAngles3D* euler, const Eigen::Matrix3d& m)
{
  Eigen::Vector3d v;
  switch(euler->rotation_order()) {
  case calin::ix::common_types::EulerAngles3D::ZXZ:
    v = m.eulerAngles(2,0,2); break;
  case calin::ix::common_types::EulerAngles3D::XYX:
    v = m.eulerAngles(0,1,0); break;
  case calin::ix::common_types::EulerAngles3D::YZY:
    v = m.eulerAngles(1,2,1); break;
  case calin::ix::common_types::EulerAngles3D::ZYZ:
    v = m.eulerAngles(2,1,2); break;
  case calin::ix::common_types::EulerAngles3D::XZX:
    v = m.eulerAngles(0,2,0); break;
  case calin::ix::common_types::EulerAngles3D::YXY:
    v = m.eulerAngles(1,0,1); break;
  default:
    throw std::runtime_error("Unsupported rotation order");
  }
  if(v(1) < 0) {
    v(0) = std::fmod(v(0) + 2*M_PI, 2*M_PI) - M_PI;
    v(1) = std::fabs(v(1));
    v(2) = std::fmod(v(2) + 2*M_PI, 2*M_PI) - M_PI;
  }
  euler->set_alpha(v(0) * 180.0/M_PI);
  euler->set_beta(v(1) * 180.0/M_PI);
  euler->set_gamma(v(2) * 180.0/M_PI);
}

void calin::math::geometry::scattering_euler(
  calin::ix::common_types::EulerAngles3D* euler, double dispersion, math::rng::RNG& rng,
  double twist_dispersion)
{
  double alpha = 0;
  double beta = 0;
  double gamma = 0;

  if(dispersion>0)
  {
    alpha = rng.uniform() * 360.0 - 180.0;
    beta = dispersion*std::sqrt(-2.0*std::log(rng.uniform()));
    gamma = -alpha;
  }

  if(twist_dispersion > 0) {
    gamma = std::fmod(std::fmod(gamma + rng.uniform()*twist_dispersion + 180, 360) + 360, 360) - 180;
  }

  euler->set_alpha(alpha);
  euler->set_beta(beta);
  euler->set_gamma(gamma);
}

calin::ix::common_types::EulerAngles3D
calin::math::geometry::scattering_euler(double dispersion, calin::math::rng::RNG& rng,
  double twist_dispersion, calin::ix::common_types::EulerAngles3D::RotationOrder rotation_order)
{
  calin::ix::common_types::EulerAngles3D euler;
  euler.set_rotation_order(rotation_order);
  scattering_euler(&euler, dispersion, rng, twist_dispersion);
  return euler;
}

bool calin::math::geometry::euler_is_zero(
  const calin::ix::common_types::EulerAngles3D& euler)
{
  return euler.alpha() == euler.beta() == euler.gamma() == 0;
}
