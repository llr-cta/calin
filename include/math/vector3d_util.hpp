/*

   calin/math/vector3d_util.hpp -- Stephen Fegan -- 2016-10-23

   Utility functions for Eigen::Vector3d

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <cmath>
#include <calin_global_definitions.hpp>
#include <common_types.pb.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <math/rng.hpp>

namespace calin { namespace math { namespace vector3d_util {

inline void reflect_from_surface(Eigen::Vector3d& v, const Eigen::Vector3d& surface_norm)
{
  assert(fabs(surface_norm.norm() - 1.0) < 1e-6);
  const double mag = 2.0*v.dot(surface_norm);
  v -= surface_norm * mag;
}

calin::ix::common_types::Vector3D* dump_as_proto(const Eigen::Vector3d& v,
  calin::ix::common_types::Vector3D* d = nullptr);

void set_from_proto(Eigen::Vector3d& v, const calin::ix::common_types::Vector3D& d);

inline Eigen::Vector3d from_proto(const calin::ix::common_types::Vector3D& d) {
  Eigen::Vector3d v;
  set_from_proto(v, d);
  return v;
}

calin::ix::common_types::Vector3D* dump_as_scaled_proto(const Eigen::Vector3d& v,
  double scale, calin::ix::common_types::Vector3D* d = nullptr);

void set_from_scaled_proto(Eigen::Vector3d& v,
  const calin::ix::common_types::Vector3D& d, double scale);

inline Eigen::Vector3d from_scaled_proto(const calin::ix::common_types::Vector3D& d, double scale) {
  Eigen::Vector3d v;
  set_from_scaled_proto(v, d, scale);
  return v;
}

void scatter_direction(Eigen::Vector3d& v, double dispersion, math::rng::RNG& rng);

} } } // namespace calin::math::vector3d_util
