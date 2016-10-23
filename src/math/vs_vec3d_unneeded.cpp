/*

   calin/math/vs_vec3d.cpp -- Stephen Fegan -- 2016-10-22

   Class for 3-vector operations. This code is a shallow layer around the
   Eigen::Vector3d class proving some level of compatibility with the Vec3D
   class from simulation code largely implemented by the author at UCLA in
   2004. It is not complient with the calin coding style.

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

#include <math/vs_vec3d.hpp>

using namespace calin::math::vs_physics;

calin::ix::common_types::Vector3D* Vec3D::
dump_as_proto(ix::common_types::Vector3D* d) const
{
  if(d == nullptr)d = new ix::common_types::Vector3D;
  d->set_x(x());
  d->set_y(y());
  d->set_z(z());
  return d;
}

void Vec3D::set_from_proto(const ix::common_types::Vector3D& d)
{
  x() = d.x();
  y() = d.y();
  z() = d.z();
}

calin::ix::common_types::Vector3D* Vec3D::
dump_scaled_as_proto(double scale, ix::common_types::Vector3D* d) const
{
  if(d == nullptr)d = new ix::common_types::Vector3D;
  d->set_x(x()*scale);
  d->set_y(y()*scale);
  d->set_z(z()*scale);
  return d;
}

void Vec3D::set_from_scaled_proto(const ix::common_types::Vector3D& d, double scale)
{
  x() = d.x()*scale;
  y() = d.y()*scale;
  z() = d.z()*scale;
}
