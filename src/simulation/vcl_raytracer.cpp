/*

   calin/simulation/vcl_raytracer.cpp -- Stephen Fegan -- 2020-11-30

   Some test functions for VCL raytracer components

   Copyright 2020, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <simulation/vcl_raytracer.hpp>

bool calin::simulation::vcl_raytracer::
test_tube_obscuration(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2, double radius,
  const Eigen::Vector3d& r, const Eigen::Vector3d& u)
{
  using Real = calin::util::vcl::VCL256DoubleReal;
  using vec3_vt = typename Real::vec3_vt;
  using real_vt = typename Real::real_vt;
  VCLTubeObscuration<Real> obs(x1, x2, radius);
  vec3_vt rv = r.template cast<real_vt>();
  vec3_vt uv = u.template cast<real_vt>();
  calin::math::ray::VCLRay<Real> ray(rv, uv);
  calin::math::ray::VCLRay<Real> ray_out;
  auto obscure = obs.doesObscure(ray, ray_out, 1.0);
  return vcl::horizontal_and(obscure);
}
