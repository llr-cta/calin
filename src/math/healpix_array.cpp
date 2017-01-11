/*

   calin/io/hex_array.cpp -- Stephen Fegan -- 2015-10-21

   Collection of functions which translate between hexagonal and Cartesian
   geometries, and provide other useful calculations for hex grids.

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

#include <chealpix.h>

#include <math/healpix_array.hpp>

unsigned calin::math::healpix_array::npixel(unsigned nside)
{
  return nside2npix(nside);
}

void calin::math::healpix_array::pixid_to_vec(unsigned nside, unsigned pixid,
  Eigen::Vector3d& vec)
{
  pix2vec_ring(nside, pixid, vec.data());
}

void calin::math::healpix_array::pixid_to_xyz(unsigned nside, unsigned pixid,
  double& x, double& y, double& z)
{
  double vec[3];
  pix2vec_ring(nside, pixid, vec);
  x = vec[0];
  y = vec[1];
  z = vec[2];
}

void calin::math::healpix_array::pixid_to_ang(unsigned nside, unsigned pixid,
  double& theta, double& phi)
{
  pix2ang_ring(nside, pixid, &theta, &phi);
}
