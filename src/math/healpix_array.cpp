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

#include <cmath>
#include <chealpix.h>
#include <math/special.hpp>
#include <math/healpix_array.hpp>

using calin::math::special::SQR;

unsigned calin::math::healpix_array::npixel(unsigned nside)
{
  return nside2npix(nside); // =12*SQR(nside);
}

double calin::math::healpix_array::cell_dimension(unsigned nside)
{
  return std::sqrt(4*M_PI/double(npixel(nside)));
}

unsigned calin::math::healpix_array::nside_for_cell_dimension(double dimension)
{
  return std::ceil(std::sqrt(4*M_PI/SQR(dimension)/12.0));
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

unsigned calin::math::healpix_array::vec_to_pixid(unsigned nside,
  const Eigen::Vector3d& vec)
{
  long pixid;
  vec2pix_ring(nside, vec.data(), &pixid);
  return pixid;
}

unsigned calin::math::healpix_array::xyz_to_pixid(unsigned nside,
  double x, double y, double z)
{
  double vec[3];
  vec[0] = x;
  vec[1] = y;
  vec[2] = z;
  long pixid;
  vec2pix_ring(nside, vec, &pixid);
  return pixid;
}

unsigned calin::math::healpix_array::ang_to_pixid(unsigned nside,
  double theta, double phi)
{
  long pixid;
  ang2pix_ring(nside, theta, phi, &pixid);
  return pixid;
}
