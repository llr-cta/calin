/*

   calin/math/healpix_array.cpp -- Stephen Fegan -- 2017-01-10

   Collection of functions which translate between HEALPix and Cartesian
   geometries, and provide other useful calculations for HEALPix grids.

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <cmath>
#include <chealpix/chealpix.h>
#include <math/special.hpp>
#include <math/healpix_array.hpp>

using calin::math::special::SQR;

uint64_t calin::math::healpix_array::npixel(unsigned nside)
{
  return nside2npix64(nside); // =12*SQR(nside);
}

uint64_t calin::math::healpix_array::
npixel_in_ring(unsigned nside, unsigned ringid)
{
  unsigned nring = 4*nside-1;
  assert(ringid < nring);
  return 4*std::min(nside, std::min(ringid+1, nring-ringid));
  // if(ringid < nside)return 4*(ringid+1);
  // else if(ringid < nring-nside)return 4*nside;
  // else return 4*(nring - ringid);
}

uint64_t calin::math::healpix_array::
npixel_contained_by_ring(unsigned nside, unsigned ringid)
{
  unsigned nring = 4*nside-1;
  assert(ringid < nring);
  unsigned npixpolar_2 = nside*(nside+1);
  unsigned ringid_conj = nring-ringid;
  return 2*(std::min((ringid+1)*(ringid+2), npixpolar_2)
    + 2*(std::min(std::max(nside, ringid+1), nring-nside) - nside)*nside
    + (npixpolar_2-std::min((ringid_conj-1)*ringid_conj, npixpolar_2)));
}

uint64_t calin::math::healpix_array::
pixid_to_ringid(unsigned nside, uint64_t pixid)
{
  assert(pixid < npixel(nside));
  unsigned npixpolar = nside*(nside+1)/2;
  pixid /= 4;
  // In the formula below sqrt(1+8*(pixid+0.5))=sqrt(5+8*pixid)
  if(pixid < npixpolar)
    return unsigned(0.5*(std::sqrt(float(5+8*pixid))-1.0));
  pixid -= npixpolar;
  if(pixid < (2*nside-1)*nside)
    return nside + pixid/nside;
  pixid -= (2*nside-1)*nside;
  pixid = npixpolar-pixid-1;
  return 4*nside - 2 - unsigned(0.5*(std::sqrt(float(5+8*pixid))-1.0));
}

double calin::math::healpix_array::cell_area(unsigned nside)
{
  return 4.0*M_PI/double(npixel(nside));
}

double calin::math::healpix_array::cell_dimension(unsigned nside)
{
  return std::sqrt(cell_area(nside));
}

unsigned calin::math::healpix_array::nside_for_cell_dimension(double dimension)
{
  return std::ceil(std::sqrt(4*M_PI/SQR(dimension)/12.0));
}

void calin::math::healpix_array::pixid_to_vec(unsigned nside, uint64_t pixid,
  Eigen::Vector3d& vec)
{
  pix2vec_ring64(nside, pixid, vec.data());
}

void calin::math::healpix_array::pixid_to_xyz(unsigned nside, uint64_t pixid,
  double& x, double& y, double& z)
{
  double vec[3];
  pix2vec_ring64(nside, pixid, vec);
  x = vec[0];
  y = vec[1];
  z = vec[2];
}

void calin::math::healpix_array::pixid_to_ang(unsigned nside, uint64_t pixid,
  double& theta, double& phi)
{
  pix2ang_ring64(nside, pixid, &theta, &phi);
}

uint64_t calin::math::healpix_array::vec_to_pixid(unsigned nside,
  const Eigen::Vector3d& vec)
{
  int64_t pixid;
  vec2pix_ring64(nside, vec.data(), &pixid);
  return pixid;
}

uint64_t calin::math::healpix_array::xyz_to_pixid(unsigned nside,
  double x, double y, double z)
{
  double vec[3];
  vec[0] = x;
  vec[1] = y;
  vec[2] = z;
  int64_t pixid;
  vec2pix_ring64(nside, vec, &pixid);
  return pixid;
}

uint64_t calin::math::healpix_array::ang_to_pixid(unsigned nside,
  double theta, double phi)
{
  int64_t pixid;
  ang2pix_ring64(nside, theta, phi, &pixid);
  return pixid;
}
