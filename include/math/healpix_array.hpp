/*

   calin/math/healpix_array.hpp -- Stephen Fegan -- 2017-01-10

   Collection of functions which translate between HEALPix and Cartesian
   geometries, and provide other useful calculations for HEALPix grids.

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <Eigen/Core>

namespace calin { namespace math { namespace healpix_array {

// *****************************************************************************
//
// General functions for the HEALPix array
//
// *****************************************************************************

uint64_t npixel(unsigned nside);
inline unsigned nring(unsigned nside) { return 4*nside-1; }
uint64_t npixel_in_ring(unsigned nside, unsigned ringid);
uint64_t npixel_contained_by_ring(unsigned nside, unsigned ringid);
uint64_t pixid_to_ringid(unsigned nside, uint64_t pixid);

double cell_area(unsigned nside);
double cell_dimension(unsigned nside);
unsigned nside_for_cell_dimension(double dimension);

void pixid_to_vec(unsigned nside, uint64_t pixid, Eigen::Vector3d& vec);
void pixid_to_xyz(unsigned nside, uint64_t pixid, double& x, double& y, double& z);
void pixid_to_ang(unsigned nside, uint64_t pixid, double& theta, double& phi);

uint64_t vec_to_pixid(unsigned nside, const Eigen::Vector3d& vec);
uint64_t xyz_to_pixid(unsigned nside, double x, double y, double z);
uint64_t ang_to_pixid(unsigned nside, double theta, double phi);

} } } // namespace calin::math::healpix_array
