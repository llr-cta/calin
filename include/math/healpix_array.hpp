/*

   calin/io/healpix_array.hpp -- Stephen Fegan -- 2017-01-10

   Collection of functions which translate between HEALPix and Cartesian
   geometries, and provide other useful calculations for HEALPix grids.

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <Eigen/Core>

namespace calin { namespace math { namespace healpix_array {

// *****************************************************************************
//
// General functions for the HEALPix array
//
// *****************************************************************************

unsigned npixel(unsigned nside);
void pixid_to_vec(unsigned nside, unsigned pixid, Eigen::Vector3d& vec);
void pixid_to_xyz(unsigned nside, unsigned pixid, double& x, double& y, double& z);
void pixid_to_ang(unsigned nside, unsigned pixid, double& theta, double& phi);

} } } // namespace calin::math::healpix_array
