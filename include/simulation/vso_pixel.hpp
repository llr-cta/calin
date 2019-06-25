/*

   calin/simulation/vso_pixel.hpp -- Stephen Fegan -- 2015-11-10

   Classes for pixels in the camera

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

//-*-mode:c++; mode:font-lock;-*-

/*! \file VSOPixel.hpp
  Pixel class header file

  \author  Stephen Fegan               \n
           UCLA                        \n
	   sfegan@astro.ucla.edu       \n

  \author  Maciej Nicewicz
           UCLA
	   nicewicz@physics.ucla.edu

  \date    11/30/2004
  \version 0.2
  \note
*/

#pragma once

#include <iostream>

#include <Eigen/Dense>
#include <simulation/vs_optics.pb.h>

namespace calin { namespace simulation { namespace vs_optics {

class VSOTelescope;

//! Class for an individual pixel in the camera
class VSOPixel
{
 public:
  // ************************************************************************
  // Constructor and Destructor
  // ************************************************************************
  VSOPixel();     //!<default constructor
  VSOPixel(const VSOTelescope* T, unsigned PID, unsigned PHID, bool REM,
           const Eigen::Vector3d& P);
  virtual ~VSOPixel();

  // ************************************************************************
  // Dump
  // ************************************************************************

#ifndef SWIG
  calin::ix::simulation::vs_optics::VSOPixelData* dump_as_proto(
    calin::ix::simulation::vs_optics::VSOPixelData* d = nullptr) const;
#else
  calin::ix::simulation::vs_optics::VSOPixelData* dump_as_proto() const;
  void dump_as_proto(calin::ix::simulation::vs_optics::VSOPixelData* d) const;
#endif

  static VSOPixel*
  create_from_proto(const ix::simulation::vs_optics::VSOPixelData& d,
                    const VSOTelescope* T);

  // ************************************************************************
  // Accessors
  // ************************************************************************
  const VSOTelescope*    telescope() const { return fTelescope; }
  bool                   removed() const { return fRemoved; }
  unsigned               id() const { return fID; }
  unsigned               hexID() const { return fHexID; }
  const Eigen::Vector3d& pos() const { return fPos; }

  Eigen::Vector3d incomingSkyVectorAtZenith(double plate_scale=1.0) const;

 private:
  const VSOTelescope*  fTelescope;         //!< Telescope
  unsigned             fID;                //!< Sequential ID
  unsigned             fHexID;             //!< Hex ID
  bool                 fRemoved;           //!< Pixel removed from camera
  Eigen::Vector3d      fPos;               //!< Position
};

} } } // namespace calin::simulation::vs_optics
