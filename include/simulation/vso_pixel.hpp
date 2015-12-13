/* 

   calin/simulation/vso_pixel.hpp -- Stephen Fegan -- 2015-11-10

   Classes for pixels in the camera

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

#include <math/vs_vec3d.hpp>
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
           const math::vs_physics::Vec3D& P);
  virtual ~VSOPixel();
        
  // ************************************************************************
  // Dump
  // ************************************************************************

  void dump_to_proto(ix::simulation::vs_optics::VSOPixelData* d) const;
  static VSOPixel*
  create_from_proto(const ix::simulation::vs_optics::VSOPixelData& d,
                    const VSOTelescope* T);

#if 0
  void dumpShort(std::ostream& stream) const;
  void dump(std::ostream& stream, unsigned l=0) const;
  static VSOPixel* createFromShortDump(std::istream& stream,
                                       const VSOTelescope* T);
#endif
  
  // ************************************************************************
  // Accessors
  // ************************************************************************
  const VSOTelescope* telescope() const { return fTelescope; }
  bool                removed() const { return fRemoved; }
  unsigned            id() const { return fID; }
  unsigned            hexID() const { return fHexID; }
  const math::vs_physics::Vec3D& pos() const { return fPos; }

  math::vs_physics::Vec3D      incomingSkyVectorAtZenith(double plate_scale=1.0) 
      const;
    
 private:
  const VSOTelescope* fTelescope;         //!< Telescope 
  unsigned            fID;                //!< Hex ID
  unsigned            fHexID;             //!< Hex ID
  bool                fRemoved;           //!< Pixel removed from camera
  math::vs_physics::Vec3D      fPos;               //!< Position
};

} } } // namespace calin::simulation::vs_optics
