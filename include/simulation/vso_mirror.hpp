/* 

   calin/simulation/vso_mirror.hpp -- Stephen Fegan -- 2015-11-05

   Class for mirror data and functionality. This code is derived from
   raytracing code largely implemented by the author at UCLA in
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

/* \file VSOMirror.hpp
  Mirror class header file

  \author   Stephen Fegan             \n
            UCLA                      \n
            sfegan@astro.ucla.edu     \n

  \author   Maciej Nicewicz           \n
            UCLA                      \n
            nicewicz@physics.ucla.edu \n

  \date    11/30/2004
*/

#pragma once

#include <iostream>

#include <math/vs_vec3d.hpp>
#include <math/vs_particle.hpp>
#include <proto/simulation_vs_optics.pb.h>

namespace calin { namespace simulation { namespace vs_optics {

class VSOTelescope;
  
// VSOMirror class for telescope, stores information about individual mirror
class VSOMirror 
{
 public:
  // ************************************************************************
  // Constructor and Destructor
  // ************************************************************************
  VSOMirror();
  VSOMirror(const VSOTelescope* T, unsigned ID, unsigned MHID, bool REM,
            const math::vs_physics::Vec3D& P, const math::vs_physics::Vec3D& A, 
            double FL, double SS, double DF);
  virtual ~VSOMirror();
    
  // ************************************************************************
  // Coordinate transformations
  // ************************************************************************
  // void reflectorToMirror(math::vs_physics::Vec3D& v) const;    //!< Transform vector from reflector to mirror
  // void mirrorToReflector(math::vs_physics::Vec3D& v) const ;   //!< Transform vector from mirror to reflector
  void reflectorToMirror(math::vs_physics::Particle& p) const; //!< Transform particle reflector to mirror
  void mirrorToReflector(math::vs_physics::Particle& p) const; //!< Transform particle mirror to reflector

  math::vs_physics::Vec3D 
  cornerInMirrorCoords(unsigned icorner, double aperture) const;
  math::vs_physics::Vec3D 
  cornerInReflectorCoords(unsigned icorner, double aperture) const;
  
  // ************************************************************************
  // Dump
  // ************************************************************************

  void dump_to_proto(ix::simulation::vs_optics::VSOMirrorData* d) const;
  static VSOMirror*
  create_from_proto(const ix::simulation::vs_optics::VSOMirrorData& d,
                    const VSOTelescope* T);
  
#if 0
  void dumpShort(std::ostream& stream) const;
  void dump(std::ostream& stream, unsigned l=0) const;
  static VSOMirror* createFromShortDump(std::istream& stream,
                                        const VSOTelescope* T);
#endif
  
  // ************************************************************************
  // Accessors
  // ************************************************************************
  const VSOTelescope* telescope() const { return fTelescope; }
  unsigned           id() const { return fID; }
  unsigned           hexID() const { return fHexID; }
  bool               removed() const { return fRemoved; }
  const math::vs_physics::Vec3D& pos() const { return fPos; }
  const math::vs_physics::Vec3D& align() const { return fAlign; }
  double             focalLength() const { return fFocalLength; }
  double             spotSize() const { return fSpotSize; }
  double             degradingFactor() const { return fDegradingFactor; }
  
  // ************************************************************************
  // Setters
  // ************************************************************************

  void realign(const math::vs_physics::Vec3D& a)
  { fAlign = a/a.Norm(); calculateRotationVector(); }
  
 private:
  const VSOTelescope* fTelescope;         //!< Telescope
  unsigned            fID;                //!< Sequential ID (starting at 0)
  unsigned            fHexID;             //!< Hex index on reflector
  bool                fRemoved;           //!< Mirror is removed from scope
  math::vs_physics::Vec3D      fPos;               //!< Position
  math::vs_physics::Vec3D      fAlign;             //!< Alignment angles
  double              fFocalLength;       //!< Focal length
  double              fSpotSize;          //!< Spot size
  double              fDegradingFactor;   //!< Degrading factor of mirror

  math::vs_physics::Vec3D      fRotationVector;
  void calculateRotationVector();
};

} } } // namespace calin::simulation::vso_optics

