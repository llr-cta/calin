/*

   calin/simulation/vso_mirror.hpp -- Stephen Fegan -- 2015-11-05

   Class for mirror data and functionality. This code is derived from
   raytracing code largely implemented by the author at UCLA in
   2004. It is not complient with the calin coding style.

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

#include <math/ray.hpp>
#include <Eigen/Dense>
#include <simulation/vs_optics.pb.h>

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
            const Eigen::Vector3d& P, const Eigen::Vector3d& A,
            double FL, double SS, double DF);
  virtual ~VSOMirror();

  // ************************************************************************
  // Coordinate transformations
  // ************************************************************************
#if 0
  void reflectorToMirror(Eigen::Vector3d& v) const;    //!< Transform vector from reflector to mirror
  void mirrorToReflector(Eigen::Vector3d& v) const ;   //!< Transform vector from mirror to reflector
#endif
  void reflectorToMirror(math::ray::Ray& r) const; //!< Transform ray reflector to mirror
  void mirrorToReflector(math::ray::Ray& r) const; //!< Transform ray mirror to reflector

  Eigen::Vector3d cornerInMirrorCoords(unsigned icorner, double aperture) const;
  Eigen::Vector3d cornerInReflectorCoords(unsigned icorner, double aperture) const;

  // ************************************************************************
  // Dump
  // ************************************************************************

#ifndef SWIG
  calin::ix::simulation::vs_optics::VSOMirrorData* dump_as_proto(
    calin::ix::simulation::vs_optics::VSOMirrorData* d = nullptr) const;
#else
  calin::ix::simulation::vs_optics::VSOMirrorData* dump_as_proto() const;
  void dump_as_proto(calin::ix::simulation::vs_optics::VSOMirrorData* d) const;
#endif

  static VSOMirror*
  create_from_proto(const ix::simulation::vs_optics::VSOMirrorData& d,
                    const VSOTelescope* T);

  // ************************************************************************
  // Accessors
  // ************************************************************************
  const VSOTelescope*    telescope() const { return fTelescope; }
  unsigned               id() const { return fID; }
  unsigned               hexID() const { return fHexID; }
  bool                   removed() const { return fRemoved; }
  const Eigen::Vector3d& pos() const { return fPos; }
  const Eigen::Vector3d& align() const { return fAlign; }
  double                 focalLength() const { return fFocalLength; }
  double                 spotSize() const { return fSpotSize; }
  double                 degradingFactor() const { return fDegradingFactor; }

  // ************************************************************************
  // Setters
  // ************************************************************************

  void realign(const Eigen::Vector3d& a)
  { fAlign = a.normalized(); calculateRotationVector(); }

 private:
  const VSOTelescope*  fTelescope;         //!< Telescope
  unsigned             fID;                //!< Sequential ID (starting at 0)
  unsigned             fHexID;             //!< Hex index on reflector
  bool                 fRemoved;           //!< Mirror is removed from scope
  Eigen::Vector3d      fPos;               //!< Position
  Eigen::Vector3d      fAlign;             //!< Alignment angles
  double               fFocalLength;       //!< Focal length
  double               fSpotSize;          //!< Spot size
  double               fDegradingFactor;   //!< Degrading factor of mirror

  Eigen::Matrix3d      fRotationVector;
  void calculateRotationVector();
};

} } } // namespace calin::simulation::vso_optics
