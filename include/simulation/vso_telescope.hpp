/*

   calin/simulation/vso_telescope.hpp -- Stephen Fegan -- 2015-11-09

   Class for telescope, a container for obscurations, mirrors and pixels

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

/*! \file VSOTelescope.hpp
  Reflector class header file

  \author  Stephen Fegan               \n
           UCLA                        \n
	   sfegan@astro.ucla.edu       \n

  \author  Maciej Nicewicz             \n
           UCLA                        \n
	   nicewicz@physics.ucla.edu   \n

  \date    11/30/2004
  \version 0.2
  \note
*/

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Dense>
//#include <math/vs_particle.hpp>
#include <math/rng.hpp>

#include <simulation/vs_optics.pb.h>
#include <simulation/vso_mirror.hpp>
#include <simulation/vso_pixel.hpp>
#include <simulation/vso_obscuration.hpp>

namespace calin { namespace simulation { namespace vs_optics {

/*! \class VSOTelescope
  \brief VSOTelescope class;
*/

class VSOTelescope
{
 public:
  VSOTelescope();
  VSOTelescope(unsigned TID, /*unsigned THID,*/ const Eigen::Vector3d&P,
               double DY, double AX, double AY, double EL, double AZ,
               const Eigen::Vector3d& T, double CR, double A,
               double FSP, double FS, double RR,
               unsigned HRN, double RIP, bool MP,
               const Eigen::Vector3d& FPT, double CD, double FOV, double D,
               double PS, double PR, double CSP, const Eigen::Vector3d& FPR,
               double CIP, bool PP,
               const std::vector<VSOObscuration*>& OBSVEC = {}
               );
  VSOTelescope(const VSOTelescope& o);
  virtual ~VSOTelescope();

#ifndef SWIG
  const VSOTelescope& operator =(const VSOTelescope& o);
#endif

  // ************************************************************************
  // Create a new telescope randomly using parameters
  // ************************************************************************

  void populateMirrorsAndPixelsRandom(
      const ix::simulation::vs_optics::IsotropicDCArrayParameters& param,
      math::rng::RNG& rng);

  void add_obscuration(VSOObscuration* obs) {
    fObscurations.push_back(obs); }
  void add_mirror(VSOMirror* mirror) {
    if(fMirrors.size()<=mirror->id())fMirrors.resize(mirror->id()+1);
    fMirrors[mirror->id()] = mirror;
    if(fMirrorsByHexID.size()<=mirror->hexID())
      fMirrorsByHexID.resize(mirror->hexID()+1);
    fMirrorsByHexID[mirror->hexID()] = mirror;
  }
  void add_pixel(VSOPixel* pix) {
    if(fPixels.size()<=pix->id())fPixels.resize(pix->id()+1);
    fPixels[pix->id()] = pix;
    if(fPixelsByHexID.size()<=pix->hexID())
      fPixelsByHexID.resize(pix->hexID()+1);
    fPixelsByHexID[pix->hexID()] = pix;
  }

  // ************************************************************************
  // Telescope repointing
  // ************************************************************************

  //! Points telescope along axis
  bool pointTelescope(const Eigen::Vector3d& v);
  bool pointTelescopeAzEl(const double az_rad, const double el_rad);

  //! Transform position-like vector from global to reflector
  void globalToReflector_pos(Eigen::Vector3d& v) const;
  //! Transform momentum-like vector from global to reflector
  void globalToReflector_mom(Eigen::Vector3d& v) const;
  //! Transform position-like vector from reflector to global
  void reflectorToGlobal_pos(Eigen::Vector3d& v) const;
  //! Transform momentum-like vector from reflector to global
  void reflectorToGlobal_mom(Eigen::Vector3d& v) const;

  //! Transform position-like vector from global to reflector
  void focalPlaneToReflector_pos(Eigen::Vector3d& v) const;
  //! Transform momentum-like vector from global to reflector
  void focalPlaneToReflector_mom(Eigen::Vector3d& v) const;
  //! Transform position-like vector from reflector to global
  void reflectorToFocalPlane_pos(Eigen::Vector3d& v) const;
  //! Transform momentum-like vector from reflector to global
  void reflectorToFocalPlane_mom(Eigen::Vector3d& v) const;

  //! Transform particle global to reflector
  void globalToReflector(calin::math::ray::Ray& r) const;
  //! Transform particle reflector to global
  void reflectorToGlobal(calin::math::ray::Ray& r) const;
  //! Transform from focal plane to reflector
  void focalPlaneToReflector(calin::math::ray::Ray& r) const;
  //! Transform from reflector to focal plane
  void reflectorToFocalPlane(calin::math::ray::Ray& r) const;

  // ************************************************************************
  // Dump and Reload
  // ************************************************************************

#ifndef SWIG
  calin::ix::simulation::vs_optics::VSOTelescopeData* dump_as_proto(
    calin::ix::simulation::vs_optics::VSOTelescopeData* d = nullptr) const;
#else
  calin::ix::simulation::vs_optics::VSOTelescopeData* dump_as_proto() const;
  void dump_as_proto(calin::ix::simulation::vs_optics::VSOTelescopeData* d) const;
#endif
  static VSOTelescope*
  create_from_proto(const ix::simulation::vs_optics::VSOTelescopeData& d);

  // ************************************************************************
  // Accessors
  // ************************************************************************

  unsigned               id() const { return fID; }
  double                 altitude() const { return fPos.z(); }
  const Eigen::Vector3d& pos() const { return fPos; }
  const Eigen::Vector3d& position() const { return fPos; }
  double                 deltaY() const { return fDeltaY; }
  double                 alphaX() const { return fAlphaX; }
  double                 alphaY() const { return fAlphaY; }
  double                 azimuth() const { return fAzimuth; }
  double                 elevation() const { return fElevation; }
  Eigen::Vector3d        opticalAxis() const;

  const Eigen::Vector3d& translation() const { return fTranslation; }
  double                 curvatureRadius() const { return fCurvatureRadius; }
  double                 aperture() const { return fAperture; }
  double                 facetSpacing() const { return  fFacetSpacing; }
  double                 facetSize() const { return  fFacetSize; }
  double                 reflectorRotation() const { return fReflectorRotation; }
  double                 cosReflectorRotation() const { return fCosReflectorRotation; }
  double                 sinReflectorRotation() const { return fSinReflectorRotation; }
  unsigned               mirrorHexRings() const { return fHexagonRingsN; }
  double                 reflectorIP() const { return fReflectorIP; }
  bool                   mirrorParity() const { return fMirrorParity; }

  const Eigen::Vector3d& focalPlanePosition() const { return fFPTranslation; }
  double                 cameraDiameter() const { return fCameraDiameter; }
  double                 fov() const { return fFieldOfView; }
  double                 cathodeDiameter() const { return fCathodeDiameter; }
  double                 pixelSpacing() const { return fPixelSpacing; }
  double                 pixelRotation() const { return fPixelRotation; }
  double                 cosPixelRotation() const { return fCosPixelRotation; }
  double                 sinPixelRotation() const { return fSinPixelRotation; }
  double                 concentratorSurvivalProb() const { return fConcSurvProb; }
  const Eigen::Vector3d& focalPlaneRotion() const { return fFPRotation; }
  double                 cameraIP() const { return fCameraIP; }
  bool                   pixelParity() const { return fPixelParity; }

  unsigned               numObscurations() const { return fObscurations.size(); }

  unsigned               numMirrors() const { return fMirrors.size(); }
  unsigned               numMirrorHexSites() const { return fMirrorsByHexID.size(); }

  unsigned               numPixels() const { return fPixels.size(); }
  unsigned               numPixelHexSites() const { return fPixelsByHexID.size(); }

  inline const VSOObscuration* obscuration(unsigned id) const;

  inline const VSOMirror* mirror(unsigned id) const;
  inline const VSOMirror* mirrorByHexID(unsigned hexID) const;

  inline const VSOPixel* pixel(unsigned id) const;
  inline const VSOPixel* pixelByHexID(unsigned hexID) const;

 private:
  // ************************************************************************
  // Telescope Parameters
  // ************************************************************************
  unsigned        fID;                //!< Sequential ID (starting at 0)
  Eigen::Vector3d fPos;               //!< Position of the telescope
  double          fDeltaY;            //!< Nonorthogonality of rotation planes
  double          fAlphaX;            //!< Alignment angles(?)
  double          fAlphaY;            //!< Alignment angles(?)
  double          fElevation;         //!< Elevation of the telescope -- angle to the x-y plane
  double          fAzimuth;           //!< Azimuthal position -- angle made in the x-y plane with the y axis

  // ************************************************************************
  // Reflector Parameters
  // ************************************************************************
  Eigen::Vector3d fTranslation;       //!< Vector in reflector r.f. to the intersection of the rotation axes
  double          fCurvatureRadius;   //!< Radius of curvature of reflector
  double          fAperture;          //!< Diameter of projection of reflector onto a plane
  double          fFacetSpacing;      //!< Size of mirror
  double          fFacetSize;         //!< Spacing between mirrors
  double          fReflectorRotation; //!< Rotation about the axis of reflector
  double          fCosReflectorRotation; //!< Rotation about the axis of reflector
  double          fSinReflectorRotation; //!< Rotation about the axis of reflector
  unsigned        fHexagonRingsN;     //!< Number of hexagon rings of mirrors
  double          fReflectorIP;       //!< Diameter of sphere embedding reflector
  bool            fMirrorParity;      //!< Parity for counting mirrors; 1 => counting toward +x-axis
  //!  in reflector r.f.; -1 => counting toward -x-axis

  // ************************************************************************
  // Camera Parameters
  // ************************************************************************
  Eigen::Vector3d fFPTranslation;
  double          fCameraDiameter;    //!< Dependent on the field of view
  double          fFieldOfView;       //!< Field of view
  double          fCathodeDiameter;   //!< Diameter of photocathode
  double          fPixelSpacing;      //!< Spacing of pixels
  double          fPixelRotation;     //!< Rotation angle of pixels wrt grid [rad]
  double          fCosPixelRotation;  //!< Cos rotation angle of pixels wrt grid
  double          fSinPixelRotation;  //!< Sin rotation angle of pixels wrt grid

  double          fConcSurvProb;      //!< Pixel properties; probability of survival
  Eigen::Vector3d fFPRotation;
  double          fCameraIP;
  bool            fPixelParity;       //!< Parity for counting pixels; 1 => counting toward +x-axis
  //!  in reflector r.f.; -1 => counting toward -x-axis

  // ************************************************************************
  // Mirrors and Pixel data
  // ************************************************************************

  std::vector<VSOObscuration*>  fObscurations;
  std::vector<VSOMirror*>       fMirrors;
  std::vector<VSOMirror*>       fMirrorsByHexID;
  std::vector<VSOPixel*>        fPixels;
  std::vector<VSOPixel*>        fPixelsByHexID;

  // ************************************************************************
  // Precalculated rotation vector -- not stored in DB
  // ************************************************************************

  Eigen::Matrix3d fFPRotationMatrix;
  Eigen::Matrix3d fRotationVector;  //!<the pre-calculated rotation vector

  void calculateFPRotationMatrix();
  void calculateRotationVector();
};

inline const VSOObscuration* VSOTelescope::obscuration(unsigned id) const
{
  if(id>=fObscurations.size())return 0;
  else return fObscurations[id];
}

inline const VSOMirror* VSOTelescope::mirror(unsigned id) const
{
  if(id>=fMirrors.size())return 0;
  else return fMirrors[id];
}

inline const VSOMirror* VSOTelescope::mirrorByHexID(unsigned hexID) const
{
  if(hexID>=fMirrorsByHexID.size())return 0;
  else return fMirrorsByHexID[hexID];
}

inline const VSOPixel* VSOTelescope::pixel(unsigned id) const
{
  if(id>fPixels.size())return 0;
  else return fPixels[id];
}

inline const VSOPixel* VSOTelescope::pixelByHexID(unsigned hexID) const
{
  if(hexID>=fPixelsByHexID.size())return 0;
  else return fPixelsByHexID[hexID];
}

} } } // namespace calin::simulation::vs_optics
