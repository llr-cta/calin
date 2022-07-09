/*

   calin/simulation/vso_telescope.hpp -- Stephen Fegan -- 2015-11-09

   Class for telescope, a container for obscurations, mirrors and pixels

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

#include <iact_data/instrument_layout.pb.h>

namespace calin { namespace simulation { namespace vs_optics {

/*! \class VSOTelescope
  \brief VSOTelescope class;
*/

#ifndef SWIG
calin::ix::iact_data::instrument_layout::TelescopeLayout*
dc_parameters_to_telescope_layout(
  const ix::simulation::vs_optics::IsotropicDCArrayParameters& param,
  unsigned telescope_id = 0, const Eigen::Vector3d& pos = Eigen::Vector3d::Zero(),
  calin::ix::iact_data::instrument_layout::TelescopeLayout* d = nullptr);
#else
calin::ix::iact_data::instrument_layout::TelescopeLayout*
dc_parameters_to_telescope_layout(
  const ix::simulation::vs_optics::IsotropicDCArrayParameters& param,
  unsigned telescope_id = 0, const Eigen::Vector3d& pos = Eigen::Vector3d::Zero());
#endif

class VSOTelescope
{
public:
  VSOTelescope();
  VSOTelescope(unsigned TID, /*unsigned THID,*/ const Eigen::Vector3d&P,
               double DY, double AX, double AY, double EL, double AZ,
               const Eigen::Vector3d& T, double CR, double A,
               double FSP, double FS, double RR, double FGSX, double FGSY,
               unsigned HRN, double RIP, const Eigen::Vector3d& RIPC, bool MP,
               const Eigen::Vector3d& FPT, double CD, double FOV, double D,
               double PS, double PR, double PGSX, double PGSZ,
               double CSP, const Eigen::Vector3d& FPR,
               double CIP, bool PP,
               double WIN_FRONT, double WIN_RAD, double WIN_THICK,
               double WIN_N,
               const std::vector<VSOObscuration*>& OBSVEC_PRE = {},
               const std::vector<VSOObscuration*>& OBSVEC_POST = {},
               const std::vector<VSOObscuration*>& OBSVEC_CAM = {}
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

  void add_pre_reflection_obscuration(VSOObscuration* obs) {
    fPreObscurations.push_back(obs); }
  void add_post_reflection_obscuration(VSOObscuration* obs) {
    fPostObscurations.push_back(obs); }
  void add_camera_obscuration(VSOObscuration* obs) {
    fCameraObscurations.push_back(obs); }
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

  //! Transform particle global to reflector
  inline void globalToReflector_pos(Eigen::Vector3d& v) const;
  //! Transform particle reflector to global
  inline void reflectorToGlobal_pos(Eigen::Vector3d& v) const;
  //! Transform from focal plane to reflector
  inline void focalPlaneToReflector_pos(Eigen::Vector3d& v) const;
  //! Transform from reflector to focal plane
  inline void reflectorToFocalPlane_pos(Eigen::Vector3d& v) const;
  //! Transform directly from focal plane to global
  inline void focalPlaneToGlobal_pos(Eigen::Vector3d& v) const;

  //! Transform particle global to reflector
  inline void globalToReflector_mom(Eigen::Vector3d& v) const;
  //! Transform particle reflector to global
  inline void reflectorToGlobal_mom(Eigen::Vector3d& v) const;
  //! Transform from focal plane to reflector
  inline void focalPlaneToReflector_mom(Eigen::Vector3d& v) const;
  //! Transform from reflector to focal plane
  inline void reflectorToFocalPlane_mom(Eigen::Vector3d& v) const;
  //! Transform directly from focal plane to global
  inline void focalPlaneToGlobal_mom(Eigen::Vector3d& v) const;

  //! Transform particle global to reflector
  inline void globalToReflector(calin::math::ray::Ray& r) const;
  //! Transform particle reflector to global
  inline void reflectorToGlobal(calin::math::ray::Ray& r) const;
  //! Transform from focal plane to reflector
  inline void focalPlaneToReflector(calin::math::ray::Ray& r) const;
  //! Transform from reflector to focal plane
  inline void reflectorToFocalPlane(calin::math::ray::Ray& r) const;
  //! Transform directly from focal plane to global
  inline void focalPlaneToGlobal(calin::math::ray::Ray& r) const;

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
  // Convert to IACT data TelescopeLayout type as best we can
  // ************************************************************************

#ifndef SWIG
  calin::ix::iact_data::instrument_layout::TelescopeLayout*
  convert_to_telescope_layout(
    calin::ix::iact_data::instrument_layout::TelescopeLayout* d = nullptr) const;
#else
  calin::ix::iact_data::instrument_layout::TelescopeLayout* convert_to_telescope_layout() const;
  void convert_to_telescope_layout(calin::ix::iact_data::instrument_layout::TelescopeLayout* d) const;
#endif

  // ************************************************************************
  // Accessors
  // ************************************************************************

  unsigned               id() const { return fID; }
  double                 altitude() const { return fPos.z(); }
  const Eigen::Vector3d& pos() const { return fPos; }
  const Eigen::Vector3d& position() const { return fPos; }
  double                 fpOffset() const { return fFPOffset; }
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
  double                 facetGridShiftX() const { return fFacetGridShiftX; }
  double                 facetGridShiftZ() const { return fFacetGridShiftZ; }
  double                 cosReflectorRotation() const { return fCosReflectorRotation; }
  double                 sinReflectorRotation() const { return fSinReflectorRotation; }
  unsigned               mirrorHexRings() const { return fHexagonRingsN; }
  double                 reflectorIP() const { return fReflectorIP; }
  const Eigen::Vector3d& reflectorIPCenter() const { return fReflectorIPCenter; }
  bool                   mirrorParity() const { return fMirrorParity; }

  const Eigen::Vector3d& focalPlanePosition() const { return fFPTranslation; }
  double                 cameraDiameter() const { return fCameraDiameter; }
  double                 fov() const { return fFieldOfView; }
  double                 cathodeDiameter() const { return fCathodeDiameter; }
  double                 pixelSpacing() const { return fPixelSpacing; }
  double                 pixelRotation() const { return fPixelRotation; }
  double                 cosPixelRotation() const { return fCosPixelRotation; }
  double                 sinPixelRotation() const { return fSinPixelRotation; }
  double                 pixelGridShiftX() const { return fPixelGridShiftX; }
  double                 pixelGridShiftZ() const { return fPixelGridShiftZ; }
  double                 concentratorSurvivalProb() const { return fConcSurvProb; }
  const Eigen::Vector3d& focalPlaneRotation() const { return fFPRotation; }
  double                 cameraIP() const { return fCameraIP; }
  bool                   pixelParity() const { return fPixelParity; }

  double                 windowFront() const { return fWindowFront; }
  double                 windowOuterRadius() const { return fWindowOuterRadius; }
  double                 windowThickness() const { return fWindowThickness; }
  double                 windowRefractiveIndex() const { return fWindowRefractiveIndex; }

  unsigned               numPreReflectionObscurations() const { return fPreObscurations.size(); }
  unsigned               numPostReflectionObscurations() const { return fPostObscurations.size(); }
  unsigned               numCameraObscurations() const { return fCameraObscurations.size(); }

  unsigned               numMirrors() const { return fMirrors.size(); }
  unsigned               numMirrorHexSites() const { return fMirrorsByHexID.size(); }

  unsigned               numPixels() const { return fPixels.size(); }
  unsigned               numPixelHexSites() const { return fPixelsByHexID.size(); }

  inline const VSOObscuration* pre_reflection_obscuration(unsigned id) const;
  std::vector<VSOObscuration*> all_pre_reflection_obscurations() { return fPreObscurations; }
#ifndef SWIG
  std::vector<const VSOObscuration*> all_pre_reflection_obscurations() const {
    return std::vector<const VSOObscuration*>(fPreObscurations.begin(), fPreObscurations.end()); }
#endif

  inline const VSOObscuration* post_reflection_obscuration(unsigned id) const;
  std::vector<VSOObscuration*> all_post_reflection_obscurations() { return fPostObscurations; }
#ifndef SWIG
  std::vector<const VSOObscuration*> all_post_reflection_obscurations() const {
    return std::vector<const VSOObscuration*>(fPostObscurations.begin(), fPostObscurations.end()); }
#endif

  inline const VSOObscuration* camera_obscuration(unsigned id) const;
  std::vector<VSOObscuration*> all_camera_obscurations() { return fCameraObscurations; }
#ifndef SWIG
  std::vector<const VSOObscuration*> all_camera_obscurations() const {
    return std::vector<const VSOObscuration*>(fCameraObscurations.begin(), fCameraObscurations.end()); }
#endif

  inline const VSOMirror* mirror(unsigned id) const;
  inline const VSOMirror* mirrorByHexID(unsigned hexID) const;
  std::vector<VSOMirror*> all_mirrors() { return fMirrors; }

  inline const VSOPixel* pixel(unsigned id) const;
  inline const VSOPixel* pixelByHexID(unsigned hexID) const;
  std::vector<VSOPixel*> all_pixels() { return fPixels; }

  const Eigen::Matrix3d& rotationGlobalToReflector() const { return rot_global_to_reflector_; }
  const Eigen::Vector3d& translationGlobalToReflector() const { return off_global_to_reflector_; }

  const Eigen::Matrix3d& rotationReflectorToFP() const { return rot_reflector_to_camera_; }
  bool hasFPRotation() const { return fHasFPRotation; }

private:
  // ************************************************************************
  // Telescope Parameters
  // ************************************************************************
  unsigned        fID;                //!< Sequential ID (starting at 0)
  Eigen::Vector3d fPos;               //!< Position of the telescope
  double          fFPOffset;            //!< Nonorthogonality of rotation planes
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
  double          fFacetGridShiftX;   //!< X-shift of facet grid with respect to refector frame
  double          fFacetGridShiftZ;   //!< Z-shift of facet grid with respect to refector frame
  unsigned        fHexagonRingsN;     //!< Number of hexagon rings of mirrors
  double          fReflectorIP;       //!< Diameter of minimum sphere embedding reflector
  Eigen::Vector3d fReflectorIPCenter; //!< Center of minimum sphere embedding reflector in reflector r.f.
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
  double          fPixelGridShiftX;   //!< X-shift of pixel grid with respect to camera frame
  double          fPixelGridShiftZ;   //!< Z-shift of pixel grid with respect to camera frame

  double          fConcSurvProb;      //!< Pixel properties; probability of survival
  Eigen::Vector3d fFPRotation;
  double          fCameraIP;
  bool            fPixelParity;       //!< Parity for counting pixels; 1 => counting toward +x-axis
  //!  in reflector r.f.; -1 => counting toward -x-axis

  // ************************************************************************
  // Window Parameters
  // ************************************************************************

  double          fWindowFront;       //!< Front of outer window sphere
  double          fWindowOuterRadius; //!< Radius of window outer sphere
  double          fWindowThickness;   //!< Thickness of window
  double          fWindowRefractiveIndex; //!< Refractive index of window

  // ************************************************************************
  // Mirrors and Pixel data
  // ************************************************************************

  std::vector<VSOObscuration*>  fPreObscurations;
  std::vector<VSOObscuration*>  fPostObscurations;
  std::vector<VSOObscuration*>  fCameraObscurations;
  std::vector<VSOMirror*>       fMirrors;
  std::vector<VSOMirror*>       fMirrorsByHexID;
  std::vector<VSOPixel*>        fPixels;
  std::vector<VSOPixel*>        fPixelsByHexID;

  // ************************************************************************
  // Precalculated rotation vector -- not stored in DB
  // ************************************************************************

  bool fHasFPRotation = false;
  Eigen::Matrix3d rot_camera_to_reflector_;
  Eigen::Matrix3d rot_reflector_to_camera_;

  Eigen::Matrix3d rot_reflector_to_global_;  //!<the pre-calculated rotation vector
  Eigen::Matrix3d rot_global_to_reflector_;
  Eigen::Vector3d off_global_to_reflector_;

  Eigen::Matrix3d rot_camera_to_global_;
  Eigen::Vector3d off_global_to_camera_;

  void calculateFPRotationMatrix();
  void calculateRotationVector();
};

inline const VSOObscuration* VSOTelescope::pre_reflection_obscuration(unsigned id) const
{
  if(id>=fPreObscurations.size())return 0;
  else return fPreObscurations[id];
}

inline const VSOObscuration* VSOTelescope::post_reflection_obscuration(unsigned id) const
{
  if(id>=fPostObscurations.size())return 0;
  else return fPostObscurations[id];
}

inline const VSOObscuration* VSOTelescope::camera_obscuration(unsigned id) const
{
  if(id>=fCameraObscurations.size())return 0;
  else return fCameraObscurations[id];
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

//! Transform particle global to reflector
inline void VSOTelescope::globalToReflector_pos(Eigen::Vector3d& v) const
{
  // First: Translate from global directly to reflector
  v -= off_global_to_reflector_;
  // Second: Rotate coordinate system to reflector orientation
  v = rot_global_to_reflector_ * v;
}

//! Transform particle reflector to global
inline void VSOTelescope::reflectorToGlobal_pos(Eigen::Vector3d& v) const
{
  // First: Rotate coordinate system to ground based
  v = rot_reflector_to_global_ * v;
  // Second: Translate from drive axes intersection to center of array
  v += off_global_to_reflector_;
}

//! Transform from focal plane to reflector
inline void VSOTelescope::focalPlaneToReflector_pos(Eigen::Vector3d& v) const
{
  // First: Rotate coordinate system
  if(fHasFPRotation) v = rot_camera_to_reflector_ * v;
  // Second: Translate from center of Focal Plane
  v += fFPTranslation;
}

//! Transform from reflector to focal plane
inline void VSOTelescope::reflectorToFocalPlane_pos(Eigen::Vector3d& v) const
{
  // First: Translate to center of Focal Plane
  v -= fFPTranslation;
  // Second: Rotate coordinate system
  if(fHasFPRotation) v = rot_reflector_to_camera_ * v;
}

//! Transform directly from focal plane to global
inline void VSOTelescope::focalPlaneToGlobal_pos(Eigen::Vector3d& v) const
{
  v = rot_camera_to_global_ * v;
  v += off_global_to_camera_;
}

//! Transform particle global to reflector
inline void VSOTelescope::globalToReflector_mom(Eigen::Vector3d& v) const
{
  v = rot_global_to_reflector_ * v;
}

//! Transform particle reflector to global
inline void VSOTelescope::reflectorToGlobal_mom(Eigen::Vector3d& v) const
{
  v = rot_reflector_to_global_ * v;
}

//! Transform from focal plane to reflector
inline void VSOTelescope::focalPlaneToReflector_mom(Eigen::Vector3d& v) const
{
  if(fHasFPRotation) v = rot_camera_to_reflector_ * v;
}

//! Transform from reflector to focal plane
inline void VSOTelescope::reflectorToFocalPlane_mom(Eigen::Vector3d& v) const
{
  if(fHasFPRotation) v = rot_reflector_to_camera_ * v;
}

//! Transform directly from focal plane to global
inline void VSOTelescope::focalPlaneToGlobal_mom(Eigen::Vector3d& v) const
{
  v = rot_camera_to_global_ * v;
}

inline void VSOTelescope::globalToReflector(calin::math::ray::Ray& r) const
{
  // First: Translate from global directly to reflector
  r.translate_origin(off_global_to_reflector_);
  // Second: Rotate coordinate system to reflector orientation
  r.rotate(rot_global_to_reflector_);
}

inline void VSOTelescope::reflectorToGlobal(calin::math::ray::Ray& r) const
{
  // First: Rotate coordinate system to ground based
  r.rotate(rot_reflector_to_global_);
  // Second: Translate from drive axes intersection to center of array
  r.untranslate_origin(off_global_to_reflector_);
}

inline void VSOTelescope::focalPlaneToReflector(calin::math::ray::Ray& r) const
{
  // First: Rotate coordinate system
  if(fHasFPRotation) r.rotate(rot_camera_to_reflector_);
  // Second: Translate from center of Focal Plane
  r.untranslate_origin(fFPTranslation);
}

inline void VSOTelescope::reflectorToFocalPlane(calin::math::ray::Ray& r) const
{
  // First: Translate to center of Focal Plane
  r.translate_origin(fFPTranslation);
  // Second: Rotate coordinate system
  if(fHasFPRotation)r.rotate(rot_reflector_to_camera_);
}

//! Transform directly from focal plane to global
inline void VSOTelescope::focalPlaneToGlobal(calin::math::ray::Ray& r) const
{
  r.rotate(rot_camera_to_global_);
  r.untranslate_origin(off_global_to_camera_);
}

} } } // namespace calin::simulation::vs_optics
