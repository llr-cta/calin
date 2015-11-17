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

#include <math/vs_vec3d.hpp>
#include <math/vs_particle.hpp>
//#include <RandomNumbers.hpp>
//#include <xytohex.hpp>

//#include "VSOArrayParameters.hpp"
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
  VSOTelescope(unsigned TID, unsigned THID, const math::vs_physics::Vec3D&P, 
               double DY, double AX, double AY, double EL, double AZ,
               const math::vs_physics::Vec3D& T, double CR, double A, 
               double FSP, double FS, double RR,
               unsigned HRN, double RIP, bool MP, 
               const math::vs_physics::Vec3D& FPT, double CD, double FOV, double D,
               double PS, double CSP, const math::vs_physics::Vec3D& FPR, 
               double CIP, bool PP,
#if 0
               bool SEC, double RI, double RI1, double RI2, double C10, 
               double C11, double C12, double C13, double C14, double C20,
               double C21, double C22,double C23, double C24, 
#endif
               const std::vector<VSOObscuration*>& OBSVEC = {}
               );
  VSOTelescope(const VSOTelescope& o);
  virtual ~VSOTelescope();
    
  const VSOTelescope& operator =(const VSOTelescope& o);
    
  // ************************************************************************
  // Create a new telescope randomly using parameters
  // ************************************************************************

#if 0
  void populateMirrorsAndPixelsRandom(const VSOArrayParameters&,
                                      RandomNumbers& rng);
#endif

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
  bool pointTelescope(const math::vs_physics::Vec3D& v); 
  bool pointTelescopeAzEl(const double az_rad, const double el_rad);
  
  //! Transform position-like vector from global to reflector
  void globalToReflector_pos(math::vs_physics::Vec3D& v) const; 
  //! Transform momentum-like vector from global to reflector
  void globalToReflector_mom(math::vs_physics::Vec3D& v) const; 
  //! Transform position-like vector from reflector to global
  void reflectorToGlobal_pos(math::vs_physics::Vec3D& v) const; 
  //! Transform momentum-like vector from reflector to global
  void reflectorToGlobal_mom(math::vs_physics::Vec3D& v) const; 

  //! Transform particle global to reflector
  void globalToReflector(math::vs_physics::Particle& p) const; 
  //! Transform particle reflector to global
  void reflectorToGlobal(math::vs_physics::Particle& p) const; 
  //! Transform from focal plane to reflector
  void focalPlaneToReflector(math::vs_physics::Particle& p) const; 
  //! Transform from reflector to focal plane
  void reflectorToFocalPlane(math::vs_physics::Particle& p) const; 
  
  // ************************************************************************
  // Dump and Reload
  // ************************************************************************

  void dump_to_proto(ix::simulation::vs_optics::VSOTelescopeData* d) const;
  static VSOTelescope*
  create_from_proto(const ix::simulation::vs_optics::VSOTelescopeData& d);
  
#if 0
  void dumpShort(std::ostream& stream) const;
  void dump(std::ostream& stream, unsigned l=0) const;
  static VSOTelescope* createFromShortDump(std::istream& stream);
#endif
  // ************************************************************************
  // Accessors
  // ************************************************************************

  unsigned       id() const { return fID; }
  unsigned       hexID() const { return fTelescopeHexID; }
  double         altitude() const { return fPos.z; }
  const math::vs_physics::Vec3D& pos() const { return fPos; }
  const math::vs_physics::Vec3D& position() const { return fPos; }
  double         deltaY() const { return fDeltaY; }
  double         alphaX() const { return fAlphaX; } 
  double         alphaY() const { return fAlphaY; }
  double         azimuth() const { return fAzimuth; }
  double         elevation() const { return fElevation; }
  math::vs_physics::Vec3D opticalAxis() const;
    
  const math::vs_physics::Vec3D& translation() const { return fTranslation; }
  double         curvatureRadius() const { return fCurvatureRadius; }
  double         aperture() const { return fAperture; }
  double         facetSpacing() const { return  fFacetSpacing; }
  double         facetSize() const { return  fFacetSize; }
  double         reflectorRotation() const { return fReflectorRotation; }
  unsigned       mirrorHexRings() const { return fHexagonRingsN; }
  double         reflectorIP() const { return fReflectorIP; }
  bool           mirrorParity() const { return fMirrorParity; }
    
  const math::vs_physics::Vec3D& focalPlanePosition() const { return fFPTranslation; }
  double         cameraDiameter() const { return fCameraDiameter; }
  double         fov() const { return fFieldOfView; }
  double         cathodeDiameter() const { return fCathodeDiameter; }
  double         pixelSpacing() const { return fPixelSpacing; }
  double         concentratorSurvivalProb() const { return fConcSurvProb; }
  const math::vs_physics::Vec3D& focalPlaneRotion() const { return fFPRotation; }
  double         cameraIP() const { return fCameraIP; }
  bool           pixelParity() const { return fPixelParity; }
    
#if 0
  bool           hasSecondary() const { return fHasSecondary; }
  double         refractiveIndex() const { return fRefractiveIndex; }
  double         refractiveIndex1() const { return fRefractiveIndex1; }
  double         refractiveIndex2() const { return fRefractiveIndex2; }
  double         ce1Parameter0() const { return fCE1Parameter0; }
  double         ce1Parameter2() const { return fCE1Parameter2; }
  double         ce1Parameter3() const { return fCE1Parameter3; }
  double         ce1Parameter4() const { return fCE1Parameter4; }
  double         ce1Parameter5() const { return fCE1Parameter5; }
  double         ce2Parameter0() const { return fCE2Parameter0; }
  double         ce2Parameter2() const { return fCE2Parameter2; }
  double         ce2Parameter3() const { return fCE2Parameter3; }
  double         ce2Parameter4() const { return fCE2Parameter4; }
  double         ce2Parameter5() const { return fCE2Parameter5; }
#endif

  unsigned       numObscurations() const { return fObscurations.size(); }

  unsigned       numMirrors() const { return fMirrors.size(); }
  unsigned       numMirrorHexSites() const { return fMirrorsByHexID.size(); }
    
  unsigned       numPixels() const { return fPixels.size(); }
  unsigned       numPixelHexSites() const { return fPixelsByHexID.size(); }
    
  inline const VSOObscuration* obscuration(unsigned id) const;    

  inline const VSOMirror* mirror(unsigned id) const;
  inline const VSOMirror* mirrorByHexID(unsigned hexID) const;
    
  inline const VSOPixel* pixel(unsigned id) const;
  inline const VSOPixel* pixelByHexID(unsigned hexID) const;
    
 private:
  // ************************************************************************
  // Telescope Parameters
  // ************************************************************************
  unsigned       fID;                //!< Sequential ID (starting at 0)
  unsigned       fTelescopeHexID;    //!< Hex ID on the honeycomb (from 1)
  math::vs_physics::Vec3D fPos;               //!< Position of the telescope
  double         fDeltaY;            //!< Nonorthogonality of rotation planes
  double         fAlphaX;            //!< Alignment angles(?)
  double         fAlphaY;            //!< Alignment angles(?)
  double         fElevation;         //!< Elevation of the telescope -- angle to the x-y plane
  double         fAzimuth;           //!< Azimuthal position -- angle made in the x-y plane with the y axis
  
  // ************************************************************************
  // Reflector Parameters
  // ************************************************************************
  math::vs_physics::Vec3D fTranslation;       //!< Vector in reflector r.f. to the intersection of the rotation axes
  double         fCurvatureRadius;   //!< Radius of curvature of reflector
  double         fAperture;          //!< Diameter of projection of reflector onto a plane
  double         fFacetSpacing;      //!< Size of mirror
  double         fFacetSize;         //!< Spacing between mirrors
  double         fReflectorRotation; //!< Rotation about the axis of reflector
  unsigned       fHexagonRingsN;     //!< Number of hexagon rings of mirrors
  double         fReflectorIP;       //!< Diameter of sphere embedding reflector 
  bool           fMirrorParity;      //!< Parity for counting mirrors; 1 => counting toward +x-axis
  //!  in reflector r.f.; -1 => counting toward -x-axis
    
  // ************************************************************************
  // Camera Parameters
  // ************************************************************************
  math::vs_physics::Vec3D fFPTranslation;     
  double         fCameraDiameter;    //!< Dependent on the field of view
  double         fFieldOfView;       //!< Field of view
  double         fCathodeDiameter;   //!< Diameter of photocathode
  double         fPixelSpacing;      //!< Spacing of pixels
  double         fConcSurvProb;      //!< Pixel properties; probability of survival
  math::vs_physics::Vec3D fFPRotation; 
  double         fCameraIP;
  bool           fPixelParity;       //!< Parity for counting pixels; 1 => counting toward +x-axis
  //!  in reflector r.f.; -1 => counting toward -x-axis

#if 0    
  // ************************************************************************
  // Secondary Optics Parameters
  // ************************************************************************
  bool           fHasSecondary;
  double         fRefractiveIndex;   //!<secondary optics parameters
  double         fRefractiveIndex1;  //!<secondary optics parameters
  double         fRefractiveIndex2;  //!<secondary optics parameters
  double         fCE1Parameter0;
  double         fCE1Parameter2;
  double         fCE1Parameter3;
  double         fCE1Parameter4;
  double         fCE1Parameter5;
  double         fCE2Parameter0;
  double         fCE2Parameter2;
  double         fCE2Parameter3;
  double         fCE2Parameter4;
  double         fCE2Parameter5;
#endif

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
    
  math::vs_physics::Vec3D fRotationVector;  //!<the pre-calculated rotation vector
    
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
  if((hexID<1)||(hexID>fMirrorsByHexID.size()))return 0;
  else return fMirrorsByHexID[hexID-1];
}
  
inline const VSOPixel* VSOTelescope::pixel(unsigned id) const
{
  if(id>fPixels.size())return 0;
  else return fPixels[id];
}
  
inline const VSOPixel* VSOTelescope::pixelByHexID(unsigned hexID) const
{
  if((hexID<1)||(hexID>fPixelsByHexID.size()))return 0;
  else return fPixelsByHexID[hexID-1];
}

} } } // namespace calin::simulation::vs_optics

