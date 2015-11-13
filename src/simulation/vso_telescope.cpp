/* 

   calin/simulation/vso_telescope.hpp -- Stephen Fegan -- 2015-11-11

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

#include <sstream>
#include <algorithm>

#include <simulation/vso_telescope.hpp>
#include <proto/simulation_vs_optics.pb.h>

using namespace calin::math::vs_physics;
using namespace calin::simulation::vs_optics;

VSOTelescope::VSOTelescope():
  fID(), fTelescopeHexID(), fPos(), 
  fDeltaY(), fAlphaX(), fAlphaY(),
  fElevation(), fAzimuth(), 
  fTranslation(), fCurvatureRadius(), fAperture(), 
  fFacetSpacing(), fFacetSize(), 
  fReflectorRotation(), fHexagonRingsN(),
  fReflectorIP(), fMirrorParity(), fFPTranslation(), fCameraDiameter(),
  fFieldOfView(), fCathodeDiameter(), fPixelSpacing(), fConcSurvProb(),
  fFPRotation(), fCameraIP(), fPixelParity(),
#if 0
  fHasSecondary(false),
  fRefractiveIndex(), fRefractiveIndex1(), fRefractiveIndex2(), 
  fCE1Parameter0(), fCE1Parameter2(), fCE1Parameter3(), fCE1Parameter4(), 
  fCE1Parameter5(), fCE2Parameter0(), fCE2Parameter2(), fCE2Parameter3(),
  fCE2Parameter4(), fCE2Parameter5(), 
#endif
  fObscurations(),
  fMirrors(), fMirrorsByHexID(), fPixels(), fPixelsByHexID(), fRotationVector()
{
  calculateRotationVector();
}

VSOTelescope::
VSOTelescope(unsigned TID, unsigned THID, const Vec3D&P, 
	     double DY, double AX, double AY, double EL, double AZ,
	     const Vec3D& T, double CR, double A, double FSP, double FS, 
	     double RR, unsigned HRN, double RIP, bool MP, 
	     const Vec3D& FPT, double CD, double FOV, double D, double PS, 
	     double CSP, const Vec3D& FPR, double CIP, bool PP,
#if 0
	     bool SEC, double RI, double RI1, double RI2, double C10, 
	     double C12, double C13, double C14, double C15, double C20,
	     double C22, double C23,double C24, double C25, 
#endif
	     const std::vector<VSOObscuration*>& OBSVEC
	     ):
  fID(TID), fTelescopeHexID(THID), fPos(P),
  fDeltaY(DY), fAlphaX(AX), fAlphaY(AY), fElevation(EL), fAzimuth(AZ), 
  fTranslation(T), fCurvatureRadius(CR), fAperture(A), fFacetSpacing(FSP), 
  fFacetSize(FS), fReflectorRotation(RR), 
  fHexagonRingsN(HRN), fReflectorIP(RIP), fMirrorParity(MP), 
  fFPTranslation(FPT),  fCameraDiameter(CD), fFieldOfView(FOV), 
  fCathodeDiameter(D), fPixelSpacing(PS), fConcSurvProb(CSP),
  fFPRotation(FPR), fCameraIP(CIP), fPixelParity(PP),
#if 0
  fHasSecondary(SEC),
  fRefractiveIndex(RI), fRefractiveIndex1(RI1), fRefractiveIndex2(RI2),
  fCE1Parameter0(C10), fCE1Parameter2(C12), fCE1Parameter3(C13), 
  fCE1Parameter4(C14), fCE1Parameter5(C15), fCE2Parameter0(C20), 
  fCE2Parameter2(C22), fCE2Parameter3(C23), fCE2Parameter4(C24),
  fCE2Parameter5(C25), 
#endif
  fObscurations(OBSVEC),
  fMirrors(), fMirrorsByHexID(), fPixels(), fPixelsByHexID(), fRotationVector()
{
  calculateRotationVector();
}

VSOTelescope::VSOTelescope(const VSOTelescope& o):
  fID(o.fID), fTelescopeHexID(o.fTelescopeHexID),
  fPos(o.fPos), fDeltaY(o.fDeltaY), fAlphaX(o.fAlphaX), fAlphaY(o.fAlphaY),
  fElevation(o.fElevation), fAzimuth(o.fAzimuth), fTranslation(o.fTranslation),
  fCurvatureRadius(o.fCurvatureRadius), fAperture(o.fAperture), 
  fFacetSpacing(o.fFacetSpacing), fFacetSize(o.fFacetSize), 
  fReflectorRotation(o.fReflectorRotation),
  fHexagonRingsN(o.fHexagonRingsN),
  fReflectorIP(o.fReflectorIP), fMirrorParity(o.fMirrorParity),
  fFPTranslation(o.fFPTranslation), fCameraDiameter(o.fCameraDiameter),
  fFieldOfView(o.fFieldOfView), fCathodeDiameter(o.fCathodeDiameter),
  fPixelSpacing(o.fPixelSpacing), fConcSurvProb(o.fConcSurvProb),
  fFPRotation(o.fFPRotation), fCameraIP(o.fCameraIP),
  fPixelParity(o.fPixelParity), 
#if 0
  fHasSecondary(o.fHasSecondary),
  fRefractiveIndex(o.fRefractiveIndex),
  fRefractiveIndex1(o.fRefractiveIndex1), 
  fRefractiveIndex2(o.fRefractiveIndex2),
  fCE1Parameter0(o.fCE1Parameter0), fCE1Parameter2(o.fCE1Parameter2),
  fCE1Parameter3(o.fCE1Parameter3), fCE1Parameter4(o.fCE1Parameter4),
  fCE1Parameter5(o.fCE1Parameter5), fCE2Parameter0(o.fCE2Parameter0),
  fCE2Parameter2(o.fCE2Parameter2), fCE2Parameter3(o.fCE2Parameter3),
  fCE2Parameter4(o.fCE2Parameter4), fCE2Parameter5(o.fCE2Parameter5),
#endif
  fObscurations(),
  fMirrors(), fMirrorsByHexID(), fPixels(), fPixelsByHexID(),
  fRotationVector()
{
  fMirrors.resize(o.fMirrors.size());
  fMirrorsByHexID.resize(o.fMirrorsByHexID.size());
  for(std::vector<VSOMirror*>::const_iterator i=o.fMirrors.begin(); 
      i!=o.fMirrors.end(); i++)
    {
      VSOMirror* mirror = new VSOMirror(**i);
      fMirrors[(*i)->id()]=mirror;
      fMirrorsByHexID[(*i)->hexID()]=mirror;
    }

  fPixels.resize(o.fPixels.size());
  fPixelsByHexID.resize(o.fPixelsByHexID.size());
  for(std::vector<VSOPixel*>::const_iterator i=o.fPixels.begin(); 
      i!=o.fPixels.end(); i++)
    {
      VSOPixel* pixel = new VSOPixel(**i);
      fPixels[(*i)->id()]=pixel;
      fPixelsByHexID[(*i)->hexID()]=pixel;
    }

  for(std::vector<VSOObscuration*>::const_iterator i=o.fObscurations.begin();
      i!=o.fObscurations.end(); i++)fObscurations.push_back((*i)->clone());

  calculateRotationVector();
}

VSOTelescope::~VSOTelescope()
{
  for(std::vector<VSOMirror*>::iterator i=fMirrors.begin();
      i!=fMirrors.end(); i++)delete *i;
  for(std::vector<VSOPixel*>::iterator i=fPixels.begin();
      i!=fPixels.end(); i++)delete *i;
  for(std::vector<VSOObscuration*>::iterator i=fObscurations.begin();
      i!=fObscurations.end(); i++)delete *i;
}

// Stroustrup third edition sec 11.3.4 recommends default copy assignment
const VSOTelescope& VSOTelescope::operator =(const VSOTelescope& o)
{
  fID                = o.fID;
  fTelescopeHexID    = o.fTelescopeHexID;
  fPos               = o.fPos;
  fDeltaY            = o.fDeltaY;
  fAlphaX            = o.fAlphaX;
  fAlphaY            = o.fAlphaY;
  fElevation         = o.fElevation;
  fAzimuth           = o.fAzimuth;
  fTranslation       = o.fTranslation;
  fCurvatureRadius   = o.fCurvatureRadius;
  fAperture          = o.fAperture;
  fFacetSpacing      = o.fFacetSpacing;
  fFacetSize         = o.fFacetSize;
  fReflectorRotation = o.fReflectorRotation;
  fHexagonRingsN     = o.fHexagonRingsN;
  fReflectorIP       = o.fReflectorIP;
  fMirrorParity      = o.fMirrorParity;
  fFPTranslation     = o.fFPTranslation;
  fCameraDiameter    = o.fCameraDiameter;
  fFieldOfView       = o.fFieldOfView;
  fCathodeDiameter   = o.fCathodeDiameter;
  fPixelSpacing      = o.fPixelSpacing;
  fConcSurvProb      = o.fConcSurvProb;
  fFPRotation        = o.fFPRotation;
  fCameraIP          = o.fCameraIP;
  fPixelParity       = o.fPixelParity;
#if 0
  fHasSecondary      = o.fHasSecondary,
  fRefractiveIndex   = o.fRefractiveIndex;
  fRefractiveIndex1  = o.fRefractiveIndex1; 
  fRefractiveIndex2  = o.fRefractiveIndex2;
  fCE1Parameter0     = o.fCE1Parameter0; 
  fCE1Parameter2     = o.fCE1Parameter2;
  fCE1Parameter3     = o.fCE1Parameter3; 
  fCE1Parameter4     = o.fCE1Parameter4;
  fCE1Parameter5     = o.fCE1Parameter5;
  fCE2Parameter0     = o.fCE2Parameter0;
  fCE2Parameter2     = o.fCE2Parameter2;
  fCE2Parameter3     = o.fCE2Parameter3;
  fCE2Parameter4     = o.fCE2Parameter4;
  fCE2Parameter5     = o.fCE2Parameter5;
#endif

  for(std::vector<VSOMirror*>::iterator i=fMirrors.begin();
      i!=fMirrors.end(); i++)delete *i;

  fMirrors.clear();
  fMirrorsByHexID.clear();

  fMirrors.resize(o.fMirrors.size());
  fMirrorsByHexID.resize(o.fMirrorsByHexID.size());

  for(std::vector<VSOMirror*>::const_iterator i=o.fMirrors.begin(); 
      i!=o.fMirrors.end(); i++)
    {
      VSOMirror* mirror = new VSOMirror(**i);
      fMirrors[(*i)->id()]=mirror;
      fMirrorsByHexID[(*i)->hexID()]=mirror;
    }

  for(std::vector<VSOPixel*>::iterator i=fPixels.begin();
      i!=fPixels.end(); i++)delete *i;

  fPixels.clear();
  fPixelsByHexID.clear();

  fPixels.resize(o.fPixels.size());
  fPixelsByHexID.resize(o.fPixelsByHexID.size());

  for(std::vector<VSOPixel*>::const_iterator i=o.fPixels.begin(); 
      i!=o.fPixels.end(); i++)
    {
      VSOPixel* pixel = new VSOPixel(**i);
      fPixels[(*i)->id()]=pixel;
      fPixelsByHexID[(*i)->hexID()]=pixel;
    }

  calculateRotationVector();

  return *this;
}

// ****************************************************************************
// Accessor
// ****************************************************************************

Vec3D VSOTelescope::opticalAxis() const
{
  return Vec3D(cos(fElevation)*sin(fAzimuth),
	       cos(fElevation)*cos(fAzimuth),
	       sin(fElevation));
}

// ****************************************************************************
// Repoint the telescope along a vector
// ****************************************************************************

bool VSOTelescope::pointTelescope(const Vec3D& v)
{
  if(v.Norm2()==0)return false;
  fElevation = atan2(v.z,sqrt(v.x*v.x + v.y*v.y));
  fAzimuth = fmod(atan2(v.x,v.y)+2.0*M_PI, 2.0*M_PI);
  calculateRotationVector();
  return true;
}

bool VSOTelescope::pointTelescopeAzEl(const double az_rad, const double el_rad)
{
  fElevation = fmod(fmod(el_rad,2.0*M_PI)+2.0*M_PI, 2.0*M_PI);
  fAzimuth = fmod(fmod(az_rad,2.0*M_PI)+2.0*M_PI, 2.0*M_PI);
  calculateRotationVector();
  return true;
}

// ****************************************************************************
// Calculate rotation vector and map between Global and Telescope coordinates
// ***************************************************************************

void VSOTelescope::calculateRotationVector()
{
  // Rotation vector maps from Reflector to Global
  fRotationVector = 
    Vec3D(1,0,0)*fElevation &
    Vec3D(0,1,0)*fDeltaY &
    Vec3D(0,0,-1)*fAzimuth &
    Vec3D(0,1,0)*fAlphaX &
    Vec3D(1,0,0)*fAlphaY;
}

void VSOTelescope::globalToReflector_pos(Vec3D& v) const
{
  // First: Translate from center of array to drive axes intersection
  v -= fPos;
  // Second: Rotate coordinate system to reflector orientation
  v.Rotate(-fRotationVector);
  // Third: Translate from intersection of drive axes to reflector
  v += fTranslation;
}

void VSOTelescope::globalToReflector_mom(Vec3D& v) const
{
  // Rotate coordinate system to reflector orientation
  v.Rotate(-fRotationVector);
}

void VSOTelescope::reflectorToGlobal_pos(Vec3D& v) const
{
  // First: Translate from reflector to intersection of drive axes
  v -= fTranslation;
  // Second: Rotate coordinate system to ground based
  v.Rotate(fRotationVector);
  // Third: Translate from drive axes intersection to center of array
  v += fPos;
}

void VSOTelescope::reflectorToGlobal_mom(Vec3D& v) const
{
  // Rotate coordinate system to ground based
  v.Rotate(fRotationVector);
}

void VSOTelescope::globalToReflector(Particle& p) const
{
  // First: Translate from center of array to drive axes intersection
  p.TranslateOrigin(Vec4D(0,fPos));
  // Second: Rotate coordinate system to reflector orientation
  p.Rotate(-fRotationVector);
  // Third: Translate from intersection of drive axes to reflector
  p.TranslateOrigin(Vec4D(0,-fTranslation));
}

void VSOTelescope::reflectorToGlobal(Particle& p) const
{
  // First: Translate from reflector to intersection of drive axes
  p.TranslateOrigin(Vec4D(0,fTranslation));
  // Second: Rotate coordinate system to ground based
  p.Rotate(fRotationVector);
  // Third: Translate from drive axes intersection to center of array
  p.TranslateOrigin(Vec4D(0,-fPos));
}

void VSOTelescope::focalPlaneToReflector(Particle& p) const
{
  // First: Rotate coordinate system
  p.Rotate(fFPRotation);
  // Second: Translate from center of Focal Plane
  p.TranslateOrigin(Vec4D(0,-fFPTranslation));
}

void VSOTelescope::reflectorToFocalPlane(Particle& p) const
{
  // First: Translate to center of Focal Plane
  p.TranslateOrigin(Vec4D(0,fFPTranslation));
  // Second: Rotate coordinate system
  p.Rotate(-fFPRotation);
}

// ****************************************************************************
// Fill the mirrors and pixels tables up with randomly generated objects
// ****************************************************************************

#if 0
void VSOTelescope::
populateMirrorsAndPixelsRandom(const VSOArrayParameters& param, 
			       RandomNumbers& rng)
{
  // **************************************************************************
  // Clear the MIRRORs table and repopulate it with randomly generated mirrors
  // **************************************************************************
  double reflector_r2 = fAperture*fAperture/4.0;
  double reflector_c2 = fCurvatureRadius*fCurvatureRadius;
  Vec3D reflector_center(0, fCurvatureRadius, 0);

  std::set<unsigned> mirrors_missing;
  tokenize(param.MirrorMissingList,mirrors_missing);

  for(std::vector<VSOMirror*>::iterator i=fMirrors.begin();
      i!=fMirrors.end(); i++)delete *i;
  fMirrors.clear();
  fMirrorsByHexID.clear();  

  int num_hex_mirror_sites = 3*fHexagonRingsN*(fHexagonRingsN+1)+1;

  unsigned id = 0;
  for(int i=0; i<num_hex_mirror_sites; i++)
    {
      int hexid = i+1;

      if(mirrors_missing.find(hexid) != mirrors_missing.end())
	{
	  fMirrorsByHexID.push_back(0);
	  continue; // skip mirror if on the mirring list
	}

      Vec3D nominal_position;

      // Compute the mirror's nominal position
      nh_to_xy(&hexid, &nominal_position.x, &nominal_position.z); 
      if(fMirrorParity)nominal_position.x=-nominal_position.x;
      nominal_position.x *= fFacetSpacing;
      nominal_position.z *= fFacetSpacing;

      double X2 = nominal_position.x*nominal_position.x;
      double Z2 = nominal_position.z*nominal_position.z;

      if( (reflector_r2 - X2 - Z2) <= 0 )
	{
	  fMirrorsByHexID.push_back(0);
	  continue; // skip mirror if projected hex position too far out
	}

      nominal_position.y = fCurvatureRadius-sqrt(reflector_c2 - X2 - Z2);

      // Add Gaussian normal/tangenetial position error
      Vec3D reflector_normal(reflector_center-nominal_position);
      reflector_normal /= reflector_normal.Norm();

      Vec3D position(nominal_position);
      position += reflector_normal*(rng.Normal()*param.MirrorPosNormalDisp);
      position -= reflector_center;
      position.ScatterDirection(param.MirrorPosTangentDisp/fCurvatureRadius,rng);
      position += reflector_center;
      
      // Rotate by global rotation angle
      position.Rotate(Vec3D(0,fReflectorRotation,0));

      // Get the (perturbed) alignment angle of mirror
      Vec3D alignment;

      if(param.ReflectorAlignMode == 0)
	{
	  Vec3D alignment_pt(param.ReflectorAlignPtX,
			     param.ReflectorAlignPtY,
			     param.ReflectorAlignPtY);
	  if(alignment_pt.y == 0)
	    {
	      double F = param.ReflectorAlignImagePlane;
	      if(F==0)F = fFPTranslation.y;
	      double d = (Vec3D(0,F,0)-nominal_position).Norm();
	      alignment_pt.y = F+d;	      
	    }

	  // Standard DC alignment to a fixed point in space (with scatter)
	  alignment = (alignment_pt-position);
	  alignment /= alignment.Norm();

	  double align_disp = 
	    param.MirrorAlignTangentDisp/alignment_pt.Norm();
	  alignment.ScatterDirection(align_disp,rng);
	}
      else
	{
	  // Alignment to a point in the focal plane
	  
	  double stheta = sin(param.ReflectorFPAlignTheta);
	  double ctheta = cos(param.ReflectorFPAlignTheta);
	  double ttheta = stheta/ctheta;

	  double sphi = sin(param.ReflectorFPAlignPhi);
	  double cphi = cos(param.ReflectorFPAlignPhi);

	  double y_fp = param.ReflectorAlignImagePlane; // Image plane
	  if(y_fp==0)y_fp = fFPTranslation.y;
	  Vec3D r_fp(y_fp*ttheta*sphi ,y_fp, -y_fp*ttheta*cphi);

	  Vec3D e_in(-stheta*sphi, ctheta, stheta*cphi);
	  if((param.ReflectorAlignObjectPlane!=0) &&
	     (!std::isinf(param.ReflectorAlignObjectPlane)))
	    {
	      double y_em = param.ReflectorAlignObjectPlane;
	      Vec3D emission_pt(-y_em*ttheta*sphi, y_em, y_em*ttheta*cphi);
	      e_in = emission_pt-position;
	      e_in /= e_in.Norm();
	    }

	  Vec3D e_out = r_fp-position;
	  e_out /= e_out.Norm();

	  alignment = e_in;

	  Vec3D e_rot = e_in^e_out;
	  double strot2 = e_rot.Norm();
	  if(strot2 != 0)
	    {
	      double ctrot2 = e_in*e_out;
	      double trot2 = atan2(strot2,ctrot2);
	      e_rot *= 0.5*trot2/strot2;
	      alignment.Rotate(e_rot);
	    }
	}

      double focal_length = 
	param.MirrorFLength + rng.Normal()*param.MirrorFLengthDisp;

      double spot_size;
      if(param.MirrorSpotSizeDisp > 0)
	spot_size = rng.GammaByMeanAndStdDev(param.MirrorSpotSize,
					     param.MirrorSpotSizeDisp);
      else spot_size = param.MirrorSpotSize;

      // If param.MirrorSpotSizePhotonFraction not set -- assume FWHM
      if((param.MirrorSpotSizePhotonFraction>0.0)&&
	 (param.MirrorSpotSizePhotonFraction<1.0))
	spot_size /= 
	  2.0*sqrt(log(1.0/(1.0-param.MirrorSpotSizePhotonFraction)));
      else spot_size /= 2.0*sqrt(log(2.0));
      
      VSOMirror* mirror = 
	new VSOMirror(this, id, hexid, false, position, alignment, 
		     focal_length, spot_size, param.MirrorDegradingFactor);
      
      fMirrors.push_back(mirror);
      fMirrorsByHexID.push_back(mirror);

      for(unsigned icorner=0;icorner<6;icorner++)
	{
	  Vec3D rc = mirror->cornerInReflectorCoords(icorner,fFacetSize);
	  double r_ip = rc.Norm();
	  fReflectorIP = std::max(fReflectorIP, 2.0*r_ip);
	}

      id++;
    }
  
  // **************************************************************************
  // Clear the PIXELSs table and repopulate it with randomly generated pixels
  // **************************************************************************

  for(std::vector<VSOPixel*>::iterator i=fPixels.begin();
      i!=fPixels.end(); i++)delete *i;
  fPixels.clear();
  fPixelsByHexID.clear();  

  std::set<unsigned> pixels_missing;
  tokenize(param.PixelMissingList,pixels_missing);

  unsigned num_hex_pixel_rings = 
    unsigned(floor((fCameraDiameter/2.0)/fPixelSpacing/cos(30.0/180.0*M_PI)))+2;
  unsigned num_hex_pixel_sites = 
    3*num_hex_pixel_rings*(num_hex_pixel_rings+1)+1;

  id = 0;
  for(unsigned i=0; i<num_hex_pixel_sites; i++)
    {
      int hexid = i+1;

      if(pixels_missing.find(hexid) != pixels_missing.end())
	{
	  fPixelsByHexID.push_back(0);
	  continue; // skip pixel if on the missing list
	}

      Vec3D nominal_position;

      // Compute the pixel's nominal position
      nh_to_xy(&hexid, &nominal_position.x, &nominal_position.z); 
      if(fPixelParity)nominal_position.x=-nominal_position.x;
      nominal_position.x *= fPixelSpacing;
      nominal_position.z *= fPixelSpacing;
      nominal_position.y = 0;

      if(nominal_position.Norm() > fCameraDiameter/2.0)
	{
	  fPixelsByHexID.push_back(0);
	  continue; // skip pixel if position too far out
	}

      nominal_position.Rotate(fFPRotation);

      VSOPixel* pixel = new VSOPixel(this, id, hexid, false, nominal_position);
      
      fPixels.push_back(pixel);
      fPixelsByHexID.push_back(pixel);

      id++;
    }
}

void VSOTelescope::dumpShort(std::ostream& stream) const
{
  stream
    << "TELESCOPE "
    << VSDataConverter::toString(fID) << ' '
    << VSDataConverter::toString(fTelescopeHexID) << ' '
    << VSDataConverter::toString(fMirrors.size()) << ' '
    << VSDataConverter::toString(fPixels.size()) << ' '
    << VSDataConverter::toString(fPos.x) << ' '

    << VSDataConverter::toString(fPos.y) << ' '
    << VSDataConverter::toString(fPos.z) << ' '
    << VSDataConverter::toString(fDeltaY) << ' '
    << VSDataConverter::toString(fAlphaX) << ' '
    << VSDataConverter::toString(fAlphaY) << ' '

    << VSDataConverter::toString(fElevation) << ' '
    << VSDataConverter::toString(fAzimuth) << ' '
    << VSDataConverter::toString(fTranslation.x) << ' '
    << VSDataConverter::toString(fTranslation.y) << ' '
    << VSDataConverter::toString(fTranslation.z) << ' '

    << VSDataConverter::toString(fCurvatureRadius) << ' '
    << VSDataConverter::toString(fAperture) << ' '
    << VSDataConverter::toString(fFacetSpacing) << ' '
    << VSDataConverter::toString(fFacetSize) << ' '
    << VSDataConverter::toString(fReflectorRotation) << ' '

    << VSDataConverter::toString(fHexagonRingsN) << ' '
    << VSDataConverter::toString(fReflectorIP) << ' '
    << VSDataConverter::toString(fMirrorParity) << ' '
    << VSDataConverter::toString(fFPTranslation.x) << ' '
    << VSDataConverter::toString(fFPTranslation.y) << ' '

    << VSDataConverter::toString(fFPTranslation.z) << ' '
    << VSDataConverter::toString(fCameraDiameter) << ' '
    << VSDataConverter::toString(fFieldOfView) << ' '
    << VSDataConverter::toString(fCathodeDiameter) << ' '
    << VSDataConverter::toString(fPixelSpacing) << ' '

    << VSDataConverter::toString(fConcSurvProb) << ' '
    << VSDataConverter::toString(fFPRotation.x) << ' '
    << VSDataConverter::toString(fFPRotation.y) << ' '
    << VSDataConverter::toString(fFPRotation.z) << ' '
    << VSDataConverter::toString(fCameraIP) << ' '

    << VSDataConverter::toString(fPixelParity) << ' '
#if 0
    << VSDataConverter::toString(fHasSecondary) << ' '
    << VSDataConverter::toString(fRefractiveIndex) << ' '
    << VSDataConverter::toString(fRefractiveIndex1) << ' '
    << VSDataConverter::toString(fRefractiveIndex2) << ' '

    << VSDataConverter::toString(fCE1Parameter0) << ' '
    << VSDataConverter::toString(fCE1Parameter2) << ' '
    << VSDataConverter::toString(fCE1Parameter3) << ' '
    << VSDataConverter::toString(fCE1Parameter4) << ' '
    << VSDataConverter::toString(fCE1Parameter5) << ' '

    << VSDataConverter::toString(fCE1Parameter0) << ' '
    << VSDataConverter::toString(fCE1Parameter2) << ' '
    << VSDataConverter::toString(fCE1Parameter3) << ' '
    << VSDataConverter::toString(fCE1Parameter4) << ' '
    << VSDataConverter::toString(fCE1Parameter5) 
#endif
    << std::endl;

  for(std::vector<VSOMirror*> ::const_iterator i = fMirrors.begin();
      i!=fMirrors.end(); i++)
    (*i)->dumpShort(stream);

  for(std::vector<VSOPixel*> ::const_iterator i = fPixels.begin();
      i!=fPixels.end(); i++)
    (*i)->dumpShort(stream);
}

void VSOTelescope::dump(std::ostream& stream, unsigned l) const
{
  stream
    << FDAV("Telescope ID", fID, "", 30, l) << std::endl
    << FDAV("Telescope Hex ID", fTelescopeHexID, "", 30, l) << std::endl
    << FDAV("Num Mirrors", fMirrors.size(), "", 30, l) << std::endl
    << FDAV("Num Pixels", fPixels.size(), "", 30, l) << std::endl
    << FDAV("Position X", fPos.x, "cm", 30, l) << std::endl

    << FDAV("Position Y", fPos.y, "cm", 30, l) << std::endl
    << FDAV("Position Z", fPos.z, "cm", 30, l) << std::endl
    << FDAV("Delta Y", fDeltaY, "rad", 30, l) << std::endl
    << FDAV("Alpha X", fAlphaX, "rad", 30, l) << std::endl
    << FDAV("Alpha Y", fAlphaY, "rad", 30, l) << std::endl

    << FDAV("Elevation", fElevation, "rad", 30, l) << std::endl
    << FDAV("Azimuth", fAzimuth, "rad", 30, l) << std::endl
    << FDAV("Translation X", fTranslation.x, "cm", 30, l) << std::endl
    << FDAV("Translation Y", fTranslation.y, "cm", 30, l) << std::endl
    << FDAV("Translation Z", fTranslation.z, "cm", 30, l) << std::endl

    << FDAV("Curvature Radius", fCurvatureRadius, "cm", 30, l) << std::endl
    << FDAV("Aperture", fAperture, "cm", 30, l) << std::endl
    << FDAV("Facet Spacing", fFacetSpacing, "cm", 30, l) << std::endl
    << FDAV("Facet Size", fFacetSize, "cm", 30, l) << std::endl
    << FDAV("Reflector Rotation", fReflectorRotation, "rad", 30, l) << std::endl 

    << FDAV("Num Mirror Hexagon Rings", fHexagonRingsN, "", 30, l) << std::endl
    << FDAV("Reflector IP", fReflectorIP, "cm", 30, l) << std::endl
    << FDAV("Mirror Parity", fMirrorParity, "", 30, l) << std::endl
    << FDAV("FP Translation X", fFPTranslation.x, "cm", 30, l) << std::endl
    << FDAV("FP Translation Y", fFPTranslation.y, "cm", 30, l) << std::endl

    << FDAV("FP Translation Z", fFPTranslation.z, "cm", 30, l) << std::endl
    << FDAV("CameraDiameter", fCameraDiameter, "cm", 30, l) << std::endl
    << FDAV("Field Of View", fFieldOfView, "deg", 30, l) << std::endl
    << FDAV("Pixel Diameter", fCathodeDiameter, "cm", 30, l) << std::endl
    << FDAV("Pixel Spacing", fPixelSpacing, "cm", 30, l) << std::endl

    << FDAV("Conc Surv Prob", fConcSurvProb, "", 30, l) << std::endl
    << FDAV("FP Rotation X", fFPRotation.x, "rad", 30, l) << std::endl
    << FDAV("FP Rotation Y", fFPRotation.y, "rad", 30, l) << std::endl
    << FDAV("FP Rotation Z", fFPRotation.z, "rad", 30, l) << std::endl
    << FDAV("Camera IP", fCameraIP, "cm", 30, l) << std::endl

    << FDAV("Pixel Parity", fPixelParity, "", 30, l) << std::endl;
#if 0
    << FDAV("Has Secondary", fHasSecondary, "", 30, l) << std::endl
    << FDAV("Refractive Index", fRefractiveIndex, "", 30, l) << std::endl
    << FDAV("Refractive Index 1", fRefractiveIndex1, "", 30, l) << std::endl
    << FDAV("Refractive Index 2", fRefractiveIndex2, "", 30, l) << std::endl

    << FDAV("CE1 Parameter 0", fCE1Parameter0, "", 30, l) << std::endl
    << FDAV("CE1 Parameter 2", fCE1Parameter2, "", 30, l) << std::endl
    << FDAV("CE1 Parameter 3", fCE1Parameter3, "", 30, l) << std::endl
    << FDAV("CE1 Parameter 4", fCE1Parameter4, "", 30, l) << std::endl
    << FDAV("CE1 Parameter 5", fCE1Parameter5, "", 30, l) << std::endl

    << FDAV("CE2 Parameter 0", fCE1Parameter0, "", 30, l) << std::endl
    << FDAV("CE2 Parameter 2", fCE1Parameter2, "", 30, l) << std::endl
    << FDAV("CE2 Parameter 3", fCE1Parameter3, "", 30, l) << std::endl
    << FDAV("CE2 Parameter 4", fCE1Parameter4, "", 30, l) << std::endl
    << FDAV("CE2 Parameter 5", fCE1Parameter5, "", 30, l) << std::endl;
#endif

  for(std::vector<VSOMirror*> ::const_iterator i = fMirrors.begin();
      i!=fMirrors.end(); i++)
    {
      stream << std::endl;
      (*i)->dump(stream,l+1);
    }

  for(std::vector<VSOPixel*> ::const_iterator i = fPixels.begin();
      i!=fPixels.end(); i++)
    {
      stream << std::endl;
      (*i)->dump(stream,l+1);
    }
}

VSOTelescope* VSOTelescope::createFromShortDump(std::istream& stream)
{
  std::string line;
  std::getline(stream,line);
  if(line.empty())return 0;

  std::istringstream linestream(line);

  VSOTelescope* telescope = new VSOTelescope;

  std::string keyword;
  linestream >> keyword;
  assert(keyword==std::string("TELESCOPE"));

  unsigned mirrors_size;
  unsigned pixels_size;

  linestream
    >> telescope->fID
    >> telescope->fTelescopeHexID
    >> mirrors_size
    >> pixels_size
    >> telescope->fPos.x
    
    >> telescope->fPos.y
    >> telescope->fPos.z
    >> telescope->fDeltaY
    >> telescope->fAlphaX
    >> telescope->fAlphaY

    >> telescope->fElevation
    >> telescope->fAzimuth
    >> telescope->fTranslation.x
    >> telescope->fTranslation.y
    >> telescope->fTranslation.z

    >> telescope->fCurvatureRadius
    >> telescope->fAperture
    >> telescope->fFacetSpacing
    >> telescope->fFacetSize
    >> telescope->fReflectorRotation

    >> telescope->fHexagonRingsN
    >> telescope->fReflectorIP
    >> telescope->fMirrorParity
    >> telescope->fFPTranslation.x
    >> telescope->fFPTranslation.y

    >> telescope->fFPTranslation.z
    >> telescope->fCameraDiameter
    >> telescope->fFieldOfView
    >> telescope->fCathodeDiameter
    >> telescope->fPixelSpacing

    >> telescope->fConcSurvProb
    >> telescope->fFPRotation.x
    >> telescope->fFPRotation.y
    >> telescope->fFPRotation.z
    >> telescope->fCameraIP

    >> telescope->fPixelParity
#if 0
    >> telescope->fHasSecondary
    >> telescope->fRefractiveIndex
    >> telescope->fRefractiveIndex1
    >> telescope->fRefractiveIndex2

    >> telescope->fCE1Parameter0
    >> telescope->fCE1Parameter2
    >> telescope->fCE1Parameter3
    >> telescope->fCE1Parameter4
    >> telescope->fCE1Parameter5

    >> telescope->fCE1Parameter0
    >> telescope->fCE1Parameter2
    >> telescope->fCE1Parameter3
    >> telescope->fCE1Parameter4
    >> telescope->fCE1Parameter5
#endif
    ;  
  
  if(!linestream)
    {
      delete telescope;
      return 0;
    }

  for(unsigned i=0; i<mirrors_size; i++)
    {
      VSOMirror* mirror = VSOMirror::createFromShortDump(stream, telescope);
      if(mirror==0)
	{
	  delete telescope;
	  return 0;
	}
      
      if(mirror->id() >= telescope->fMirrors.size())
	telescope->fMirrors.resize(mirror->id()+1);
      telescope->fMirrors[mirror->id()]=mirror;
      
      if(mirror->hexID() > telescope->fMirrorsByHexID.size())
	telescope->fMirrorsByHexID.resize(mirror->hexID());
      telescope->fMirrorsByHexID[mirror->hexID()-1]=mirror;
    }

  for(unsigned i=0; i<pixels_size; i++)
    {
      VSOPixel* pixel = VSOPixel::createFromShortDump(stream, telescope);
      if(pixel==0)
	{
	  delete telescope;
	  return 0;
	}
      
      if(pixel->id() >= telescope->fPixels.size())
	telescope->fPixels.resize(pixel->id()+1);
      telescope->fPixels[pixel->id()]=pixel;
      
      if(pixel->hexID() > telescope->fPixelsByHexID.size())
	telescope->fPixelsByHexID.resize(pixel->hexID());
      telescope->fPixelsByHexID[pixel->hexID()-1]=pixel;
    }

  // Recalculate rotation vector
  telescope->calculateRotationVector();
  
  return telescope;
}

#endif
