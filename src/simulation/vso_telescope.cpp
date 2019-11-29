/*

   calin/simulation/vso_telescope.cpp -- Stephen Fegan -- 2015-11-11

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

#include <sstream>
#include <algorithm>

#include <util/log.hpp>
#include <math/special.hpp>
#include <math/hex_array.hpp>
#include <math/regular_grid.hpp>
#include <math/vector3d_util.hpp>
#include <simulation/vso_telescope.hpp>
#include <simulation/vs_optics.pb.h>
#include <iact_data/instrument_layout.hpp>

using namespace calin::util::log;
using namespace calin::simulation::vs_optics;
using calin::math::special::SQR;

VSOTelescope::VSOTelescope():
    fID(), fPos(),
    fDeltaY(), fAlphaX(), fAlphaY(),
    fElevation(), fAzimuth(),
    fTranslation(), fCurvatureRadius(), fAperture(),
    fFacetSpacing(), fFacetSize(),
    fReflectorRotation(), fCosReflectorRotation(1.0), fSinReflectorRotation(),
    fFacetGridShiftX(), fFacetGridShiftZ(),
    fHexagonRingsN(), fReflectorIP(), fReflectorIPCenter(),
    fMirrorParity(), fFPTranslation(), fCameraDiameter(),
    fFieldOfView(), fCathodeDiameter(), fPixelSpacing(),
    fPixelRotation(), fCosPixelRotation(1.0), fSinPixelRotation(0.0),
    fPixelGridShiftX(), fPixelGridShiftZ(), fConcSurvProb(),
    fFPRotation(), fCameraIP(), fPixelParity(),
    fWindowFront(), fWindowOuterRadius(), fWindowThickness(), fWindowRefractiveIndex(),
    fPreObscurations(), fPostObscurations(), fCameraObscurations(),
    fMirrors(), fMirrorsByHexID(), fPixels(), fPixelsByHexID(), rot_reflector_to_global_()
{
  calculateFPRotationMatrix();
  calculateRotationVector();
}

VSOTelescope::
VSOTelescope(unsigned TID, const Eigen::Vector3d&P,
	     double DY, double AX, double AY, double EL, double AZ,
	     const Eigen::Vector3d& T, double CR, double A, double FSP, double FS,
	     double RR, double FGSX, double FGSZ,
       unsigned HRN, double RIP, const Eigen::Vector3d& RIPC, bool MP,
	     const Eigen::Vector3d& FPT, double CD, double FOV,
       double D, double PS, double PR, double PGSX, double PGSZ,
	     double CSP, const Eigen::Vector3d& FPR, double CIP, bool PP,
       double WIN_FRONT, double WIN_RAD, double WIN_THICK, double WIN_N,
	     const std::vector<VSOObscuration*>& OBSVEC_PRE,
       const std::vector<VSOObscuration*>& OBSVEC_POST,
       const std::vector<VSOObscuration*>& OBSVEC_CAM
	     ):
    fID(TID), /*fTelescopeHexID(THID),*/ fPos(P),
    fDeltaY(DY), fAlphaX(AX), fAlphaY(AY), fElevation(EL), fAzimuth(AZ),
    fTranslation(T), fCurvatureRadius(CR), fAperture(A), fFacetSpacing(FSP),
    fFacetSize(FS),
    fReflectorRotation(RR), fCosReflectorRotation(std::cos(RR)), fSinReflectorRotation(std::sin(RR)),
    fFacetGridShiftX(FGSX), fFacetGridShiftZ(FGSZ),
    fHexagonRingsN(HRN), fReflectorIP(RIP), fReflectorIPCenter(RIPC), fMirrorParity(MP),
    fFPTranslation(FPT),  fCameraDiameter(CD), fFieldOfView(FOV),
    fCathodeDiameter(D), fPixelSpacing(PS), fPixelRotation(PR),
    fCosPixelRotation(std::cos(PR)), fSinPixelRotation(std::sin(PR)),
    fPixelGridShiftX(PGSX), fPixelGridShiftZ(PGSZ),
    fConcSurvProb(CSP),
    fFPRotation(FPR), fCameraIP(CIP), fPixelParity(PP),
    fWindowFront(WIN_FRONT), fWindowOuterRadius(WIN_RAD),
    fWindowThickness(WIN_THICK), fWindowRefractiveIndex(WIN_N),
    fPreObscurations(OBSVEC_PRE), fPostObscurations(OBSVEC_POST), fCameraObscurations(OBSVEC_CAM),
    fMirrors(), fMirrorsByHexID(), fPixels(), fPixelsByHexID(), rot_reflector_to_global_()
{
  calculateFPRotationMatrix();
  calculateRotationVector();
}

VSOTelescope::VSOTelescope(const VSOTelescope& o):
    fID(o.fID),
    fPos(o.fPos), fDeltaY(o.fDeltaY), fAlphaX(o.fAlphaX), fAlphaY(o.fAlphaY),
    fElevation(o.fElevation), fAzimuth(o.fAzimuth), fTranslation(o.fTranslation),
    fCurvatureRadius(o.fCurvatureRadius), fAperture(o.fAperture),
    fFacetSpacing(o.fFacetSpacing), fFacetSize(o.fFacetSize),
    fReflectorRotation(o.fReflectorRotation),
    fCosReflectorRotation(o.fCosReflectorRotation),
    fSinReflectorRotation(o.fSinReflectorRotation),
    fFacetGridShiftX(o.fFacetGridShiftX),
    fFacetGridShiftZ(o.fFacetGridShiftZ),
    fHexagonRingsN(o.fHexagonRingsN),
    fReflectorIP(o.fReflectorIP), fReflectorIPCenter(o.fReflectorIPCenter),
    fMirrorParity(o.fMirrorParity),
    fFPTranslation(o.fFPTranslation), fCameraDiameter(o.fCameraDiameter),
    fFieldOfView(o.fFieldOfView), fCathodeDiameter(o.fCathodeDiameter),
    fPixelSpacing(o.fPixelSpacing), fPixelRotation(o.fPixelRotation),
    fCosPixelRotation(o.fCosPixelRotation),
    fSinPixelRotation(o.fSinPixelRotation),
    fPixelGridShiftX(o.fPixelGridShiftX),
    fPixelGridShiftZ(o.fPixelGridShiftZ),
    fConcSurvProb(o.fConcSurvProb),
    fFPRotation(o.fFPRotation), fCameraIP(o.fCameraIP),
    fPixelParity(o.fPixelParity),
    fWindowFront(o.fWindowFront), fWindowOuterRadius(o.fWindowOuterRadius),
    fWindowThickness(o.fWindowThickness), fWindowRefractiveIndex(o.fWindowRefractiveIndex),
    fPreObscurations(), fPostObscurations(), fCameraObscurations(),
    fMirrors(), fMirrorsByHexID(), fPixels(), fPixelsByHexID(),
    rot_reflector_to_global_()
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

  for(const auto *i : o.fPreObscurations)fPreObscurations.push_back(i->clone());
  for(const auto *i : o.fPostObscurations)fPostObscurations.push_back(i->clone());
  for(const auto *i : o.fCameraObscurations)fCameraObscurations.push_back(i->clone());

  calculateFPRotationMatrix();
  calculateRotationVector();
}

VSOTelescope::~VSOTelescope()
{
  for(std::vector<VSOMirror*>::iterator i=fMirrors.begin();
      i!=fMirrors.end(); i++)delete *i;
  for(std::vector<VSOPixel*>::iterator i=fPixels.begin();
      i!=fPixels.end(); i++)delete *i;
  for(auto* i : fPreObscurations)delete i;
  for(auto* i : fPostObscurations)delete i;
  for(auto* i : fCameraObscurations)delete i;
}

// Stroustrup third edition sec 11.3.4 recommends default copy assignment
const VSOTelescope& VSOTelescope::operator =(const VSOTelescope& o)
{
  fID                = o.fID;
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
  fCosReflectorRotation = o.fCosReflectorRotation;
  fSinReflectorRotation = o.fSinReflectorRotation;
  fFacetGridShiftX   = o.fFacetGridShiftX;
  fFacetGridShiftZ   = o.fFacetGridShiftZ;
  fHexagonRingsN     = o.fHexagonRingsN;
  fReflectorIP       = o.fReflectorIP;
  fReflectorIPCenter = o.fReflectorIPCenter;
  fMirrorParity      = o.fMirrorParity;
  fFPTranslation     = o.fFPTranslation;
  fCameraDiameter    = o.fCameraDiameter;
  fFieldOfView       = o.fFieldOfView;
  fCathodeDiameter   = o.fCathodeDiameter;
  fPixelSpacing      = o.fPixelSpacing;
  fPixelRotation     = o.fPixelRotation;
  fCosPixelRotation  = o.fCosPixelRotation;
  fSinPixelRotation  = o.fSinPixelRotation;
  fPixelGridShiftX   = o.fPixelGridShiftX;
  fPixelGridShiftZ   = o.fPixelGridShiftZ;
  fConcSurvProb      = o.fConcSurvProb;
  fFPRotation        = o.fFPRotation;
  fCameraIP          = o.fCameraIP;
  fPixelParity       = o.fPixelParity;
  fWindowFront      = o.fWindowFront;
  fWindowOuterRadius = o.fWindowOuterRadius;
  fWindowThickness   = o.fWindowThickness;
  fWindowRefractiveIndex = o.fWindowRefractiveIndex;

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

  for(auto* i : fPreObscurations)delete i;
  fPreObscurations.clear();
  for(auto* i : o.fPreObscurations)fPreObscurations.push_back(i->clone());

  for(auto* i : fPostObscurations)delete i;
  fPostObscurations.clear();
  for(auto* i : o.fPostObscurations)fPostObscurations.push_back(i->clone());

  for(auto* i : fCameraObscurations)delete i;
  fCameraObscurations.clear();
  for(auto* i : o.fCameraObscurations)fCameraObscurations.push_back(i->clone());

  calculateFPRotationMatrix();
  calculateRotationVector();

  return *this;
}

// ****************************************************************************
// Accessor
// ****************************************************************************

Eigen::Vector3d VSOTelescope::opticalAxis() const
{
  return Eigen::Vector3d(cos(fElevation)*sin(fAzimuth),
	       cos(fElevation)*cos(fAzimuth),
	       sin(fElevation));
}

// ****************************************************************************
// Repoint the telescope along a vector
// ****************************************************************************

bool VSOTelescope::pointTelescope(const Eigen::Vector3d& v)
{
  if(v.squaredNorm()==0)return false;
  fElevation = atan2(v.z(),sqrt(v.x()*v.x() + v.y()*v.y()));
  fAzimuth = fmod(atan2(v.x(),v.y())+2.0*M_PI, 2.0*M_PI);
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

void VSOTelescope::calculateFPRotationMatrix()
{
  double fp_rot_norm = fFPRotation.norm();
  if(std::fabs(fp_rot_norm) > 0)
  {
    fHasFPRotation = false;
    rot_camera_to_reflector_ = Eigen::Matrix3d::Identity();
    rot_reflector_to_camera_ = Eigen::Matrix3d::Identity();
  }
  else
  {
    fHasFPRotation = false;
    rot_camera_to_reflector_ = Eigen::AngleAxisd(fp_rot_norm,
      fFPRotation.normalized()).toRotationMatrix();
    rot_reflector_to_camera_ = rot_camera_to_reflector_.transpose();
  }
}

void VSOTelescope::calculateRotationVector()
{
  // Rotation vector maps from Reflector to Global
  rot_reflector_to_global_ =
    Eigen::AngleAxisd(fAlphaY,     Eigen::Vector3d::UnitX()) *
    Eigen::AngleAxisd(fAlphaX,     Eigen::Vector3d::UnitY()) *
    Eigen::AngleAxisd(-fAzimuth,   Eigen::Vector3d::UnitZ()) *
    Eigen::AngleAxisd(fDeltaY,     Eigen::Vector3d::UnitY()) *
    Eigen::AngleAxisd(fElevation,  Eigen::Vector3d::UnitX());
  rot_global_to_reflector_ = rot_reflector_to_global_.transpose();
  off_global_to_reflector_ = fPos - rot_reflector_to_global_ * fTranslation;

  rot_camera_to_global_ = rot_reflector_to_global_ * rot_camera_to_reflector_;
  off_global_to_camera_ =
    fPos - rot_reflector_to_global_ * (fTranslation - fFPTranslation);

#if 0
  rot_reflector_to_global_ =
      Eigen::Vector3d(1,0,0)*fElevation &
      Eigen::Vector3d(0,1,0)*fDeltaY &
      Eigen::Vector3d(0,0,-1)*fAzimuth &
      Eigen::Vector3d(0,1,0)*fAlphaX &
      Eigen::Vector3d(1,0,0)*fAlphaY;
#endif
}

// ****************************************************************************
// Fill the mirrors and pixels tables up with randomly generated objects
// ****************************************************************************

void VSOTelescope::
populateMirrorsAndPixelsRandom(
    const ix::simulation::vs_optics::IsotropicDCArrayParameters& param,
    math::rng::RNG& rng)
{
  // **************************************************************************
  // Clear the MIRRORs table and repopulate it with randomly generated mirrors
  // **************************************************************************
  const double reflector_r2 = fAperture*fAperture/4.0;
  const double reflector_c2 = fCurvatureRadius*fCurvatureRadius;
  Eigen::Vector3d reflector_center(0, fCurvatureRadius, 0);

  const double cos_reflector_rot = cosReflectorRotation();
  const double sin_reflector_rot = sinReflectorRotation();

  std::set<unsigned> mirrors_missing;
  for(auto iid : param.reflector().facet_missing_list())
    mirrors_missing.insert(iid);

  for(std::vector<VSOMirror*>::iterator i=fMirrors.begin();
      i!=fMirrors.end(); i++)delete *i;
  fMirrors.clear();
  fMirrorsByHexID.clear();

  int num_hex_rings = fHexagonRingsN;
  if(num_hex_rings <= 0)
  {
    num_hex_rings = -1;
    bool has_mirror = true;
    while(has_mirror)
    {
      ++num_hex_rings;
      has_mirror = false;
      for(int irun=0;irun<(num_hex_rings+1);irun++)
      {
        auto hexid = math::hex_array::
          positive_ringid_segid_runid_to_hexid(num_hex_rings+1,0,irun);
        double x;
        double z;
        math::hex_array::hexid_to_xy_trans(hexid, x, z, fMirrorParity,
          cos_reflector_rot, sin_reflector_rot, fFacetSpacing,
          facetGridShiftX(), facetGridShiftZ());
        if(SQR(x)+SQR(z) <= reflector_r2)
        {
          has_mirror = true;
          break;
        }
      }
    }
  }

  int num_hex_mirror_sites =
      math::hex_array::ringid_to_nsites_contained(num_hex_rings);

  unsigned id = 0;
  for(int hexid=0; hexid<num_hex_mirror_sites; hexid++)
  {
    if(mirrors_missing.find(hexid) != mirrors_missing.end())
    {
      fMirrorsByHexID.push_back(0);
      continue; // skip mirror if on the mirring list
    }

    Eigen::Vector3d nominal_position(0,0,0);

    // Compute the mirror's nominal position
    math::hex_array::hexid_to_xy_trans(hexid,
        nominal_position.x(), nominal_position.z(), fMirrorParity,
        cos_reflector_rot, sin_reflector_rot, fFacetSpacing,
        facetGridShiftX(), facetGridShiftZ());

    double X2 = SQR(nominal_position.x());
    double Z2 = SQR(nominal_position.z());

    if( (reflector_r2 - X2 - Z2) <= 0 )
    {
      fMirrorsByHexID.push_back(0);
      continue; // skip mirror if projected hex position too far out
    }

    nominal_position.y() = fCurvatureRadius-sqrt(reflector_c2 - X2 - Z2);

    // Add Gaussian normal/tangenetial position error
    Eigen::Vector3d reflector_normal(reflector_center-nominal_position);
    reflector_normal.normalize();

    Eigen::Vector3d position(nominal_position);
    position += reflector_normal*(rng.normal()*
                             param.reflector().facet_pos_normal_dispersion());
    position -= reflector_center;
    calin::math::vector3d_util::scatter_direction(position,
      param.reflector().facet_pos_tangent_dispersion()/fCurvatureRadius, rng);
    position += reflector_center;

    // Get the (perturbed) alignment angle of mirror
    Eigen::Vector3d alignment(0,0,0);

    if(param.reflector().alignment_case() ==
       ix::simulation::vs_optics::HexDCReflectorParameters::kNormAlign)
    {
      Eigen::Vector3d alignment_pt;
      calin::math::vector3d_util::set_from_proto(alignment_pt,
        param.reflector().norm_align().alignment_pt());
      if(alignment_pt.y() == 0)
      {
        double F = param.reflector().alignment_image_plane();
        if(F==0)F = fFPTranslation.y();
        double d = (Eigen::Vector3d(0,F,0)-nominal_position).norm();
        alignment_pt.y() = F+d;
      }

      // Standard DC alignment to a fixed point in space (with scatter)
      alignment = (alignment_pt-position).normalized();

      double align_disp =
        param.reflector().facet_alignment_dispersion()/alignment_pt.norm();
      calin::math::vector3d_util::scatter_direction(alignment, align_disp, rng);
    }
    else if(param.reflector().alignment_case() ==
            ix::simulation::vs_optics::HexDCReflectorParameters::kPsfAlign)
    {
      // Alignment to a point in the focal plane

      double theta =
          param.reflector().psf_align().image_pt_theta()*M_PI/180.0;
      double stheta = sin(theta);
      double ctheta = cos(theta);
      double ttheta = stheta/ctheta;

      double phi =
          param.reflector().psf_align().image_pt_phi()*M_PI/180.0;
      double sphi = sin(phi);
      double cphi = cos(phi);

      double y_fp = param.reflector().alignment_image_plane(); // Image plane
      if(y_fp==0)y_fp = fFPTranslation.y();
      Eigen::Vector3d r_fp(y_fp*ttheta*sphi, y_fp, -y_fp*ttheta*cphi);

      Eigen::Vector3d e_in(-stheta*sphi, ctheta, stheta*cphi);
      if((param.reflector().psf_align().object_plane()) and
         (!std::isinf(param.reflector().psf_align().object_plane())))
      {
        double y_em = param.reflector().psf_align().object_plane();
        Eigen::Vector3d emission_pt(-y_em*ttheta*sphi, y_em, y_em*ttheta*cphi);
        e_in = (emission_pt-position).normalized();
      }

      Eigen::Vector3d e_out = (r_fp-position).normalized();

      alignment = e_in;

      Eigen::Vector3d e_rot = e_in.cross(e_out);
      double strot2 = e_rot.norm();
      if(strot2 != 0)
      {
        double ctrot2 = e_in.dot(e_out);
        double trot2 = std::atan2(strot2,ctrot2);
        alignment = Eigen::AngleAxisd(0.5*trot2, e_rot/strot2) * alignment;
      }
    }
    else
    {
      throw std::runtime_error("Unknown alignment scheme");
    }

    double focal_length = param.reflector().facet_focal_length();
    if(focal_length == 0)
      focal_length = (nominal_position-fFPTranslation).norm();
    focal_length += rng.normal()*param.reflector().facet_focal_length_dispersion();

    double spot_size;
    if(param.reflector().facet_spot_size_dispersion() > 0)
      spot_size = rng.gamma_by_mean_and_sigma(
          param.reflector().facet_spot_size(),
          param.reflector().facet_spot_size_dispersion());
    else spot_size = param.reflector().facet_spot_size();

    // If param.MirrorSpotSizePhotonFraction not set -- assume FWHM
    // ... spot_size comes in as diameter at some fraction, goes out as
    //     diameter that corresponds to dispersion
    if((param.reflector().facet_spot_size_probability()>0.0)&&
       (param.reflector().facet_spot_size_probability()<1.0))
      spot_size =
        spot_size/std::sqrt(2.0*std::log(1.0/(1.0-
                          param.reflector().facet_spot_size_probability())));
    else
      spot_size =
        spot_size/std::sqrt(2.0*log(2.0));

    VSOMirror* mirror =
	    new VSOMirror(this, id, hexid, false, position, alignment,
                      focal_length, spot_size,
                      param.reflector().weathering_factor());

    fMirrors.push_back(mirror);
    fMirrorsByHexID.push_back(mirror);

    for(unsigned icorner=0;icorner<6;icorner++)
    {
      Eigen::Vector3d rc = mirror->cornerInReflectorCoords(icorner,fFacetSize);
      double d_ip = 2.0*std::sqrt(SQR(rc.x()) + SQR(rc.z()));
      if(d_ip > fReflectorIP)
      {
        fReflectorIP = d_ip;
        fReflectorIPCenter = Eigen::Vector3d(0,rc.y(),0);
      }
    }

    id++;
  }

  // **************************************************************************
  // Clear the PIXELSs table and repopulate it
  // **************************************************************************

  for(std::vector<VSOPixel*>::iterator i=fPixels.begin();
      i!=fPixels.end(); i++)delete *i;
  fPixels.clear();
  fPixelsByHexID.clear();

  std::set<unsigned> modules_missing;
  for(auto iid : param.pixel().module_missing_list())
    modules_missing.insert(iid);

  unsigned module_size = param.pixel().hex_module_size();
  bool use_a_config = true;
  if(param.pixel().hex_module_layout_use_b_configuration())use_a_config = false;
  std::set<unsigned> pixel_hexids;
  if(modules_missing.find(0) == modules_missing.end())
    for(auto id : math::hex_array::cluster_hexid_to_member_hexid(0, module_size, use_a_config))
      pixel_hexids.insert(id);
  unsigned module_ring_id = 1;
  bool module_ring_in_camera = true;
  if(param.pixel().module_num_hex_rings()==0 and fCameraDiameter<=0)
    module_ring_in_camera = false;
  while(module_ring_in_camera
        and (param.pixel().module_num_hex_rings()==0
             or module_ring_id<=param.pixel().module_num_hex_rings()))
  {
    module_ring_in_camera = false;
    for(unsigned iseg=0;iseg<6;iseg++)
      for(unsigned irun=0;irun<module_ring_id;irun++)
      {
        unsigned module_id =
            math::hex_array::
            ringid_segid_runid_to_hexid(module_ring_id, iseg, irun);
        if(modules_missing.find(module_id) != modules_missing.end())
          continue;
        int u, v;
        math::hex_array::cluster_hexid_to_center_uv(module_id, module_size, u, v, use_a_config);
        double x, z;
        math::hex_array::uv_to_xy_trans(u, v, x, z,
          fCosPixelRotation, fSinPixelRotation, fPixelSpacing,
          fPixelGridShiftX, fPixelGridShiftZ);
        if(fCameraDiameter>0 and
            std::sqrt(x*x+z*z)>fCameraDiameter/2.0)
	        continue; // skip module if position too far out
        module_ring_in_camera = true;
        for(auto id : math::hex_array::
            cluster_hexid_to_member_hexid(module_id, module_size, use_a_config))
          pixel_hexids.insert(id);
      }
    module_ring_id++;
  }

  double cos_fp_rot = fCosPixelRotation;
  double sin_fp_rot = fSinPixelRotation;
  unsigned pixelid = 0;
  for(auto hexid : pixel_hexids)
  {
    while(fPixelsByHexID.size()<hexid)fPixelsByHexID.push_back(0);
    Eigen::Vector3d nominal_position;

    // Compute the pixel's nominal position
    math::hex_array::hexid_to_xy_trans(hexid,
      nominal_position.x(), nominal_position.z(), fPixelParity,
      cos_fp_rot, sin_fp_rot, fPixelSpacing, fPixelGridShiftX, fPixelGridShiftZ);
    nominal_position.y() = 0;
    nominal_position = rot_camera_to_reflector_ * nominal_position;

    VSOPixel* pixel =
        new VSOPixel(this, pixelid, hexid, false, nominal_position);

    fPixels.push_back(pixel);
    fPixelsByHexID.push_back(pixel);

    pixelid++;
  }
}

calin::ix::simulation::vs_optics::VSOTelescopeData* VSOTelescope::
dump_as_proto(calin::ix::simulation::vs_optics::VSOTelescopeData* d) const
{
  if(d == nullptr)d = new calin::ix::simulation::vs_optics::VSOTelescopeData;

  d->Clear();
  d->set_id(fID);
  calin::math::vector3d_util::dump_as_proto(fPos, d->mutable_pos());
  d->set_delta_y(fDeltaY/M_PI*180.0);
  d->set_alpha_x(fAlphaX/M_PI*180.0);
  d->set_alpha_y(fAlphaY/M_PI*180.0);
  d->mutable_alt_az()->set_altitude(fElevation/M_PI*180.0);
  d->mutable_alt_az()->set_azimuth(fAzimuth/M_PI*180.0);
  calin::math::vector3d_util::dump_as_proto(fTranslation, d->mutable_translation());
  d->set_curvature_radius(fCurvatureRadius);
  d->set_aperture(fAperture);
  d->set_facet_spacing(fFacetSpacing);
  d->set_facet_size(fFacetSize);
  d->set_optic_axis_rotation(fReflectorRotation/M_PI*180.0);
  d->set_facet_grid_shift_x(fFacetGridShiftX);
  d->set_facet_grid_shift_z(fFacetGridShiftZ);
  d->set_hexagon_rings_n(fHexagonRingsN);
  d->set_reflector_ip(fReflectorIP);
  calin::math::vector3d_util::dump_as_proto(fReflectorIPCenter, d->mutable_reflector_ip_center());
  d->set_facet_labeling_parity(fMirrorParity);
  calin::math::vector3d_util::dump_as_proto(fFPTranslation, d->mutable_fp_translation());
  d->set_camera_diameter(fCameraDiameter);
  d->set_field_of_view(fFieldOfView);
  d->set_cathode_diameter(fCathodeDiameter);
  d->set_pixel_spacing(fPixelSpacing);
  d->set_pixel_rotation(fPixelRotation/M_PI*180.0);
  d->set_pixel_grid_shift_x(fPixelGridShiftX);
  d->set_pixel_grid_shift_z(fPixelGridShiftZ);
  d->set_conc_survival_prob(fConcSurvProb);
  calin::math::vector3d_util::dump_as_scaled_proto(fFPRotation, 180.0/M_PI, d->mutable_fp_rotation());
  d->set_camera_ip(fCameraIP);
  d->set_pixel_labeling_parity(fPixelParity);
  d->set_window_front(fWindowFront);
  d->set_window_outer_radius(fWindowOuterRadius);
  d->set_window_thickness(fWindowThickness);
  d->set_window_refractive_index(fWindowRefractiveIndex);

  for(auto iobs : fPreObscurations)
    iobs->dump_as_proto(d->add_pre_reflection_obscuration());
  for(auto iobs : fPostObscurations)
    iobs->dump_as_proto(d->add_post_reflection_obscuration());
  for(auto iobs : fCameraObscurations)
    iobs->dump_as_proto(d->add_camera_obscuration());
  for(auto imir : fMirrors)
    if(imir != nullptr)imir->dump_as_proto(d->add_mirror());
  for(auto ipix : fPixels)
    if(ipix != nullptr)ipix->dump_as_proto(d->add_pixel());

  return d;
}

VSOTelescope* VSOTelescope::
create_from_proto(const ix::simulation::vs_optics::VSOTelescopeData& d)
{
  VSOTelescope* scope = new VSOTelescope(
    d.id(), // TID
    calin::math::vector3d_util::from_proto(d.pos()), // P
    d.delta_y()*M_PI/180.0, // DY
    d.alpha_x()*M_PI/180.0, // AX
    d.alpha_y()*M_PI/180.0, // AY
    d.alt_az().altitude()*M_PI/180.0, // EL
    d.alt_az().azimuth()*M_PI/180.0,  // AZ
    calin::math::vector3d_util::from_proto(d.translation()), // T
    d.curvature_radius(),  // CR
    d.aperture(),  // A
    d.facet_spacing(), // FSP
    d.facet_size(), // FS
    d.optic_axis_rotation()*M_PI/180.0, // RR
    d.facet_grid_shift_x(), // FGSX
    d.facet_grid_shift_z(), // FGSZ
    d.hexagon_rings_n(), // HRN
    d.reflector_ip(), // RIP
    calin::math::vector3d_util::from_proto(d.reflector_ip_center()), // RIPC
    d.facet_labeling_parity(), // MP
    calin::math::vector3d_util::from_proto(d.fp_translation()), // FPT
    d.camera_diameter(), // CD
    d.field_of_view(), // FOV
    d.cathode_diameter(), // D
    d.pixel_spacing(), // PS
    d.pixel_grid_shift_x(), // PGSX
    d.pixel_grid_shift_z(), // PGSZ
    d.pixel_rotation()*M_PI/180.0, // PR
    d.conc_survival_prob(), // CSP
    calin::math::vector3d_util::from_scaled_proto(d.fp_rotation(), M_PI/180.0), // FPR
    d.camera_ip(), // CIP
    d.pixel_labeling_parity(), // PP
    d.window_front(),
    d.window_outer_radius(),
    d.window_thickness(),
    d.window_refractive_index()
    );

  for(int i=0; i<d.pre_reflection_obscuration_size(); i++)
    scope->add_pre_reflection_obscuration(VSOObscuration::create_from_proto(d.pre_reflection_obscuration(i)));
  for(int i=0; i<d.post_reflection_obscuration_size(); i++)
    scope->add_post_reflection_obscuration(VSOObscuration::create_from_proto(d.post_reflection_obscuration(i)));

  for(int i=0; i<d.mirror_size(); i++)
    scope->add_mirror(VSOMirror::create_from_proto(d.mirror(i),scope));
  for(int i=0; i<d.pixel_size(); i++)
    scope->add_pixel(VSOPixel::create_from_proto(d.pixel(i),scope));

  return scope;
}

calin::ix::iact_data::instrument_layout::TelescopeLayout*
VSOTelescope::convert_to_telescope_layout(
  calin::ix::iact_data::instrument_layout::TelescopeLayout* d) const
{
  if(d == nullptr)d = new calin::ix::iact_data::instrument_layout::TelescopeLayout;

  d->Clear();

  d->set_telescope_type(calin::ix::iact_data::instrument_layout::TelescopeLayout::NO_TELESCOPE);
  d->set_telescope_index(0);
  calin::math::vector3d_util::dump_as_proto(fPos, d->mutable_position());
  d->set_effective_focal_length(fFPTranslation.y());

  auto* c = d->mutable_camera();

  c->set_camera_type(calin::ix::iact_data::instrument_layout::CameraLayout::NO_CAMERA);
  c->set_camera_number(fID);

  c->set_pixel_grid_layout(calin::ix::iact_data::instrument_layout::CameraLayout::HEX_GRID);
  c->set_pixel_grid_spacing(fPixelSpacing);
  c->set_pixel_grid_rotation(fPixelRotation/M_PI*180.0);
  c->set_pixel_grid_cos_rotation(fCosPixelRotation);
  c->set_pixel_grid_sin_rotation(fSinPixelRotation);
  c->set_pixel_grid_offset_x(fPixelGridShiftX);
  c->set_pixel_grid_offset_y(fPixelGridShiftZ);
  c->set_pixel_grid_geometric_area(calin::math::hex_array::cell_area(fPixelSpacing));
  // What do we do with : fPixelParity

  math::regular_grid::HexGrid grid(fPixelSpacing, fPixelRotation,
    fPixelGridShiftX, fPixelGridShiftZ);
  for(const auto* ipix : fPixels)
  {
    auto* ch = c->add_channel();
    ch->set_channel_index(ipix->id());
    ch->set_pixel_index(ipix->id());
    ch->set_pixel_spiral_index(ipix->id());
    ch->set_pixel_grid_index(ipix->hexID());

    c->add_pixel_channel_index(ipix->id());
    c->add_pixel_spiral_channel_index(ipix->id());

    ch->set_channel_set_index(0);
    ch->set_module_index(-1);
    ch->set_module_channel_index(-1);
    double x;
    double y;
    grid.gridid_to_xy(ipix->hexID(), x, y);
    ch->set_x(x);
    ch->set_y(y);
    ch->set_diameter(fPixelSpacing);
    ch->set_geometric_area(grid.cell_area(ipix->hexID()));
    auto nbr = grid.gridid_to_neighbour_gridids(ipix->hexID());
    ch->set_is_boundary_pixel(false);
    for(auto inbr : nbr) {
      const VSOPixel* npix = this->pixelByHexID(inbr);
      if(npix)ch->add_neighbour_channel_indexes(npix->id());
      else ch->set_is_boundary_pixel(true);
    }
    Eigen::VectorXd xv;
    Eigen::VectorXd yv;
    grid.gridid_to_vertexes_xy(ipix->hexID(), xv, yv);
    for(unsigned i=0; i<xv.size(); i++)
      ch->add_outline_polygon_vertex_x(xv(i));
    for(unsigned i=0; i<yv.size(); i++)
      ch->add_outline_polygon_vertex_y(yv(i));
    ch->add_outline_polygon_vertex_index(xv.size());

    if(ch->is_boundary_pixel()) {
      c->add_boundary_pixel_channel_index(ipix->id());
    }
  }

  calin::iact_data::instrument_layout::compute_camera_and_module_outlines(c);

  return d;
}

calin::ix::iact_data::instrument_layout::TelescopeLayout*
calin::simulation::vs_optics::dc_parameters_to_telescope_layout(
  const ix::simulation::vs_optics::IsotropicDCArrayParameters& param,
  unsigned telescope_id, const Eigen::Vector3d& pos,
  calin::ix::iact_data::instrument_layout::TelescopeLayout* d)
{
  if(d == nullptr)d = new calin::ix::iact_data::instrument_layout::TelescopeLayout;

  d->Clear();

  d->set_telescope_type(calin::ix::iact_data::instrument_layout::TelescopeLayout::NO_TELESCOPE);
  d->set_telescope_index(telescope_id);
  calin::math::vector3d_util::dump_as_proto(pos, d->mutable_position());
  d->set_effective_focal_length(param.focal_plane().translation().y());

  auto* c = d->mutable_camera();
  c->set_camera_type(calin::ix::iact_data::instrument_layout::CameraLayout::NO_CAMERA);
  c->set_camera_number(telescope_id);

  c->set_pixel_grid_layout(
    calin::ix::iact_data::instrument_layout::CameraLayout::HEX_GRID);
  c->set_pixel_grid_spacing(param.pixel().spacing());
  c->set_pixel_grid_rotation(param.pixel().grid_rotation());
  c->set_pixel_grid_cos_rotation(std::cos(param.pixel().grid_rotation()/180.0*M_PI));
  c->set_pixel_grid_sin_rotation(std::sin(param.pixel().grid_rotation()/180.0*M_PI));
  c->set_pixel_grid_offset_x(param.pixel().pixel_grid_shift_x());
  c->set_pixel_grid_offset_y(param.pixel().pixel_grid_shift_z());
  c->set_pixel_grid_geometric_area(
    calin::math::hex_array::cell_area(param.pixel().spacing()));

  std::set<unsigned> modules_missing;
  for(auto iid : param.pixel().module_missing_list())
    modules_missing.insert(iid);

  unsigned module_size = param.pixel().hex_module_size();
  bool use_a_config = true;
  if(param.pixel().hex_module_layout_use_b_configuration())use_a_config = false;
  std::map<unsigned, std::pair<unsigned, unsigned> > pixel_hexids;
  unsigned module_index = 0;
  if(modules_missing.find(0) == modules_missing.end())
  {
    auto* m = c->add_module();
    m->set_module_index(0);
    m->set_module_grid_index(0);
    unsigned modchan_id = 0;
    for(auto id : math::hex_array::cluster_hexid_to_member_hexid(0, module_size, use_a_config)) {
      pixel_hexids[id] = std::make_pair(module_index, modchan_id);
      m->add_channels_in_module(id); // temporarily add hexid
      ++modchan_id;
    }
    ++module_index;
  }
  unsigned module_ring_id = 1;
  bool module_ring_in_camera = true;
  if(param.pixel().module_num_hex_rings()==0 and param.focal_plane().camera_diameter()<=0)
    module_ring_in_camera = false;
  while(module_ring_in_camera
        and (param.pixel().module_num_hex_rings()==0
             or module_ring_id<=param.pixel().module_num_hex_rings()))
  {
    module_ring_in_camera = false;
    for(unsigned iseg=0;iseg<6;iseg++)
      for(unsigned irun=0;irun<module_ring_id;irun++)
      {
        unsigned module_id =
            math::hex_array::
            ringid_segid_runid_to_hexid(module_ring_id, iseg, irun);
        if(modules_missing.find(module_id) != modules_missing.end())
          continue;
        int u, v;
        math::hex_array::cluster_hexid_to_center_uv(module_id, module_size, u, v, use_a_config);
        double x, z;
        math::hex_array::uv_to_xy(u, v, x, z);
        if(param.focal_plane().camera_diameter()>0 and
            std::sqrt(x*x+z*z)*param.pixel().spacing()>param.focal_plane().camera_diameter()/2.0)
	        continue; // skip module if position too far out
        module_ring_in_camera = true;
        unsigned modchan_id = 0;
        auto* m = c->add_module();
        m->set_module_index(module_index);
        m->set_module_grid_index(module_id);
        for(auto id : math::hex_array::
            cluster_hexid_to_member_hexid(module_id, module_size, use_a_config)) {
          m->add_channels_in_module(id); // temporarily add hexid
          pixel_hexids[id] = std::make_pair(module_index, modchan_id);
          ++modchan_id;
        }
        ++module_index;
      }
    module_ring_id++;
  }

  math::regular_grid::HexGrid grid(c->pixel_grid_spacing(),
    c->pixel_grid_rotation()/180.0*M_PI,
    c->pixel_grid_offset_x(), c->pixel_grid_offset_y());
  unsigned pixelid = 0;
  for(auto pixel : pixel_hexids)
  {
    unsigned hexid = pixel.first;
    unsigned modid = pixel.second.first;
    unsigned modchanid = pixel.second.second;

    auto* ch = c->add_channel();
    ch->set_channel_index(pixelid);
    ch->set_pixel_index(pixelid);
    ch->set_pixel_grid_index(hexid);
    ch->set_channel_set_index(0);
    ch->set_module_index(modid);
    ch->set_module_channel_index(modchanid);
    double x;
    double y;
    grid.gridid_to_xy(hexid, x, y);
    ch->set_x(x);
    ch->set_y(y);
    ch->set_diameter(c->pixel_grid_spacing());
    ch->set_geometric_area(grid.cell_area(hexid));
    auto nbr = grid.gridid_to_neighbour_gridids(hexid);
    for(auto inbr : nbr) {
      auto inpix = pixel_hexids.find(inbr);
      if(inpix != pixel_hexids.end())
        ch->add_neighbour_channel_indexes(
          std::distance(pixel_hexids.begin(), inpix));
    }
    Eigen::VectorXd xv;
    Eigen::VectorXd yv;
    grid.gridid_to_vertexes_xy(hexid, xv, yv);
    for(unsigned i=0; i<xv.size(); i++)
      ch->add_outline_polygon_vertex_x(xv(i));
    for(unsigned i=0; i<yv.size(); i++)
      ch->add_outline_polygon_vertex_y(yv(i));
    ch->add_outline_polygon_vertex_index(xv.size());
    c->mutable_module(modid)->set_channels_in_module(modchanid, pixelid);
    pixelid++;
  }

  calin::iact_data::instrument_layout::compute_camera_and_module_outlines(c);

  return d;
}
