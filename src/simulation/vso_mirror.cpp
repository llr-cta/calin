/*

   calin/simulation/vso_mirror.cpp -- Stephen Fegan -- 2015-11-05

   Class for mirror data and functionality. This code is derived from
   raytracing code largely implemented by the author at UCLA in
   2004. It is not complient with the calin coding style.

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

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
#include <cassert>

#include <math/vector3d_util.hpp>
#include <simulation/vso_mirror.hpp>
#include <simulation/vso_telescope.hpp>

using namespace calin::simulation::vs_optics;

VSOMirror::VSOMirror():
  fTelescope(0), fID(0), fHexID(0), fRemoved(false), fPos(), fAlign(),
  fFocalLength(0), fSpotSize(0), fDegradingFactor(0), fRotationVector()
{
  // nothing to see here
}

VSOMirror::VSOMirror(const VSOTelescope* T,
		     unsigned ID, unsigned MHID, bool REM,
		     const Eigen::Vector3d& P, const Eigen::Vector3d& A,
		     double FL, double SS, double DF):
  fTelescope(T), fID(ID), fHexID(MHID), fRemoved(REM), fPos(P), fAlign(A),
  fFocalLength(FL), fSpotSize(SS), fDegradingFactor(DF), fRotationVector()
{
  calculateRotationVector();
}

VSOMirror::~VSOMirror()
{
  // nothing to see here
}

void VSOMirror::calculateRotationVector()
{
  if(fTelescope == nullptr) {
    fRotationVector = Eigen::Matrix3d::Identity();
    return;
  }

  Eigen::Matrix3d rrot =
    Eigen::AngleAxisd(-fTelescope->reflectorRotation(), Eigen::Vector3d::UnitY()).toRotationMatrix();
  Eigen::Vector3d align = rrot * fAlign;
  Eigen::Vector3d normalize = align.cross(Eigen::Vector3d::UnitY());
  double sintheta = normalize.norm();
  if(sintheta == 0) {
    fRotationVector = rrot;
  } else {
    double costheta = align.y(); // align.dot(Eigen::Vector3d::UnitY());
    normalize.normalize();
    fRotationVector =
      Eigen::AngleAxisd(std::atan(sintheta/costheta), normalize) * rrot;
  }

#if 0
  Eigen::Vector3d rrot(0,-fTelescope->reflectorRotation(),0);
  Eigen::Vector3d align(fAlign);
  align.Rotate(rrot);
  Eigen::Vector3d normalize = align^Eigen::Vector3d(0,1,0);
  double sintheta = normalize.Norm();
  if(sintheta == 0)
    {
      fRotationVector = rrot;
    }
  else
    {
      double costheta = align*Eigen::Vector3d(0,1,0);
      normalize *= atan(sintheta/costheta)/sintheta;
      fRotationVector = rrot & normalize;
    }
#endif
}

void VSOMirror::reflectorToMirror(math::ray::Ray& r) const
{
  //  std::cerr << fRotationVector << std::endl;
  // First: Translate from reflector to mirror mount point
  r.translate_origin(fPos);
  // Second: Rotate coordinate system so mirror normal is y-axis
  r.rotate(fRotationVector);
}

void VSOMirror::mirrorToReflector(math::ray::Ray& r) const
{
  // First: Rotate coordinate system back to the reflector
  r.derotate(fRotationVector);
  // Third: Translate from mirror mount point to the reflector
  r.untranslate_origin(fPos);
}

Eigen::Vector3d VSOMirror::
cornerInMirrorCoords(unsigned icorner, double aperture) const
{
  static const double sin30 = 1.0/2.0;
  static const double cos30 = sqrt(3.0)/2.0;

  double r = aperture/cos30/2.0;
  double x;
  double z;
  switch(icorner%6)
  {
  case 0: x =  cos30; z =  sin30; break;
  case 1: x =    0.0; z =    1.0; break;
  case 2: x = -cos30; z =  sin30; break;
  case 3: x = -cos30; z = -sin30; break;
  case 4: x =    0.0; z =   -1.0; break;
  case 5: x =  cos30; z = -sin30; break;
  }
  double y = 2.0*fFocalLength-std::sqrt(4.0*fFocalLength*fFocalLength-r*r);
  return Eigen::Vector3d(r*x, y, r*z);
}

Eigen::Vector3d VSOMirror::
cornerInReflectorCoords(unsigned icorner, double aperture) const
{
  Eigen::Vector3d r = cornerInMirrorCoords(icorner, aperture);
  r = fRotationVector.transpose() * r;
  r += fPos;
  return r;
}

calin::ix::simulation::vs_optics::VSOMirrorData*
VSOMirror::dump_as_proto(ix::simulation::vs_optics::VSOMirrorData* d) const
{
  if(d == nullptr)d = new calin::ix::simulation::vs_optics::VSOMirrorData;
  d->set_id(fID);
  d->set_hex_id(fHexID);
  d->set_removed(fRemoved);
  calin::math::vector3d_util::dump_as_proto(fPos,d->mutable_pos());
  calin::math::vector3d_util::dump_as_proto(fAlign,d->mutable_align());
  d->set_focal_length(fFocalLength);
  d->set_spot_size(fSpotSize);
  d->set_degrading_factor(fDegradingFactor);
  return d;
}

VSOMirror* VSOMirror::
create_from_proto(const ix::simulation::vs_optics::VSOMirrorData& d,
                  const VSOTelescope* T)
{
  VSOMirror* mirror = new VSOMirror;

  mirror->fTelescope       = T;
  mirror->fID              = d.id();
  mirror->fHexID           = d.hex_id();
  mirror->fRemoved         = d.removed();
  calin::math::vector3d_util::set_from_proto(mirror->fPos, d.pos());
  calin::math::vector3d_util::set_from_proto(mirror->fAlign, d.align());
  mirror->fFocalLength     = d.focal_length();
  mirror->fSpotSize        = d.spot_size();
  mirror->fDegradingFactor = d.degrading_factor();

  mirror->calculateRotationVector();

  return mirror;
}
