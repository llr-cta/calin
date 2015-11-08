/* 

   calin/simulation/vso_mirror.cpp -- Stephen Fegan -- 2015-11-05

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

#include <sstream>
#include <cassert>

#include <vso_mirror.hpp>
#include <vso_telescope.hpp>

using namespace calin::math::vs_physics;
using namespace calin::simulation::vs_optics;

VSOMirror::VSOMirror():
  fTelescope(0), fID(0), fHexID(0), fRemoved(false), fPos(), fAlign(), 
  fFocalLength(0), fSpotSize(0), fDegradingFactor(0), fRotationVector()
{
  // nothing to see here
}

VSOMirror::VSOMirror(const VSOTelescope* T, 
		     unsigned ID, unsigned MHID, bool REM,
		     const Physics::Vec3D& P, const Physics::Vec3D& A, 
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
  Vec3D rrot(0,-fTelescope->reflectorRotation(),0);
  Vec3D align(fAlign);
  align.Rotate(rrot);
  Vec3D normalize = align^Vec3D(0,1,0);
  double sintheta = normalize.Norm();
  if(sintheta == 0)
    {
      fRotationVector = rrot;
    }
  else 
    {
      double costheta = align*Vec3D(0,1,0);
      normalize *= atan(sintheta/costheta)/sintheta;
      fRotationVector = rrot & normalize;
    }
}

void VSOMirror::reflectorToMirror(Physics::Particle& p) const
{
  //  std::cerr << fRotationVector << std::endl;
  // First: Translate from reflector to mirror mount point
  p.TranslateOrigin(Vec4D(0,fPos));
  // Second: Rotate coordinate system so mirror normal is y-axis
  p.Rotate(fRotationVector);
  // Third: Fix parity
  // if(fTelescope->mirrorParity())mom3.x=-mom3.x;
  //  std::cout << "A: " << p.Momenta().r << ' ' << fRotationVector << ' ';
  //  std::cout << p.Momenta().r << std::endl;
#warning clean me up
}

void VSOMirror::mirrorToReflector(Physics::Particle& p) const
{
#warning Fix parity
  // First: Fix parity
  // if(fTelescope->mirrorParity())v.x=-v.x;
  // Second: Rotate coordinate system back to the reflector
  p.Rotate(-fRotationVector);
  // Third: Translate from mirror mount point to the reflector
  p.TranslateOrigin(Vec4D(0,-fPos));
}

Physics::Vec3D VSOMirror::
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
  return Vec3D(r*x, y, r*z);
}

Physics::Vec3D VSOMirror::
cornerInReflectorCoords(unsigned icorner, double aperture) const
{
  Vec3D r = cornerInMirrorCoords(icorner, aperture);
  r.Rotate(-fRotationVector);
  r += fPos;
  return r;
}

void VSOMirror::dumpShort(std::ostream& stream) const
{
  stream
    << "MIRROR "
    << VSDataConverter::toString(fTelescope->id()) << ' '
    << VSDataConverter::toString(fID) << ' '
    << VSDataConverter::toString(fHexID) << ' '
    << VSDataConverter::toString(fRemoved) << ' '
    << VSDataConverter::toString(fPos.x) << ' '

    << VSDataConverter::toString(fPos.y) << ' '
    << VSDataConverter::toString(fPos.z) << ' '
    << VSDataConverter::toString(fAlign.x) << ' '
    << VSDataConverter::toString(fAlign.y) << ' '
    << VSDataConverter::toString(fAlign.z) << ' '

    << VSDataConverter::toString(fFocalLength) << ' '
    << VSDataConverter::toString(fSpotSize) << ' '
    << VSDataConverter::toString(fDegradingFactor) << std::endl;
}

void VSOMirror::dump(std::ostream& stream, unsigned l) const
{
  stream
    << FDAV("Telescope ID", fTelescope->id(), "", 30, l) << std::endl
    << FDAV("Mirror ID", fID, "", 30, l) << std::endl
    << FDAV("Mirror Hex ID", fHexID, "", 30, l) << std::endl
    << FDAV("Removed", fRemoved, "", 30, l) << std::endl
    << FDAV("Position X", fPos.x, "", 30, l) << std::endl

    << FDAV("Position Y", fPos.y, "", 30, l) << std::endl
    << FDAV("Position Z", fPos.z, "", 30, l) << std::endl
    << FDAV("AlignmentX", fAlign.x, "", 30, l) << std::endl
    << FDAV("AlignmentY", fAlign.y, "", 30, l) << std::endl
    << FDAV("AlignmentZ", fAlign.z, "", 30, l) << std::endl

    << FDAV("FocalLength", fFocalLength, "", 30, l) << std::endl
    << FDAV("SpotSize", fSpotSize, "", 30, l) << std::endl
    << FDAV("DegradingFactor", fDegradingFactor, "", 30, l) << std::endl;
}

VSOMirror* VSOMirror::createFromShortDump(std::istream& stream,
					  const VSOTelescope* T)
{
  std::string line;
  std::getline(stream,line);
  if(line.empty())return 0;

  std::istringstream linestream(line);

  VSOMirror* mirror = new VSOMirror;
  mirror->fTelescope=T;

  std::string keyword;
  linestream >> keyword;
  assert(keyword==std::string("MIRROR"));

  unsigned TID;

  linestream
    >> TID
    >> mirror->fID
    >> mirror->fHexID
    >> mirror->fRemoved
    >> mirror->fPos.x

    >> mirror->fPos.y
    >> mirror->fPos.z
    >> mirror->fAlign.x
    >> mirror->fAlign.y
    >> mirror->fAlign.z

    >> mirror->fFocalLength
    >> mirror->fSpotSize
    >> mirror->fDegradingFactor;

  if(!linestream)
    {
      delete mirror;
      return 0;
    }

  mirror->calculateRotationVector();

  return mirror;
}
