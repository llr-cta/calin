/*

   calin/simulation/vso_pixel.cpp -- Stephen Fegan -- 2015-11-10

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

#include <sstream>
#include <cassert>

#include <simulation/vso_pixel.hpp>
#include <simulation/vso_telescope.hpp>

using namespace calin::math::vs_physics;
using namespace calin::simulation::vs_optics;

VSOPixel::VSOPixel():
  fTelescope(0), fID(0), fHexID(0), fRemoved(false), fPos()
{
  // nothing to see here
}

VSOPixel::VSOPixel(const VSOTelescope* T,
		   unsigned PID, unsigned PHID, bool REM,
		   const Vec3D& P):
  fTelescope(T), fID(PID), fHexID(PHID), fRemoved(REM), fPos(P)
{
  // nothing to see here
}

VSOPixel::~VSOPixel()
{
  // nothing to see here
}

Vec3D VSOPixel::incomingSkyVectorAtZenith(double plate_scale) const
{
  double x = fPos.x*plate_scale;
  double y = -fPos.z*plate_scale;
  double z = -(fTelescope->focalPlanePosition().y+fPos.y);

  Vec3D p(x,y,z);
  p/=p.Norm();

  return p;
}

calin::ix::simulation::vs_optics::VSOPixelData*
VSOPixel::dump_as_proto(calin::ix::simulation::vs_optics::VSOPixelData* d) const
{
  if(d == nullptr)d = new calin::ix::simulation::vs_optics::VSOPixelData;
  d->set_id(fID);
  d->set_hex_id(fHexID);
  d->set_removed(fRemoved);
  fPos.dump_as_proto(d->mutable_pos());
  return d;
}

VSOPixel*
VSOPixel::create_from_proto(const ix::simulation::vs_optics::VSOPixelData& d,
                            const VSOTelescope* T)
{
  VSOPixel* pixel = new VSOPixel;

  pixel->fTelescope       = T;
  pixel->fID              = d.id();
  pixel->fHexID           = d.hex_id();
  pixel->fRemoved         = d.removed();
  pixel->fPos.set_from_proto(d.pos());

  return pixel;
}

#if 0
void VSOPixel::dumpShort(std::ostream& stream) const
{
  stream
    << "PIXEL "

    << VSDataConverter::toString(fTelescope->id()) << ' '
    << VSDataConverter::toString(fID) << ' '
    << VSDataConverter::toString(fHexID) << ' '
    << VSDataConverter::toString(fRemoved) << ' '
    << VSDataConverter::toString(fPos.x) << ' '

    << VSDataConverter::toString(fPos.y) << ' '
    << VSDataConverter::toString(fPos.z) << std::endl;
}

void VSOPixel::dump(std::ostream& stream, unsigned l) const
{
  stream
    << FDAV("Telescope ID", fTelescope->id(), "", 30, l) << std::endl
    << FDAV("Pixel ID", fID, "", 30, l) << std::endl
    << FDAV("Pixel Hex ID", fHexID, "", 30, l) << std::endl
    << FDAV("Removed", fRemoved, "", 30, l) << std::endl
    << FDAV("Position X", fPos.x, "", 30, l) << std::endl

    << FDAV("Position Y", fPos.y, "", 30, l) << std::endl
    << FDAV("Position Z", fPos.z, "", 30, l) << std::endl;
}

VSOPixel* VSOPixel::createFromShortDump(std::istream& stream,
					  const VSOTelescope* T)
{
  std::string line;
  std::getline(stream,line);
  if(line.empty())return 0;

  std::istringstream linestream(line);

  VSOPixel* pixel = new VSOPixel;
  pixel->fTelescope=T;

  std::string keyword;
  linestream >> keyword;
  assert(keyword==std::string("PIXEL"));

  unsigned TID;

  linestream
    >> TID
    >> pixel->fID
    >> pixel->fHexID
    >> pixel->fRemoved
    >> pixel->fPos.x

    >> pixel->fPos.y
    >> pixel->fPos.z;

  if(!linestream)
    {
      delete pixel;
      return 0;
    }

  return pixel;
}
#endif
