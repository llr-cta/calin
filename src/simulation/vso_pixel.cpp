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
#include <math/vector3d_util.hpp>

using namespace calin::simulation::vs_optics;

VSOPixel::VSOPixel():
  fTelescope(0), fID(0), fHexID(0), fRemoved(false), fPos(0,0,0)
{
  // nothing to see here
}

VSOPixel::VSOPixel(const VSOTelescope* T,
		   unsigned PID, unsigned PHID, bool REM,
		   const Eigen::Vector3d& P):
  fTelescope(T), fID(PID), fHexID(PHID), fRemoved(REM), fPos(P)
{
  // nothing to see here
}

VSOPixel::~VSOPixel()
{
  // nothing to see here
}

Eigen::Vector3d VSOPixel::incomingSkyVectorAtZenith(double plate_scale) const
{
  double x = fPos.x()*plate_scale;
  double y = -fPos.z()*plate_scale;
  double z = -(fTelescope->focalPlanePosition().y()+fPos.y());

  Eigen::Vector3d p(x,y,z);
  p.normalize();

  return p;
}

calin::ix::simulation::vs_optics::VSOPixelData*
VSOPixel::dump_as_proto(calin::ix::simulation::vs_optics::VSOPixelData* d) const
{
  if(d == nullptr)d = new calin::ix::simulation::vs_optics::VSOPixelData;
  d->set_id(fID);
  d->set_hex_id(fHexID);
  d->set_removed(fRemoved);
  calin::math::vector3d_util::dump_as_proto(fPos,d->mutable_pos());
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
  calin::math::vector3d_util::set_from_proto(pixel->fPos,d.pos());

  return pixel;
}
