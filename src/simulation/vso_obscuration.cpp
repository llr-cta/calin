/*

   calin/simulation/vso_obscuration.cpp -- Stephen Fegan -- 2015-11-09

   Classes for obscurations (such as the camera and telescope arms)

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

/*! \file VSOObscuration.hpp
  Obscuration class header file

  \author  Stephen Fegan               \n
           UCLA                        \n
	   sfegan@llr.in2p3.fr         \n

  \date    07/04/2013
  \version 0.1
  \note
*/

#include <math/vector3d_util.hpp>
#include <math/ray.hpp>
#include <simulation/vso_obscuration.hpp>

using namespace calin::math::ray;
using namespace calin::simulation::vs_optics;

VSOObscuration::~VSOObscuration()
{
  // nothing to see here
}

VSOObscuration* VSOObscuration::
create_from_proto(const ix::simulation::vs_optics::VSOObscurationData& d)
{
  if(d.has_disk())return VSODiskObscuration::create_from_proto(d.disk());
  else if(d.has_tube())return VSOTubeObscuration::create_from_proto(d.tube());
  else if(d.has_aligned_box())return VSOAlignedBoxObscuration::create_from_proto(d.aligned_box());
  else {
    throw std::runtime_error("VSOObscuration::create_from_proto: unknown obscuration type");
    return 0;
  }
}

VSODiskObscuration::~VSODiskObscuration()
{
  // nothing to see here
}

bool VSODiskObscuration::doesObscure(const Ray& r_in, Ray& r_out) const
{
  if(fICO && r_in.direction().y()>0)return false;
  r_out = r_in;
  if(r_out.propagate_to_plane(fN, -fD0, false)
     && (r_out.position()-fX0).norm()<=fR) return true;
   return false;
}

calin::ix::simulation::vs_optics::VSOObscurationData* VSODiskObscuration::
dump_as_proto(ix::simulation::vs_optics::VSOObscurationData* d) const
{
  if(d == nullptr)d = new ix::simulation::vs_optics::VSOObscurationData;
  auto* dd = d->mutable_disk();
  calin::math::vector3d_util::dump_as_proto(fX0, dd->mutable_center_pos());
  calin::math::vector3d_util::dump_as_proto(fN, dd->mutable_normal());
  dd->set_diameter(2.0*fR);
  dd->set_incoming_only(fICO);
  return d;
}

VSODiskObscuration* VSODiskObscuration::
create_from_proto(const ix::simulation::vs_optics::VSODiskObscurationData& d)
{
  return new VSODiskObscuration(
    calin::math::vector3d_util::from_proto(d.center_pos()),
    calin::math::vector3d_util::from_proto(d.normal()),
    d.diameter()/2.0, d.incoming_only());
}

VSODiskObscuration* VSODiskObscuration::clone() const
{
  return new VSODiskObscuration(*this);
}


VSOTubeObscuration::~VSOTubeObscuration()
{
  // nothing to see here
}

bool VSOTubeObscuration::doesObscure(const Ray& r_in, Ray& r_out) const
{
  if(fICO && r_in.direction().y()>0)return false;

  r_out = r_in;
  Ray::IPOut ipo;
  ipo = r_out.propagate_to_cylinder(fX1, fN, fR, Ray::IP_NEXT, false);

  if(ipo != Ray::IPO_NONE)
  {
    double Dc = r_out.position().dot(fN);
    if((std::fabs(Dc-fD1)<=fD)&&(std::fabs(Dc-fD2)<=fD))return true;
    if(ipo == Ray::IPO_SECOND)return false;
    ipo = r_out.propagate_to_cylinder(fX1, fN, fR, Ray::IP_LATEST, false);
    if(ipo != Ray::IPO_SECOND)return false;
    Dc = r_out.position().dot(fN);
    if((std::fabs(Dc-fD1)<=fD)&&(std::fabs(Dc-fD2)<=fD))return true;
  }

  return false;
}

calin::ix::simulation::vs_optics::VSOObscurationData* VSOTubeObscuration::
dump_as_proto(ix::simulation::vs_optics::VSOObscurationData* d) const
{
  if(d == nullptr)d = new calin::ix::simulation::vs_optics::VSOObscurationData;
  auto* dd = d->mutable_tube();
  calin::math::vector3d_util::dump_as_proto(fX1, dd->mutable_end1_pos());
  calin::math::vector3d_util::dump_as_proto(fX2, dd->mutable_end2_pos());
  dd->set_diameter(2.0*fR);
  dd->set_incoming_only(fICO);
  return d;
}

VSOTubeObscuration* VSOTubeObscuration::
create_from_proto(const ix::simulation::vs_optics::VSOTubeObscurationData& d)
{
  return new VSOTubeObscuration(
    calin::math::vector3d_util::from_proto(d.end1_pos()),
    calin::math::vector3d_util::from_proto(d.end2_pos()),
    d.diameter()/2.0, d.incoming_only());
}

VSOTubeObscuration* VSOTubeObscuration::clone() const
{
  return new VSOTubeObscuration(*this);
}

VSOAlignedBoxObscuration::~VSOAlignedBoxObscuration()
{
  // nothing to see here
}

bool VSOAlignedBoxObscuration::doesObscure(const Ray& r_in, Ray& r_out) const
{
  // See: https://tavianator.com/fast-branchless-raybounding-box-intersections/
  // and: https://tavianator.com/fast-branchless-raybounding-box-intersections-part-2-nans/

  if(incoming_only_ && r_in.direction().y()>0)return false;

  // Normalized direction vector
  const double vx = 1.0 / r_in.direction().x();
  const double vy = 1.0 / r_in.direction().y();
  const double vz = 1.0 / r_in.direction().z();

  Eigen::Vector3d min_rel = min_corner_ - r_in.position();
  Eigen::Vector3d max_rel = min_corner_ - r_in.position();

  const double tx1 = min_rel.x() * vx;
  const double tx2 = max_rel.x() * vx;
  double tmin = std::min(tx1, tx2);
  double tmax = std::max(tx1, tx2);

  const double ty1 = min_rel.y() * vy;
  const double ty2 = max_rel.y() * vy;
  tmin = std::max(tmin, std::min(std::min(ty1, ty2), tmax));
  tmax = std::min(tmax, std::max(std::max(ty1, ty2), tmin));

  const double tz1 = min_rel.z() * vz;
  const double tz2 = max_rel.z() * vz;
  tmin = std::max(tmin, std::min(std::min(tz1, tz2), tmax));
  tmax = std::min(tmax, std::max(std::max(tz1, tz2), tmin));

  if(tmax > std::max(tmin, 0.0)) {
    r_out = r_in;
    if(tmin > 0)r_out.propagate_dist(tmin);
    return true;
  }

  return false;
}

VSOAlignedBoxObscuration* VSOAlignedBoxObscuration::clone() const
{
  return new VSOAlignedBoxObscuration(*this);
}

calin::ix::simulation::vs_optics::VSOObscurationData*
VSOAlignedBoxObscuration::dump_as_proto(
  calin::ix::simulation::vs_optics::VSOObscurationData* d) const
{
  if(d == nullptr)d = new calin::ix::simulation::vs_optics::VSOObscurationData;
  auto* dd = d->mutable_aligned_box();
  calin::math::vector3d_util::dump_as_proto(max_corner_, dd->mutable_max_corner());
  calin::math::vector3d_util::dump_as_proto(min_corner_, dd->mutable_min_corner());
  dd->set_incoming_only(incoming_only_);
  return d;
}

VSOAlignedBoxObscuration* VSOAlignedBoxObscuration::create_from_proto(
  const ix::simulation::vs_optics::VSOAlignedBoxObscurationData& d)
{
  return new VSOAlignedBoxObscuration(
    calin::math::vector3d_util::from_proto(d.max_corner()),
    calin::math::vector3d_util::from_proto(d.min_corner()),
    d.incoming_only());
}
