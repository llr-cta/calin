/*

   calin/simulation/vso_obscuration.cpp -- Stephen Fegan -- 2015-11-09

   Classes for obscurations (such as the camera and telescope arms)

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
#include <math/geometry.hpp>
#include <math/special.hpp>
#include <math/hex_array.hpp>

using namespace calin::math::ray;
using namespace calin::simulation::vs_optics;
using calin::math::special::SQR;
using namespace calin::util::log;

VSOObscuration::~VSOObscuration()
{
  // nothing to see here
}

VSOObscuration* VSOObscuration::
create_from_proto(const ix::simulation::vs_optics::VSOObscurationData& d)
{
  VSOObscuration* obs = nullptr;
  if(d.has_disk()) {
    obs = VSODiskObscuration::create_from_proto(d.disk());
  } else if(d.has_tube()) {
    obs = VSOTubeObscuration::create_from_proto(d.tube());
  } else if(d.has_aligned_box()) {
    obs = VSOAlignedBoxObscuration::create_from_proto(d.aligned_box());
  } else if(d.has_rectangular_aperture()) {
    obs = VSOAlignedRectangularAperture::create_from_proto(d.rectangular_aperture());
  } else if(d.has_circular_aperture()) {
    obs = VSOAlignedCircularAperture::create_from_proto(d.circular_aperture());
  } else if(d.has_hexagonal_aperture()) {
    obs = VSOAlignedHexagonalAperture::create_from_proto(d.hexagonal_aperture());
  } else if(d.has_tile_aperture()) {
    obs = VSOAlignedTileAperture::create_from_proto(d.tile_aperture());
  } else if(d.has_box_collection()) {
    obs = VSOBoxCollectionObscuration::create_from_proto(d.box_collection());
  } else {
    throw std::runtime_error("VSOObscuration::create_from_proto: unknown obscuration type");
    return 0;
  }
  return obs;
}

// *****************************************************************************
// *****************************************************************************
//
// VSODiskObscuration
//
// *****************************************************************************
// *****************************************************************************

VSODiskObscuration::~VSODiskObscuration()
{
  // nothing to see here
}

bool VSODiskObscuration::doesObscure(const Ray& r_in, Ray& r_out, double n) const
{
  r_out = r_in;
  if(r_out.propagate_to_plane(fN, -fD0, false, n)
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
  dd->set_identification(identification_);
  return d;
}

VSODiskObscuration* VSODiskObscuration::
create_from_proto(const ix::simulation::vs_optics::VSODiskObscurationData& d)
{
  return new VSODiskObscuration(
    calin::math::vector3d_util::from_proto(d.center_pos()),
    calin::math::vector3d_util::from_proto(d.normal()),
    d.diameter()/2.0, d.identification());
}

VSODiskObscuration* VSODiskObscuration::clone() const
{
  return new VSODiskObscuration(*this);
}

// *****************************************************************************
// *****************************************************************************
//
// VSOTubeObscuration
//
// *****************************************************************************
// *****************************************************************************

VSOTubeObscuration::~VSOTubeObscuration()
{
  // nothing to see here
}

bool VSOTubeObscuration::doesObscure(const Ray& r_in, Ray& r_out, double n) const
{
  r_out = r_in;
  Ray::IPOut ipo;
  ipo = r_out.propagate_to_cylinder(fX1, fN, fR, Ray::IP_NEXT, false, n);

  if(ipo != Ray::IPO_NONE)
  {
    double Dc = r_out.position().dot(fN);
    if((std::fabs(Dc-fD1)<=fD)&&(std::fabs(Dc-fD2)<=fD))return true;
    if(ipo == Ray::IPO_SECOND)return false;
    ipo = r_out.propagate_to_cylinder(fX1, fN, fR, Ray::IP_LATEST, false, n);
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
  dd->set_identification(identification_);
  return d;
}

VSOTubeObscuration* VSOTubeObscuration::
create_from_proto(const ix::simulation::vs_optics::VSOTubeObscurationData& d)
{
  return new VSOTubeObscuration(
    calin::math::vector3d_util::from_proto(d.end1_pos()),
    calin::math::vector3d_util::from_proto(d.end2_pos()),
    d.diameter()/2.0, d.identification());
}

VSOTubeObscuration* VSOTubeObscuration::clone() const
{
  return new VSOTubeObscuration(*this);
}

// *****************************************************************************
// *****************************************************************************
//
// VSOAlignedBoxObscuration
//
// *****************************************************************************
// *****************************************************************************

VSOAlignedBoxObscuration::~VSOAlignedBoxObscuration()
{
  // nothing to see here
}

bool VSOAlignedBoxObscuration::
doesObscure(const Ray& r_in, Ray& r_out, double n) const
{
  double tmin;
  double tmax;
  if(calin::math::geometry::box_has_future_intersection(tmin, tmax,
    min_corner_, max_corner_, r_in.position(), r_in.direction()))
  {
    r_out = r_in;
    if(tmin > 0)r_out.propagate_dist(tmin, n);
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
  dd->set_identification(identification_);
  return d;
}

VSOAlignedBoxObscuration* VSOAlignedBoxObscuration::create_from_proto(
  const ix::simulation::vs_optics::VSOAlignedBoxObscurationData& d)
{
  return new VSOAlignedBoxObscuration(
    calin::math::vector3d_util::from_proto(d.max_corner()),
    calin::math::vector3d_util::from_proto(d.min_corner()),
    d.identification());
}

// *****************************************************************************
// *****************************************************************************
//
// VSOAlignedRectangularAperture
//
// *****************************************************************************
// *****************************************************************************

VSOAlignedRectangularAperture::~VSOAlignedRectangularAperture()
{
  // nothing to see here
}


bool VSOAlignedRectangularAperture::
doesObscure(const calin::math::ray::Ray& r_in, calin::math::ray::Ray& r_out, double n) const
{
  r_out = r_in;
  if(r_out.propagate_to_y_plane(-center_.y(), false, n))
  {
    const double De =
      std::max(fabs(r_out.x()-center_.x())-flat_to_flat_x_2_,
               fabs(r_out.z()-center_.z())-flat_to_flat_z_2_);
    if(inverted_) {
      return De <= 0;
    } else {
      return De > 0;
    }
  }

  return false;
}

VSOAlignedRectangularAperture* VSOAlignedRectangularAperture::clone() const
{
  return new VSOAlignedRectangularAperture(center_,2*flat_to_flat_x_2_,2*flat_to_flat_z_2_);
}

calin::ix::simulation::vs_optics::VSOObscurationData*
VSOAlignedRectangularAperture::dump_as_proto(
  calin::ix::simulation::vs_optics::VSOObscurationData* d) const
{
  if(d == nullptr)d = new calin::ix::simulation::vs_optics::VSOObscurationData;
  auto* dd = d->mutable_rectangular_aperture();
  calin::math::vector3d_util::dump_as_proto(center_, dd->mutable_center_pos());
  dd->set_flat_to_flat_x(2*flat_to_flat_x_2_);
  dd->set_flat_to_flat_z(2*flat_to_flat_z_2_);
  dd->set_invert(inverted_);
  dd->set_identification(identification_);
  return d;
}

VSOAlignedRectangularAperture* VSOAlignedRectangularAperture::create_from_proto(
  const calin::ix::simulation::vs_optics::VSOAlignedRectangularApertureData& d)
{
  return new VSOAlignedRectangularAperture(
    calin::math::vector3d_util::from_proto(d.center_pos()),
    d.flat_to_flat_x(), d.flat_to_flat_z(), d.invert(), d.identification());
}

// *****************************************************************************
// *****************************************************************************
//
// VSOAlignedHexagonalAperture
//
// *****************************************************************************
// *****************************************************************************

VSOAlignedHexagonalAperture::~VSOAlignedHexagonalAperture()
{
  // nothing to see here
}

bool VSOAlignedHexagonalAperture::
doesObscure(const calin::math::ray::Ray& r_in, calin::math::ray::Ray& r_out, double n) const
{
  r_out = r_in;
  if(r_out.propagate_to_y_plane(-center_.y(), false, n))
  {
    constexpr double cos60 = 0.5;
    constexpr double sin60 = 0.5*CALIN_HEX_ARRAY_SQRT3;
    const double x = r_out.x()-center_.x();
    const double z = r_out.z()-center_.z();
    const double x0 = std::fabs(x);
    const double xp = std::fabs(cos60*x - sin60*z);
    const double xn = std::fabs(cos60*x + sin60*z);
    const double De = std::max(x0,std::max(xp,xn)) - flat_to_flat_2_;
    return De > 0;
  }

  return false;
}

VSOAlignedHexagonalAperture* VSOAlignedHexagonalAperture::clone() const
{
  return new VSOAlignedHexagonalAperture(center_,2*flat_to_flat_2_);
}

calin::ix::simulation::vs_optics::VSOObscurationData*
VSOAlignedHexagonalAperture::dump_as_proto(
  calin::ix::simulation::vs_optics::VSOObscurationData* d) const
{
  if(d == nullptr)d = new calin::ix::simulation::vs_optics::VSOObscurationData;
  auto* dd = d->mutable_hexagonal_aperture();
  calin::math::vector3d_util::dump_as_proto(center_, dd->mutable_center_pos());
  dd->set_flat_to_flat(2*flat_to_flat_2_);
  dd->set_identification(identification_);
  return d;
}

VSOAlignedHexagonalAperture* VSOAlignedHexagonalAperture::create_from_proto(
  const calin::ix::simulation::vs_optics::VSOAlignedHexagonalApertureData& d)
{
  return new VSOAlignedHexagonalAperture(
    calin::math::vector3d_util::from_proto(d.center_pos()), d.flat_to_flat(),
    d.identification());
}

// *****************************************************************************
// *****************************************************************************
//
// VSOAlignedCircularAperture
//
// *****************************************************************************
// *****************************************************************************

VSOAlignedCircularAperture::~VSOAlignedCircularAperture()
{
  // nothing to see here
}

bool VSOAlignedCircularAperture::
doesObscure(const calin::math::ray::Ray& r_in, calin::math::ray::Ray& r_out, double n) const
{
  r_out = r_in;
  if(r_out.propagate_to_y_plane(-center_.y(), false, n))
  {
    const double D2e =
      SQR(r_out.x()-center_.x())+SQR(r_out.z()-center_.z())-radius_sq_;
    if(inverted_) {
      return D2e <= 0;
    } else {
      return D2e > 0;
    }
  }

  return false;
}

VSOAlignedCircularAperture* VSOAlignedCircularAperture::clone() const
{
  return new VSOAlignedCircularAperture(center_,2*sqrt(radius_sq_));
}

calin::ix::simulation::vs_optics::VSOObscurationData*
VSOAlignedCircularAperture::dump_as_proto(
  calin::ix::simulation::vs_optics::VSOObscurationData* d) const
{
  if(d == nullptr)d = new calin::ix::simulation::vs_optics::VSOObscurationData;
  auto* dd = d->mutable_circular_aperture();
  calin::math::vector3d_util::dump_as_proto(center_, dd->mutable_center_pos());
  dd->set_diameter(2*sqrt(radius_sq_));
  dd->set_invert(inverted_);
  dd->set_identification(identification_);
  return d;
}

VSOAlignedCircularAperture* VSOAlignedCircularAperture::create_from_proto(
  const calin::ix::simulation::vs_optics::VSOAlignedCircularApertureData& d)
{
  return new VSOAlignedCircularAperture(
    calin::math::vector3d_util::from_proto(d.center_pos()),d.diameter(),d.invert(),
    d.identification());
}

// *****************************************************************************
// *****************************************************************************
//
// VSOAlignedTileAperture
//
// *****************************************************************************
// *****************************************************************************

VSOAlignedTileAperture::~VSOAlignedTileAperture()
{
  // nothing to see here
}

bool VSOAlignedTileAperture::
doesObscure(const calin::math::ray::Ray& r_in, calin::math::ray::Ray& r_out, double n) const
{
  r_out = r_in;
  if(r_out.propagate_to_y_plane(-center_.y(), false, n))
  {
    double x = r_out.x() - center_.x();
    double z = r_out.z() - center_.z();

    x -= edge_x_;
    z -= edge_z_;
    x *= pitch_x_inv_;
    z *= pitch_z_inv_;
    x -= std::round(x);
    z -= std::round(z);

    const double Dm =
      std::min(std::fabs(x)-opaque_frac_x_2_, std::fabs(z)-opaque_frac_z_2_);
    return (Dm < 0);
  }

  return false;
}

VSOAlignedTileAperture* VSOAlignedTileAperture::clone() const
{
  return new VSOAlignedTileAperture(*this);
  // center_,2*flat_to_flat_x_2_,2*flat_to_flat_z_2_);
  //   1.0/pitch_x_inv_, 1.0/pitch_z_inv_
  //     double center_x, double center_z,
  //     double support_width_x, double support_width_z):
  //   VSOObscuration(), center_(center),
  //   flat_to_flat_x_2_(0.5*flat_to_flat_x), flat_to_flat_z_2_(0.5*flat_to_flat_z),
  //   pitch_x_inv_(1.0/pitch_x), pitch_z_inv_(1.0/pitch_z),
  //   edge_x_(center_x - 0.5*pitch_x),
  //   edge_z_(center_z - 0.5*pitch_z),
  //   opaque_frac_x_2_(0.5*support_width_x*pitch_x_inv_),
  //   opaque_frac_z_2_(0.5*support_width_z*pitch_z_inv_),
}

calin::ix::simulation::vs_optics::VSOObscurationData*
VSOAlignedTileAperture::dump_as_proto(
  calin::ix::simulation::vs_optics::VSOObscurationData* d) const
{
  if(d == nullptr)d = new calin::ix::simulation::vs_optics::VSOObscurationData;
  auto* dd = d->mutable_tile_aperture();
  calin::math::vector3d_util::dump_as_proto(center_, dd->mutable_center_pos());
  dd->set_pitch_x(1.0/pitch_x_inv_);
  dd->set_pitch_z(1.0/pitch_z_inv_);
  dd->set_center_x(edge_x_ + 0.5*dd->pitch_x());
  dd->set_center_z(edge_z_ + 0.5*dd->pitch_z());
  dd->set_support_width_x(2.0*opaque_frac_x_2_*dd->pitch_x());
  dd->set_support_width_z(2.0*opaque_frac_z_2_*dd->pitch_z());
  dd->set_identification(identification_);
  return d;
}

VSOAlignedTileAperture* VSOAlignedTileAperture::create_from_proto(
  const calin::ix::simulation::vs_optics::VSOAlignedTileApertureData& d)
{
  return new VSOAlignedTileAperture(
    calin::math::vector3d_util::from_proto(d.center_pos()),
    d.pitch_x(), d.pitch_z(),
    d.center_x(), d.center_z(),
    d.support_width_x(), d.support_width_z(), d.identification());
}

// *****************************************************************************
// *****************************************************************************
//
// VSOBoxCollectionObscuration
//
// *****************************************************************************
// *****************************************************************************

VSOBoxCollectionObscuration::
VSOBoxCollectionObscuration(const std::vector<VSOTubeObscuration*>& tubes,
    const std::string& identification, bool adopt_obscurations):
  VSOObscuration(identification),
  tubes_(tubes), adopt_obscurations_(adopt_obscurations)
{
  constexpr double inf = std::numeric_limits<double>::infinity();
  double min_vol = inf;
  double min_area;
  double min_theta;
  for(double theta = 0; theta<360.0; theta+=1.0) {
    double st = std::sin(theta/180.0*M_PI);
    double ct = std::cos(theta/180.0*M_PI);
    Eigen::Vector3d min_corner(inf,inf,inf);
    Eigen::Vector3d max_corner(-inf,-inf,-inf);
    for(const auto* tube : tubes) {
      Eigen::Vector3d x1(ct*tube->end1_pos().x()+st*tube->end1_pos().z(),
        tube->end1_pos().y(), ct*tube->end1_pos().z()-st*tube->end1_pos().x());
      Eigen::Vector3d x2(ct*tube->end2_pos().x()+st*tube->end2_pos().z(),
        tube->end2_pos().y(), ct*tube->end2_pos().z()-st*tube->end2_pos().x());
      min_corner.x() = std::min(min_corner.x(), std::min(x1.x(),x2.x())-tube->radius());
      max_corner.x() = std::max(max_corner.x(), std::max(x1.x(),x2.x())+tube->radius());
      min_corner.y() = std::min(min_corner.y(), std::min(x1.y(),x2.y())-tube->radius());
      max_corner.y() = std::max(max_corner.y(), std::max(x1.y(),x2.y())+tube->radius());
      min_corner.z() = std::min(min_corner.z(), std::min(x1.z(),x2.z())-tube->radius());
      max_corner.z() = std::max(max_corner.z(), std::max(x1.z(),x2.z())+tube->radius());
    }
    Eigen::Vector3d dx = max_corner-min_corner;
    double vol = dx.x()*dx.y()*dx.z();
    if(vol < min_vol) {
      min_corner_ = min_corner;
      max_corner_ = max_corner;
      crot_ = ct;
      srot_ = st;
      min_theta = theta;
      min_area = dx.x()*dx.z();
      min_vol = vol;
    }
  }
#if 0
  LOG(INFO) << identification_ << " : ntube=" << tubes.size() << "\n  min_vol="
    << min_vol <<  " cm^3, min_area= " << min_area << " cm^2, \n  theta="
    << min_theta << "\n  min_corner=" << min_corner_.transpose() << "\n  max_corner=" << max_corner_.transpose();
#endif
}

VSOBoxCollectionObscuration::~VSOBoxCollectionObscuration()
{
  if(adopt_obscurations_) {
    for(auto* obs : tubes_)delete obs;
  }
}

bool VSOBoxCollectionObscuration::doesObscure(const calin::math::ray::Ray& p_in,
  calin::math::ray::Ray& p_out, double n) const
{
  Eigen::Vector3d pos(crot_*p_in.x()+srot_*p_in.z(), p_in.y(), crot_*p_in.z()-srot_*p_in.x());
  Eigen::Vector3d dir(crot_*p_in.ux()+srot_*p_in.uz(), p_in.uy(), crot_*p_in.uz()-srot_*p_in.ux());
  bool hits_box = calin::math::geometry::box_has_future_intersection(
    min_corner_, max_corner_, pos, dir);
  if(not hits_box) {
    // fast exit if outer box not hit
    return false;
  }

  bool hits_tube = false;
  double obs_time  = std::numeric_limits<double>::infinity();
  for(const auto* tube : tubes_) {
    math::ray::Ray ray_out;
    if(tube->doesObscure(p_in, ray_out, n)) {
      hits_tube = true;
      if(ray_out.ct() < obs_time) {
        p_out = ray_out;
        obs_time = ray_out.ct();
      }
    }
  }
  return hits_tube;
}

VSOBoxCollectionObscuration* VSOBoxCollectionObscuration::clone() const
{
  std::vector<VSOTubeObscuration*> new_tubes;
  for(auto* obs : tubes_) {
    new_tubes.push_back(obs->clone());
  }
  return new VSOBoxCollectionObscuration(new_tubes, identification_, true);
}

calin::ix::simulation::vs_optics::VSOObscurationData*
VSOBoxCollectionObscuration::dump_as_proto(
  calin::ix::simulation::vs_optics::VSOObscurationData* d) const
{
  if(d == nullptr)d = new calin::ix::simulation::vs_optics::VSOObscurationData;
  auto* dd = d->mutable_box_collection();
  for(auto* obs : tubes_) {
    calin::ix::simulation::vs_optics::VSOObscurationData obs_data;
    obs->dump_as_proto(&obs_data);
    dd->CopyFrom(obs_data.tube());
  }
  dd->set_identification(identification_);
  return d;
}

VSOBoxCollectionObscuration* VSOBoxCollectionObscuration::create_from_proto(
    const ix::simulation::vs_optics::VSOBoxCollectionObscurationData& d)
{
  std::vector<VSOTubeObscuration*> new_tubes;
  for(const auto& obs_data : d.tube_obscuration()) {
    new_tubes.push_back(VSOTubeObscuration::create_from_proto(obs_data));
  }
  return new VSOBoxCollectionObscuration(new_tubes, d.identification(), true);
}
