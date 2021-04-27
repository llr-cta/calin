/*

   calin/simulation/vso_obscuration.hpp -- Stephen Fegan -- 2015-11-09

   Class for telescope data and functionality. This code is derived
   from raytracing code largely implemented by the author at UCLA in
   2004. It is not complient with the calin coding style.

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

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <math/ray.hpp>
#include <simulation/vs_optics.pb.h>

namespace calin { namespace simulation { namespace vs_optics {

class VSOObscuration
{
 public:
  virtual ~VSOObscuration();
  virtual bool doesObscure(const calin::math::ray::Ray& p_in,
                           calin::math::ray::Ray& p_out, double n) const = 0;

  virtual VSOObscuration* clone() const = 0;

#ifndef SWIG
  virtual calin::ix::simulation::vs_optics::VSOObscurationData* dump_as_proto(
    calin::ix::simulation::vs_optics::VSOObscurationData* d = nullptr) const = 0;
#else
  virtual calin::ix::simulation::vs_optics::VSOObscurationData* dump_as_proto() const = 0;
  virtual void dump_as_proto(calin::ix::simulation::vs_optics::VSOObscurationData* d) const = 0;
#endif

  static VSOObscuration*
  create_from_proto(const ix::simulation::vs_optics::VSOObscurationData& d);
};

class VSODiskObscuration: public VSOObscuration
{
 public:
  VSODiskObscuration(const Eigen::Vector3d& center,
                     const Eigen::Vector3d& normal,
                     double radius):
      VSOObscuration(), fX0(center), fN(normal), fR(radius), fD0()
  {
    fD0 = center.dot(normal);
  }
  virtual ~VSODiskObscuration();
  bool doesObscure(const calin::math::ray::Ray& p_in,
                   calin::math::ray::Ray& p_out, double n) const override;
  VSODiskObscuration* clone() const override;

#ifndef SWIG
  calin::ix::simulation::vs_optics::VSOObscurationData* dump_as_proto(
    calin::ix::simulation::vs_optics::VSOObscurationData* d = nullptr) const override;
#else
  calin::ix::simulation::vs_optics::VSOObscurationData* dump_as_proto() const override;
  void dump_as_proto(calin::ix::simulation::vs_optics::VSOObscurationData* d) const override;
#endif

  static VSODiskObscuration*
  create_from_proto(const ix::simulation::vs_optics::VSODiskObscurationData& d);

  const Eigen::Vector3d& center_pos() const { return fX0; }
  const Eigen::Vector3d& normal() const { return fN; }
  double diameter() const { return 2.0*fR; }

 private:
  Eigen::Vector3d fX0;
  Eigen::Vector3d fN;
  double                  fR;
  double                  fD0;
};

class VSOTubeObscuration: public VSOObscuration
{
 public:
  VSOTubeObscuration(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2,
                     double radius):
      VSOObscuration(), fX1(x1), fX2(x2), fR(radius), fN(), fD1(), fD2()
  {
    fN = x2-x1;
    fN.normalize();
    fD1 = fN.dot(x1);
    fD2 = fN.dot(x2);
    fD  = std::fabs(fD2-fD1);
  }

  virtual ~VSOTubeObscuration();
  bool doesObscure(const calin::math::ray::Ray& p_in,
                   calin::math::ray::Ray& p_out, double n) const override;
  VSOTubeObscuration* clone() const override;

#ifndef SWIG
  calin::ix::simulation::vs_optics::VSOObscurationData* dump_as_proto(
    calin::ix::simulation::vs_optics::VSOObscurationData* d = nullptr) const override;
#else
  calin::ix::simulation::vs_optics::VSOObscurationData* dump_as_proto() const override;
  void dump_as_proto(calin::ix::simulation::vs_optics::VSOObscurationData* d) const override;
#endif

  static VSOTubeObscuration*
  create_from_proto(const ix::simulation::vs_optics::VSOTubeObscurationData& d);

  const Eigen::Vector3d& end1_pos() const { return fX1; }
  const Eigen::Vector3d& end2_pos() const { return fX2; }
  double radius() const { return fR; }
  double diameter() const { return 2.0*fR; }

 private:
  Eigen::Vector3d         fX1;
  Eigen::Vector3d         fX2;
  double                  fR;
  Eigen::Vector3d         fN;
  double                  fD1;
  double                  fD2;
  double                  fD;
};

class VSOAlignedBoxObscuration: public VSOObscuration
{
 public:
  VSOAlignedBoxObscuration(const Eigen::Vector3d& max_corner,
                           const Eigen::Vector3d& min_corner):
      VSOObscuration(), min_corner_(min_corner), max_corner_(max_corner)
  {
    // nothing to see here
  }

  virtual ~VSOAlignedBoxObscuration();
  bool doesObscure(const calin::math::ray::Ray& r_in,
                   calin::math::ray::Ray& r_out, double n) const override;
  VSOAlignedBoxObscuration* clone() const override;

#ifndef SWIG
  calin::ix::simulation::vs_optics::VSOObscurationData* dump_as_proto(
    calin::ix::simulation::vs_optics::VSOObscurationData* d = nullptr) const override;
#else
  calin::ix::simulation::vs_optics::VSOObscurationData* dump_as_proto() const override;
  void dump_as_proto(calin::ix::simulation::vs_optics::VSOObscurationData* d) const override;
#endif

  static VSOAlignedBoxObscuration* create_from_proto(
    const ix::simulation::vs_optics::VSOAlignedBoxObscurationData& d);

  const Eigen::Vector3d& max_corner() const { return max_corner_; }
  const Eigen::Vector3d& min_corner() const { return min_corner_; }

 private:
  Eigen::Vector3d         min_corner_;
  Eigen::Vector3d         max_corner_;
};

class VSOAlignedRectangularAperture: public VSOObscuration
{
public:
  VSOAlignedRectangularAperture(double center_y, double flat_to_flat_xy, bool inverted = false):
    VSOObscuration(), center_(0,0,center_y),
    flat_to_flat_x_2_(0.5*flat_to_flat_xy), flat_to_flat_z_2_(0.5*flat_to_flat_xy),
    inverted_(inverted)
  {
   // nothing to see here
  }

  VSOAlignedRectangularAperture(const Eigen::Vector3d& center,
                                double flat_to_flat_xy, bool inverted = false):
    VSOObscuration(), center_(center),
    flat_to_flat_x_2_(0.5*flat_to_flat_xy), flat_to_flat_z_2_(0.5*flat_to_flat_xy),
    inverted_(inverted)
  {
   // nothing to see here
  }

  VSOAlignedRectangularAperture(const Eigen::Vector3d& center,
                                double flat_to_flat_x, double flat_to_flat_z,
                                bool inverted = false):
    VSOObscuration(), center_(center),
    flat_to_flat_x_2_(0.5*flat_to_flat_x), flat_to_flat_z_2_(0.5*flat_to_flat_z),
    inverted_(inverted)
  {
   // nothing to see here
  }

  virtual ~VSOAlignedRectangularAperture();
  bool doesObscure(const calin::math::ray::Ray& p_in,
                  calin::math::ray::Ray& p_out, double n) const override;
  VSOAlignedRectangularAperture* clone() const override;

#ifndef SWIG
  calin::ix::simulation::vs_optics::VSOObscurationData* dump_as_proto(
    calin::ix::simulation::vs_optics::VSOObscurationData* d = nullptr) const override;
#else
  calin::ix::simulation::vs_optics::VSOObscurationData* dump_as_proto() const override;
  void dump_as_proto(calin::ix::simulation::vs_optics::VSOObscurationData* d) const override;
#endif

  static VSOAlignedRectangularAperture* create_from_proto(
    const ix::simulation::vs_optics::VSOAlignedRectangularApertureData& d);

  const Eigen::Vector3d& center() const { return center_; }
  double flat_to_flat_x_2() const { return flat_to_flat_x_2_; }
  double flat_to_flat_z_2() const { return flat_to_flat_z_2_; }
  bool inverted() const { return inverted_; }

private:
  Eigen::Vector3d         center_;
  double                  flat_to_flat_x_2_;
  double                  flat_to_flat_z_2_;
  bool                    inverted_;
};

class VSOAlignedTileAperture: public VSOObscuration
{
public:
  VSOAlignedTileAperture(const Eigen::Vector3d& center,
      double pitch_x, double pitch_z,
      double center_x, double center_z,
      double support_width_x, double support_width_z):
    VSOObscuration(), center_(center),
    pitch_x_inv_(1.0/pitch_x), pitch_z_inv_(1.0/pitch_z),
    edge_x_(center_x - 0.5*pitch_x),
    edge_z_(center_z - 0.5*pitch_z),
    opaque_frac_x_2_(0.5*support_width_x*pitch_x_inv_),
    opaque_frac_z_2_(0.5*support_width_z*pitch_z_inv_)
  {
   // nothing to see here
  }

  virtual ~VSOAlignedTileAperture();
  bool doesObscure(const calin::math::ray::Ray& p_in,
                  calin::math::ray::Ray& p_out, double n) const override;
  VSOAlignedTileAperture* clone() const override;

#ifndef SWIG
  calin::ix::simulation::vs_optics::VSOObscurationData* dump_as_proto(
    calin::ix::simulation::vs_optics::VSOObscurationData* d = nullptr) const override;
#else
  calin::ix::simulation::vs_optics::VSOObscurationData* dump_as_proto() const override;
  void dump_as_proto(calin::ix::simulation::vs_optics::VSOObscurationData* d) const override;
#endif

  static VSOAlignedTileAperture* create_from_proto(
    const ix::simulation::vs_optics::VSOAlignedTileApertureData& d);

private:
  Eigen::Vector3d         center_;
  double                  pitch_x_inv_;
  double                  pitch_z_inv_;
  double                  edge_x_;
  double                  edge_z_;
  double                  opaque_frac_x_2_;
  double                  opaque_frac_z_2_;
};

class VSOAlignedHexagonalAperture: public VSOObscuration
{
public:
  VSOAlignedHexagonalAperture(double center_y, double flat_to_flat):
    VSOObscuration(), center_(0,0,center_y), flat_to_flat_2_(0.5*flat_to_flat)
  {
   // nothing to see here
  }

  VSOAlignedHexagonalAperture(const Eigen::Vector3d& center, double flat_to_flat):
    VSOObscuration(), center_(center), flat_to_flat_2_(0.5*flat_to_flat)
  {
   // nothing to see here
  }

  virtual ~VSOAlignedHexagonalAperture();
  bool doesObscure(const calin::math::ray::Ray& p_in,
                  calin::math::ray::Ray& p_out, double n) const override;
  VSOAlignedHexagonalAperture* clone() const override;

#ifndef SWIG
  calin::ix::simulation::vs_optics::VSOObscurationData* dump_as_proto(
    calin::ix::simulation::vs_optics::VSOObscurationData* d = nullptr) const override;
#else
  calin::ix::simulation::vs_optics::VSOObscurationData* dump_as_proto() const override;
  void dump_as_proto(calin::ix::simulation::vs_optics::VSOObscurationData* d) const override;
#endif

  static VSOAlignedHexagonalAperture* create_from_proto(
    const ix::simulation::vs_optics::VSOAlignedHexagonalApertureData& d);

private:
  Eigen::Vector3d         center_;
  double                  flat_to_flat_2_;
};

class VSOAlignedCircularAperture: public VSOObscuration
{
public:
  VSOAlignedCircularAperture(double center_y, double diameter, bool inverted = false):
    VSOObscuration(), center_(0,0,center_y), radius_sq_(0.25*diameter*diameter),
    inverted_(inverted)
  {
   // nothing to see here
  }

  VSOAlignedCircularAperture(const Eigen::Vector3d& center, double diameter, bool inverted = false):
    VSOObscuration(), center_(center), radius_sq_(0.25*diameter*diameter),
    inverted_(inverted)
{
   // nothing to see here
  }

  virtual ~VSOAlignedCircularAperture();
  bool doesObscure(const calin::math::ray::Ray& p_in,
                  calin::math::ray::Ray& p_out, double n) const override;
  VSOAlignedCircularAperture* clone() const override;

#ifndef SWIG
  calin::ix::simulation::vs_optics::VSOObscurationData* dump_as_proto(
    calin::ix::simulation::vs_optics::VSOObscurationData* d = nullptr) const override;
#else
  calin::ix::simulation::vs_optics::VSOObscurationData* dump_as_proto() const override;
  void dump_as_proto(calin::ix::simulation::vs_optics::VSOObscurationData* d) const override;
#endif

  static VSOAlignedCircularAperture* create_from_proto(
    const ix::simulation::vs_optics::VSOAlignedCircularApertureData& d);

  const Eigen::Vector3d& center() const { return center_; }
  double radius_sq() const { return radius_sq_; }
  bool inverted() const { return inverted_; }

private:
  Eigen::Vector3d         center_;
  double                  radius_sq_;
  bool                    inverted_;
};

} } } // namespace calin::simulation::vs_optics
