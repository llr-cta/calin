/*

   calin/simulation/vso_obscuration.hpp -- Stephen Fegan -- 2015-11-09

   Class for telescope data and functionality. This code is derived
   from raytracing code largely implemented by the author at UCLA in
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
                     double radius, bool incoming_only):
      VSOObscuration(), fX0(center), fN(normal), fR(radius), fD0(),
      fICO(incoming_only)
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
  bool incoming_only() const { return fICO; }

 private:
  Eigen::Vector3d fX0;
  Eigen::Vector3d fN;
  double                  fR;
  double                  fD0;
  bool                    fICO;
};

class VSOTubeObscuration: public VSOObscuration
{
 public:
  VSOTubeObscuration(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2,
                     double radius, bool incoming_only):
      VSOObscuration(), fX1(x1), fX2(x2), fR(radius), fN(), fD1(), fD2(),
      fICO(incoming_only)
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
  double diameter() const { return 2.0*fR; }
  bool incoming_only() const { return fICO; }

 private:
  Eigen::Vector3d         fX1;
  Eigen::Vector3d         fX2;
  double                  fR;
  Eigen::Vector3d         fN;
  double                  fD1;
  double                  fD2;
  double                  fD;
  bool                    fICO;
};

class VSOAlignedBoxObscuration: public VSOObscuration
{
 public:
  VSOAlignedBoxObscuration(const Eigen::Vector3d& max_corner,
                           const Eigen::Vector3d& min_corner,
                           bool incoming_only):
      VSOObscuration(), min_corner_(min_corner), max_corner_(max_corner),
      incoming_only_(incoming_only)
  {
    // nothing to see here
  }

  virtual ~VSOAlignedBoxObscuration();
  bool doesObscure(const calin::math::ray::Ray& p_in,
                   calin::math::ray::Ray& p_out, double n) const override;
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
  bool incoming_only() const { return incoming_only_; }

 private:
  Eigen::Vector3d         min_corner_;
  Eigen::Vector3d         max_corner_;
  bool                    incoming_only_;
};

} } } // namespace calin::simulation::vs_optics
