/*

   calin/math/vs_vec3d.hpp -- Stephen Fegan -- 2016-10-22

   Class for 3-vector operations. This code is a shallow layer around the
   Eigen::Vector3d class proving some level of compatibility with the Vec3D
   class from simulation code largely implemented by the author at UCLA in
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

#pragma once

#include <iostream>
#include <cmath>
#include <calin_global_definitions.hpp>
#include <common_types.pb.h>
#include <Eigen/Dense>
#include <math/rng.hpp>

namespace calin { namespace math { namespace vs_physics {

class Vec3D: public Eigen::Vector3d
{
public:
  Vec3D(): Eigen::Vector3d(0,0,0) { /* nothing to see here */ }
  Vec3D(const Eigen::Vector3d& v):
    Eigen::Vector3d(v) { /* nothing to see here */ }
  Vec3D(double x, double y, double z):
    Eigen::Vector3d(x,y,z) { /* nothing to see here */ }
  Vec3D(const Eigen::Vector3d& v, double scale):
    Eigen::Vector3d(v*scale) { /* nothing to see here */ }

  using Eigen::Vector3d::operator=;
  using Eigen::Vector3d::operator+=;
  using Eigen::Vector3d::operator-=;
  using Eigen::Vector3d::operator*=;
  using Eigen::Vector3d::operator/=;
  using Eigen::Vector3d::operator+;
  using Eigen::Vector3d::operator-;
  using Eigen::Vector3d::operator*;
  using Eigen::Vector3d::operator/;
  using Eigen::Vector3d::operator==;
  using Eigen::Vector3d::operator!=;

  //inline void Polar(double& r, double& theta, double& phi) const;

  double Norm() const { return norm(); }
  double Norm2() const { return squaredNorm(); }
  //void Rotate(const Vec3D& axis);          //!<rotates vector around axis
  void P() { *this = -(*this); }
  void Reflect(const Vec3D& mirror_norm) {
    double n2 = mirror_norm.Norm2();
    if(n2 == 0)return;
    *this -= mirror_norm * (2.0*((*this).dot(mirror_norm))/n2);
  }

  void ScatterDirection(double dispersion, math::rng::RNG& rng);

  Vec3D& Reset(const Vec3D& v = Vec3D()) {
    *this = v; return *this; }
  Vec3D& Reset(double x, double y, double z) {
    *this = Eigen::Vector3d(x,y,z); return *this; }

  //inline bool operator == ( const Vec3D& v ) const;
  //inline bool operator != ( const Vec3D& v ) const;

  //inline Vec3D& operator = ( const Vec3D& v );  //!<assignment
  //inline Vec3D& operator += ( const Vec3D& v ); //!<assignment: addition
  //inline Vec3D& operator -= ( const Vec3D& v ); //!<assignment: subtraction
  //inline Vec3D& operator *= ( double d );       //!<assignment: multiply by scaler
  //inline Vec3D& operator /= ( double d );       //!<assignment: divide by scaler

  //Vec3D& operator &= (const Vec3D& r);   //!<assignment: composition of rotations

  //inline Vec3D  operator + (const Vec3D& v) const;  //!<addition
  //inline Vec3D  operator - (const Vec3D& v) const;  //!<subtraction
  double operator * (const Vec3D& v) const {
    return this->dot(v); }
  Vec3D& operator ^= ( const Vec3D& v ) {
    *this = this->cross(v); return *this; }
  Vec3D operator ^ (const Vec3D& v) const {
    return this->cross(v); }

  //inline Vec3D  operator & (const Vec3D& v) const;  //!<addition of rotations

  //inline operator Eigen::Vector3d() const;

  ix::common_types::Vector3D* dump_as_proto(
    ix::common_types::Vector3D* d = nullptr) const;
  void set_from_proto(const ix::common_types::Vector3D& d);

  ix::common_types::Vector3D* dump_scaled_as_proto(
    double scale, ix::common_types::Vector3D* d = nullptr) const;
  void set_from_scaled_proto(const ix::common_types::Vector3D& d, double scale);

  Vec3D(const ix::common_types::Vector3D& d) {
    set_from_proto(d); }
  Vec3D(const ix::common_types::Vector3D& d, double scale) {
    set_from_scaled_proto(d, scale); }

  //void Dump(std::ostream& stream = std::cout) const; //!<prints coordinates
  //void DumpShort(std::ostream& stream = std::cout) const; //!<prints coordinates

  //inline Vec3D DeRotate(const Vec3D& r) const;
};

} } } // namespace calin::math::vs_physics
