/* 

   calin/math/vs_vec3d.hpp -- Stephen Fegan -- 2015-11-05

   Class for 3-vector operations. This code is derived from
   simulation code largely implemented by the author at UCLA in
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

/*! \file Vec3D.hpp
  Vec3D class header file
  
  \author   Stephen Fegan             \n
            UCLA                      \n
	    sfegan@astro.ucla.edu     \n
  \author   Maciej Nicewicz           \n
            UCLA                      \n
	    nicewicz@physics.ucla.edu \n
  \author   Vladimir Vassiliev        \n
            UCLA                      \n
	    vvv@astro.ucla.edu        \n

  \date     11/09/2004
  \version  1.2
  \note
*/

#pragma once

#include <iostream>
#include <cmath>
#include <proto/calin_common_types.pb.h>

/*!  \class Vec3D
     \brief 3 dimensional vector class

     This class defines 3D vectors, set of operations 
     in vector field, scalar and vector products, 
     as well as rotation, parity transformation, negation,
     and normalization.

*/

namespace calin { namespace math {

namespace random_numbers {
class RandomNumbers;
}

namespace vs_physics {

class Vec3D
{
 public:
  inline Vec3D();                                   //!<default constructor
  inline Vec3D(const Vec3D& v);                     //!<copy constructor
  inline Vec3D( double _x, double _y, double _z );  //!<overloaded constructor
  Vec3D(const ix::Vector3D& d); // construct from protobuf
  
  inline void Polar(double& r, double& theta, double& phi) const;

  inline double Norm() const;              //!<calculates norm
  inline double Norm2() const;             //!<calculates scalar product with itself
  void Rotate(const Vec3D& axis);          //!<rotates vector around axis
  inline void P();                         //!<parity transformation
  inline void Reflect(const Vec3D& norm);  //!reflect in normal

  void ScatterDirection(double dispersion, random_numbers::RandomNumbers& rng);

  inline Vec3D& Reset(const Vec3D& v = Vec3D());
  inline Vec3D& Reset(double _x, double _y, double _z);

  inline bool operator == ( const Vec3D& v ) const;
  inline bool operator != ( const Vec3D& v ) const;
  
  inline Vec3D& operator = ( const Vec3D& v );  //!<assignment
  inline Vec3D& operator += ( const Vec3D& v ); //!<assignment: addition
  inline Vec3D& operator -= ( const Vec3D& v ); //!<assignment: subtraction 
  inline Vec3D& operator ^= ( const Vec3D& v ); //!<vector product  
  inline Vec3D& operator *= ( double d );       //!<assignment: multiply by scaler
  inline Vec3D& operator /= ( double d );       //!<assignment: divide by scaler

  Vec3D& operator &= (const Vec3D& r);   //!<assignment: composition of rotations

  inline Vec3D  operator + (const Vec3D& v) const;  //!<addition
  inline Vec3D  operator - (const Vec3D& v) const;  //!<subtraction
  inline double operator * (const Vec3D& v) const;  //!<scalar product
  inline Vec3D  operator ^ (const Vec3D& v) const;  //!<vector product  

  inline Vec3D  operator & (const Vec3D& v) const;  //!<addition of rotations

  void dump_to_proto(ix::Vector3D* d) const;
  void set_from_proto(const ix::Vector3D& d);
  
  void Dump(std::ostream& stream = std::cout) const; //!<prints coordinates
  void DumpShort(std::ostream& stream = std::cout) const; //!<prints coordinates

  inline Vec3D DeRotate(const Vec3D& r) const;
    
 public:
  double x, y, z;   //!<components

 private:
  constexpr static double SMALLEST_ROTATION_ANGLE = 1e-9;
};

inline Vec3D operator - ( const Vec3D& v );           //!<negation
inline Vec3D operator * ( double d, const Vec3D& v ); //!<left scalar mult.
inline Vec3D operator * ( const Vec3D& v, double d ); //!<right scalar mult.
inline Vec3D operator / ( const Vec3D& v, double d ); //!<right scalar division

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Default class constructor
inline Vec3D::Vec3D(): x(), y(), z()
{
  // nothing to see here
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Copy constructor
// \param o: Vector to copy
inline Vec3D::Vec3D(const Vec3D& o): x(o.x), y(o.y), z(o.z)
{
  // nothing to see here
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Overloaded class constructor
/// \param _x: x
/// \param _y: y
/// \param _z: z
inline Vec3D::Vec3D( double _x, double _y, double _z ): x(_x), y(_y), z(_z)
{
  // nothing to see here
}
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for vector norm
/// \param _r: r
/// \param _theta: theta
/// \param _phi: phi
inline void Vec3D::Polar(double& r, double& theta, double& phi) const
{
  r=Norm();
  theta = atan(sqrt(x*x+y*y)/z);
  phi = atan2(y,x);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for vector norm
inline double Vec3D::Norm() const
{
  return sqrt(x*x + y*y + z*z);
}
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for vector norm in square
inline double Vec3D::Norm2() const
{
  return x*x + y*y + z*z;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for parity transformation of vector
inline void Vec3D::P()
{
  *this=-(*this);
  return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for reflection of vector in mirror
inline void Vec3D::Reflect(const Vec3D& norm)  //!reflect in normal
{
  double n2 = norm.Norm2();
  if(n2 == 0)return;
  *this -= norm * (2.0*((*this)*norm)/n2);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Reset the vector (assignment)
/// \param v: vector
inline Vec3D& Vec3D::Reset(const Vec3D& v)
{
  return *this = v;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Reset the vector components
/// \param _x: x
/// \param _y: y
/// \param _z: z
inline Vec3D& Vec3D::Reset(double _x, double _y, double _z)
{
  x = _x;
  y = _y;
  z = _z;
  return *this;
}

inline bool Vec3D::operator == ( const Vec3D& v ) const
{
  return x==v.x && y==v.y && z==v.z;
}

inline bool Vec3D::operator != ( const Vec3D& v ) const
{
  return x!=v.x || y!=v.y || z!=v.z;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Assignment operator =
inline Vec3D& Vec3D::operator = ( const Vec3D& v )
{
  x = v.x;
  y = v.y;
  z = v.z;
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Assignment operator +=
inline Vec3D& Vec3D::operator += ( const Vec3D& v )
{
  x += v.x;
  y += v.y;
  z += v.z;
  return *this;
}
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Assignment operator -=
inline Vec3D& Vec3D::operator -= ( const Vec3D& v )
{
  x -= v.x;
  y -= v.y;
  z -= v.z;
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Assignment operator ^= vector multiplication
/// \param v: vector
/*! \note 
  Normally ^ operator is a bitwise exclusive OR, which 
  has precedence lower than * / and even lower than + -.
  Thus it is executed last if no brackets are used.
  Exmp: r=7.5*a+l*b-c*m+2*b^c=(7.5*a+l*b-c*m+2*b)^c
  See examples file exmp_Vec3D.cpp for additional
  comment(s)
*/
inline Vec3D& Vec3D::operator ^= ( const Vec3D& v )
{
  Vec3D temp(y*v.z - z*v.y,
             z*v.x - x*v.z,
             x*v.y - y*v.x);
  *this=temp;
  return *this;
}
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Assignment operator *= scalar multiplication
/// \param d: scalar
inline Vec3D& Vec3D::operator *= ( double d )
{
  x *= d;
  y *= d;
  z *= d;
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Assignment operator /= scalar division
/// \param d: scalar
inline Vec3D& Vec3D::operator /= ( double d )
{
  x /= d;
  y /= d;
  z /= d;
  return *this;
}
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator +
inline Vec3D Vec3D::operator + ( const Vec3D& v ) const
{
  Vec3D temp(*this);
  return temp += v;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator -
inline Vec3D Vec3D::operator - ( const Vec3D& v ) const
{
  Vec3D temp(*this);
  return temp -= v;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator *
inline double Vec3D::operator * ( const Vec3D& v ) const
{
  return x*v.x + y*v.y + z*v.z;
}
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator ^ of vector product
/*! \note 
  Normally ^ operator is a bitwise exclusive OR, which 
  has precedence lower than * / and even lower than + -.
  Thus it is executed last if no brackets are used.
  Exmp: r=7.5*a+l*b-c*m+2*b^c=(7.5*a+l*b-c*m+2*b)^c
  See examples file exmp_Vec3D.cpp for additional
  comment(s)
*/
inline Vec3D Vec3D::operator ^ ( const Vec3D& v ) const
{
  Vec3D temp(*this);
  return temp ^= v;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator & of rotation composition
/*! \note 
  Composition is opposite to the sense of matrix
  multiplication. If r=r1&r2 then rotation r1 happens
  first, followed by r2
  
  Normally & operator is a bitwise exclusive AND, which 
  has precedence lower than * / and even lower than + -.
  Thus it is executed last if no brackets are used.
  Exmp: r=7.5*a+l*b-c*m+2*b^c=(7.5*a+l*b-c*m+2*b)^c
  See examples file exmp_Vec3D.cpp for additional
  comment(s)
*/
inline Vec3D Vec3D::operator & ( const Vec3D& v ) const
{
  Vec3D temp(*this);
  return temp &= v;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// DeRotate a vector to a system which is parallel to another
/// vector
inline Vec3D Vec3D::DeRotate(const Vec3D& r) const
{
  Vec3D temp(*this); 
  double norm;
  double theta;
  double phi;
  r.Polar(norm,theta,phi);
  temp.Rotate(Vec3D(0,0,-phi));
  temp.Rotate(Vec3D(0,-theta,0));
  return temp;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator * for left scalar multiplication
inline Vec3D operator * ( double d, const Vec3D& v )
{
  Vec3D temp(d*v.x, d*v.y, d*v.z );
  return temp;
}
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator * for right scalar multiplication
inline Vec3D operator * ( const Vec3D& v, double d )
{
  Vec3D temp( v.x*d, v.y*d, v.z*d );
  return temp;
}
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator * for right scalar division
inline Vec3D operator / ( const Vec3D& v, double d )
{
  Vec3D temp( v.x/d, v.y/d, v.z/d );
  return temp;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Negation
inline Vec3D operator - ( const Vec3D& v )
{
  Vec3D temp( -v.x, -v.y, -v.z );
  return temp;
}

} } } // namespace calin::math::vs_physics

#ifndef SWIG
namespace std {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Stream insertion
inline ostream& operator << (ostream& stream,
                             const calin::math::vs_physics::Vec3D& v)
{
  v.DumpShort(stream);
  return stream;
}
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Stream extraction
inline istream& operator >> (istream& stream,
                             calin::math::vs_physics::Vec3D& v)
{
  char c;
  stream >> c >> v.x >> c >> v.y >> c >> v.z >> c;
  return stream;
}

} // namespace std
#endif
