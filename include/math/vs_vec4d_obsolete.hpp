/* 

   calin/math/vs_vec4d.hpp -- Stephen Fegan -- 2015-11-05

   Class for 4-vector operations. This code is derived from simulation
   code largely implemented by Vladimir Vassiliev at UCLA in 2004. It
   is not complient with the calin coding style.

   Some portions are :

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

/*! \file Vec4D.hpp
  Vec4D class header file
  
  \author   Vladimir Vassiliev        \n
            UCLA                      \n
	    vvv@astro.ucla.edu        \n

  \date     07/30/2004
  \version  1.0
  \note
*/

#pragma once

#include <cmath>

#include <math/vs_vec3d.hpp>

namespace calin { namespace math { namespace vs_physics {

/*!  \class Vec4D
  \brief 4 dimensional Lorentz vector class
  
  This class defines 4D vectors, set of operations 
  in vector field, scalar product, rotations and boosts 
  as well as parity and time transformations.
     
*/
  
class Vec4D
{
 public:
  inline Vec4D();                  //!<default constructor
  inline Vec4D( const Vec4D& v );  //!<copy constructor
  inline Vec4D( double _r0, Vec3D _r );   //!<constructor
  inline Vec4D( double _r0, double r1, double r2, double r3 );  //!<constructor
    
  inline double Norm2() const;     //!<calculates scalar product with itself

  inline void Rotate(const Vec3D& axis);    //!<rotates vector around axis
  inline bool Boost (const Vec3D& axis);    //!<boosts vector
  inline Vec3D SRFBoost() const;   //!<finds boost to special reference frame
  inline void P();                 //!<parity transformation
  inline void T();                 //!<time inversion
  inline void Reflect(const Vec3D& norm); //!< reflect vector

  inline Vec4D& Reset(const Vec4D& v = Vec4D());
  inline Vec4D& Reset(double _r0, double r1, double r2, double r3);

  inline Vec4D& operator = ( const Vec4D& v );  //!<assignment
  inline Vec4D& operator += ( const Vec4D& v ); //!<assignment: addition
  inline Vec4D& operator -= ( const Vec4D& v ); //!<assignment: subtraction 
  inline Vec4D& operator *= ( double d );       //!<assignment: multiply by scaler
  inline Vec4D& operator /= ( double d );       //!<assignment: divide by scaler

  inline Vec4D  operator + ( const Vec4D& v ) const;  //!<addition
  inline Vec4D  operator - ( const Vec4D& v ) const;  //!<subtraction
  inline double operator * ( const Vec4D& v ) const;  //!<scalar product

  void Dump(std::ostream& stream = std::cout) const; //!<prints coordinates
  void DumpShort(std::ostream& stream = std::cout) const; //!<prints coordinates

 public:
  double r0;   //!< time component
  Vec3D r;     //!< space components
};
  
inline Vec4D operator - ( const Vec4D& v );           //!<negation
inline Vec4D operator * ( double d, const Vec4D& v ); //!<left scalar mult.
inline Vec4D operator * ( const Vec4D& v, double d ); //!<right scalar mult.
inline Vec4D operator / ( const Vec4D& v, double d ); //!<right scalar division

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Default class constructor
inline Vec4D::Vec4D(): r0(), r()
{
  // nothing to see here
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Copy constructor
inline Vec4D::Vec4D(const Vec4D& v): r0(v.r0), r(v.r)
{
  // nothing to see here
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Overloaded class constructor
/// \param _r0: time component
/// \param _r: space components
inline Vec4D::Vec4D( double _r0, Vec3D _r ): r0(_r0), r(_r)
{
  // nothing to see here
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Overloaded class constructor
/// \param _r0: time component
/// \param r1: space components in x-direction
/// \param r2: space components in y-direction
/// \param r3: space components in z-direction
inline Vec4D::Vec4D( double _r0, double r1, double r2, double r3 ): 
    r0(_r0), r(r1,r2,r3)
{
  // nothing to see here
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for vector norm in square
inline double Vec4D::Norm2() const
{
  return r0*r0-r*r;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for vector rotation around axis vector.
/// \param axis: axis of rotation vector [rad] 
/*! \note 
  The modulus of rotation angle is equal to the axis 
  norm and is given in radians. The rotation of e_x 
  around e_z by PI/2 is equal to e_y.
*/
inline void Vec4D::Rotate(const Vec3D& axis)
{
  r.Rotate(axis);
  return;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for vector boost
/// \param boost: boost vector (boost,boost) < 1.0 
/*! \note 
  
 */
inline bool Vec4D::Boost( const Vec3D& boost )
{
  double s=boost.Norm();
    
  if( s >= (double)1. ){
    return false;     // boost is invalid
  }

  Vec3D e_boost;
  double gamma=1./sqrt(1.-(boost*boost));
  if( s != (double)0 ){
    e_boost=boost/s;
  } 
  else {
    return true;      // no boost required, boost=0..
  }


  double par=(r*e_boost);
  Vec4D t(gamma*(r0-(boost*r)),
          gamma*(par*e_boost-r0*boost)+r-par*e_boost);
  
  *this=t; 

  if(gamma >=1.0E+6) {
    std::cout <<"Boost: Boost accuracy may be lost " << std::endl;
    std::cout <<"       Boost is too large > 1E+6 " << std::endl;
  }
    
  return true;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for finding boost to special r.f.
inline Vec3D Vec4D::SRFBoost() const
{
  double n=Norm2();
  if( n > 0.) {            // time-like vector 
    // e.g. physical particle
    
    return r/r0;           // boost to r.f. where r=0 
    // i.e. particle's rest frame
      
  } else if ( n < 0.) {    // space-like vector 
    // e.g. tachion
      
    return r0*r/(r*r);     // boost to r.f. where r0=0
      
  }
  // light-like
  return Vec3D();          // boost=0. is returned
}
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for parity transformation of Lorentz vector
inline void Vec4D::P()
{
  r=-r;  
  return;
}
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for time inversion of Lorentz vector
inline void Vec4D::T( )
{
  r0=-r0;
  return;
}
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for reflection in mirror (defined by norm)
inline void Vec4D::Reflect(const Vec3D& norm)
{
  r.Reflect(norm);
}
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Reset the vector (assignment)
/// \param v: vector
inline Vec4D& Vec4D::Reset(const Vec4D& v)
{
  return *this = v;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Reset the vector components
/// \param _r0: time component
/// \param r1: space components in x-direction
/// \param r2: space components in y-direction
/// \param r3: space components in z-direction
inline Vec4D& Vec4D::Reset(double _r0, double r1, double r2, double r3)
{
  r0 = _r0;
  r.x = r1;
  r.y = r2;
  r.z = r3;
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator =
inline Vec4D& Vec4D::operator = ( const Vec4D& v )
{
  r0 = v.r0;
  r = v.r;
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Assignment operator +=
inline Vec4D& Vec4D::operator += ( const Vec4D& v )
{
  r0 += v.r0;
  r += v.r;
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Assignment operator -=
inline Vec4D& Vec4D::operator -= ( const Vec4D& v )
{
  r0 -= v.r0;
  r -= v.r;
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Assignment operator *= scalar multiplication
/// \param d: scalar
inline Vec4D& Vec4D::operator *= ( double d )
{
  r0 *= d;
  r *= d;
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Assignment operator *= scalar division
/// \param d: scalar
inline Vec4D& Vec4D::operator /= ( double d )
{
  r0 /= d;
  r /= d;
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator +
inline Vec4D Vec4D::operator + ( const Vec4D& v ) const
{
  Vec4D temp(*this);
  return temp+=v;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator -
inline Vec4D Vec4D::operator - ( const Vec4D& v ) const
{
  Vec4D temp(*this);
  return temp-=v;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator *
inline double Vec4D::operator * ( const Vec4D& v ) const
{
  return r0*v.r0-r*v.r;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator * for left scalar multiplication
inline Vec4D operator * ( double a, const Vec4D& v )
{    
  return Vec4D(a*v.r0, a*v.r);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator * for right scalar multiplication
inline Vec4D operator * ( const Vec4D& v, double a )
{
  return Vec4D(v.r0*a, v.r*a);
}
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator * for right scalar division
inline Vec4D operator / ( const Vec4D& v, double a )
{
  return Vec4D( v.r0/a, v.r/a );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Negation
inline Vec4D operator - ( const Vec4D& v )
{
  return Vec4D( -v.r0, -v.r );
}
 
} } } // namespace calin::math::vs_physics

#ifndef SWIG
namespace std {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Stream insertion
inline ostream& operator<< (ostream& stream,
                            const calin::math::vs_physics::Vec4D& v)
{
  v.DumpShort(stream);
  return stream;
}
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Stream extraction
inline istream& operator>> (istream& stream, calin::math::vs_physics::Vec4D& v)
{
  char c;
  stream >> c >> v.r0 >> c >> v.r.x >> c >> v.r.y >> c >> v.r.z >> c;
  return stream;
}

} // namespace std
#endif
