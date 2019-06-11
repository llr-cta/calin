/*

   calin/math/vs_vec3d_obsolete.cpp -- Stephen Fegan -- 2015-11-05

   Class for 3D vector operations. This code is derived from
   simulation code largely implemented by the author at UCLA in
   2004. It is not complient with the calin coding style.

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

/*! \file vec3D.cpp
  vec3D class implementation file

  \author   Stephen Fegan             \n
            UCLA                      \n
            sfegan@astro.ucla.edu     \n
<  \author   Maciej Nicewicz           \n
            UCLA                      \n
            nicewicz@physics.ucla.edu \n
  \author   Vladimir Vassiliev        \n
            UCLA                      \n
            vvv@astro.ucla.edu        \n

  \date     11/09/2004
  \version  1.2
  \note
*/

#include <math/vs_vec3d_obsolete.hpp>

using namespace calin::math::vs_physics::obsolete;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for vector rotation around axis vector.
/// \param axis: axis of rotation vector [rad]
/*! \note
  The modulus of rotation angle is equal to the axis
  norm and is given in radians. The rotation of e_x
  around e_z by PI/2 is equal to e_y.
*/
void Vec3D::Rotate( const Vec3D& axis )
{
  double angle = axis.Norm();

  if(angle < SMALLEST_ROTATION_ANGLE)return;

  Vec3D e_axis = axis / angle;

  double par   = (*this)*e_axis;
  Vec3D p_axis = (*this)-par*e_axis;

#if 0
  double per   = p_axis.Norm();

  if(per == 0.0)
    {
      // vector is parallel to rotation axis
      return ;
    }
  else
    {
      p_axis=p_axis/per;
    }

  *this=par*e_axis+cos(angle)*per*p_axis+sin(angle)*per*(e_axis^p_axis);
#else
  // "Rotation Formula" -- see Goldstein chap 4
  *this = par*e_axis + cos(angle)*p_axis + sin(angle)*(e_axis^p_axis);
#endif

  return;
}

void Vec3D::ScatterDirection(double dispersion, math::rng::RNG& rng)
{
  if(dispersion==0)return;
  double t = Norm();
  if(t==0)return;

  double phi = rng.uniform() * 2.0*M_PI;
  //double theta = atan(dispersion*rng.Normal());
  double theta = dispersion*sqrt(-2*log(rng.uniform()));

  Vec3D tangent_a;

  if((fabs(x)<=fabs(y))&&(fabs(x)<=fabs(z)))tangent_a = (*this)^Vec3D(1,0,0);
  else if(fabs(y)<=fabs(z))tangent_a = (*this)^Vec3D(0,1,0);
  else tangent_a = (*this)^Vec3D(0,0,1);
  tangent_a /= tangent_a.Norm();

  Vec3D tangent_b((*this)^tangent_a);
  tangent_b /= tangent_b.Norm();

#if 1
  *this = (*this)*cos(theta) +
          (tangent_a*cos(phi)+tangent_b*sin(phi))*t*sin(theta);
#else
  Vec3D axis = tangent_a*cos(phi) + tangent_b*sin(phi);
  Rotate(axis*theta);
#endif
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator &= of rotation composition
/*! \note
  Composition is opposite to the sense of matrix
  multiplication. The composition r=r1&r2 is equivalent to
  a rotation first by r1 then by r2.

  For example, for two rotations r1 and r2, and for any vector p

  Vec3D r1(?,?,?);
  Vec3D r2(?,?,?);
  Vec3D p(?,?,?);

  Vec3D r=r1&r2;

  p.Rotate(r1);   \/  is equivalent to p.Rotate(r)
  p.Rotate(r2);   /\
*/

/*
   Composition algorithm comes from consideration of rotation
   as 2-spinors... i.e. write the both rotations in terms of Pauli
   matrices, multiply and then compare terms. See any QM book or
   Goldstein chap 4.

   $\mathbold{R_1} = \mathbold{1} \cos\theta_1 / 2 +
    \hat{n}_1 . \vec{\mathbold{sigma}} \sin\theta_1 / 2 $

   $\mathbold{R_2} = \mathbold{1} \cos\theta_2 / 2 +
    \hat{n}_2 . \vec{\mathbold{sigma}} \sin\theta_2 / 2 $

    Multiply the matrices using the following, collect terms
    and compare with $R_1$ or $R_2$ to deduce the composite
    angle and axis of rotation.

   $(\vec{\mathbold{sigma}}.\vec{A})
    (\vec{\mathbold{sigma}}.\vec{B}) = \mathbold{1}\vec{A}.\vec{B}+
    i\vec{\mathbold{sigma}}(\vec{A}\cross\vec{B})$
*/

Vec3D& Vec3D::operator &= (const Vec3D& r)
{
#define R1 (*this)
#define R2 r

  double r1_theta = R1.Norm();
  double r2_theta = R2.Norm();

  if(r1_theta == 0)
    {
      // R1 is zero so set to R2
      *this = R2;
    }
  else if(r2_theta == 0)
    {
      // R2 is zero so set to R1
      *this = R1;
    }
  else
    {
      double sin_half_r1 = sin(r1_theta/2.0);
      double cos_half_r1 = cos(r1_theta/2.0);
      double sin_half_r2 = sin(r2_theta/2.0);
      double cos_half_r2 = cos(r2_theta/2.0);

      Vec3D r1_hat = R1/r1_theta;
      Vec3D r2_hat = R2/r2_theta;

      double coshalftheta =
	      cos_half_r1*cos_half_r2-
	      sin_half_r1*sin_half_r2*(r1_hat*r2_hat);

      Vec3D axis(r1_hat*sin_half_r1*cos_half_r2+
		    r2_hat*sin_half_r2*cos_half_r1-
		      (r1_hat^r2_hat)*sin_half_r1*sin_half_r2);

      double sinhalftheta = axis.Norm();

      double halftheta=atan(sinhalftheta/coshalftheta);

      *this = axis/sinhalftheta*halftheta*2.0;
    }
  return *this;
}

Vec3D::Vec3D(const ix::common_types::Vector3D& d)
{
  set_from_proto(d);
}

Vec3D::Vec3D(const ix::common_types::Vector3D& d, double scale)
{
  set_from_scaled_proto(d,scale);
}

calin::ix::common_types::Vector3D* Vec3D::
dump_as_proto(ix::common_types::Vector3D* d) const
{
  if(d == nullptr)d = new ix::common_types::Vector3D;
  d->set_x(x);
  d->set_y(y);
  d->set_z(z);
  return d;
}

void Vec3D::set_from_proto(const ix::common_types::Vector3D& d)
{
  x = d.x();
  y = d.y();
  z = d.z();
}

calin::ix::common_types::Vector3D* Vec3D::
dump_scaled_as_proto(double scale, ix::common_types::Vector3D* d) const
{
  if(d == nullptr)d = new ix::common_types::Vector3D;
  d->set_x(x*scale);
  d->set_y(y*scale);
  d->set_z(z*scale);
  return d;
}

void Vec3D::set_from_scaled_proto(const ix::common_types::Vector3D& d, double scale)
{
  x = d.x()*scale;
  y = d.y()*scale;
  z = d.z()*scale;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Print vector components and vector norm in square
void Vec3D::Dump( std::ostream& stream ) const
{
  stream << std::endl;
  stream << " .X:   " << x << std::endl;
  stream << " .Y:   " << y << std::endl;
  stream << " .Z:   " << z << std::endl;
  stream << "Norm2: " << Norm2() << std::endl;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Print vector components in short form
void Vec3D::DumpShort( std::ostream& stream ) const
{
  stream << "( " << x << ' ' << y << ' ' << z << " )";
}

#ifdef TESTMAIN

// Compile with: g++ -O3 -DTESTMAIN -o test Vec3D.cpp -lm

#include<sstream>
#include<vector>

#include"RandomNumbers.hpp"

int main()
{
  RandomNumbers rng("random.seeds");

  Vec3D a(10,20,30);
  Vec3D b(1,2,3);
  std::cout << "sizeof(Vec3D): " << sizeof(Vec3D) << std::endl;
  std::cout << "sizeof(Vec3D*): " << sizeof(Vec3D*) << std::endl << std::endl;

  std::cout << "a:" << a << " b:" << b << std::endl;
  std::cout << "a.Norm2():" << a.Norm2() << "  "
	    << "b.Norm2():" << b.Norm2() << std::endl;
  std::cout << "a.Norm():" << a.Norm() << "  "
	    << "b.Norm():" << b.Norm() << std::endl;

  std::cout << "a+=b: " << (a+=b) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a-=b: " << (a-=b) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl;

  std::cout << "a*=1.5: " << (a*=1.5) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a/=1.5: " << (a/=1.5) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a.Reset(): " << a.Reset() << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl;
  std::cout << "a.Reset(b): " << a.Reset(b) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a.Reset(1,0,0): " << a.Reset(1,0,0) << std::endl;
  std::cout << "b.Reset(0,0,pi/4): " << b.Reset(0,0,M_PI/2) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl;
  std::cout << "a.Rotate(b)" << std::endl;
  a.Rotate(b);
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a.Reset(0,1,0): " << a.Reset(0,1,0) << std::endl;
  std::cout << "b.Reset(0,0,pi/4): " << b.Reset(0,0,M_PI/2) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl;
  std::cout << "a.Rotate(b)" << std::endl;
  a.Rotate(b);
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a.Reset(0,0,1): " << a.Reset(0,0,1) << std::endl;
  std::cout << "b.Reset(0,0,pi/4): " << b.Reset(0,0,M_PI/2) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl;
  std::cout << "a.Rotate(b)" << std::endl;
  a.Rotate(b);
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a.Reset(1,-2,3): " << a.Reset(1,-2,3) << std::endl;
  std::cout << "a.P()" << std::endl;
  a.P();
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a.Reset(1,0,0): " << a.Reset(1,0,0) << std::endl;
  std::cout << "b.Reset(1,0,0): " << b.Reset(1,0,0) << std::endl;
  std::cout << "a^=b:" << (a^=b) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a.Reset(1,0,0): " << a.Reset(1,0,0) << std::endl;
  std::cout << "b.Reset(0,1,0): " << b.Reset(0,1,0) << std::endl;
  std::cout << "a^=b:" << (a^=b) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a.Reset(1,0,0): " << a.Reset(1,0,0) << std::endl;
  std::cout << "b.Reset(0,0,1): " << b.Reset(0,0,1) << std::endl;
  std::cout << "a^=b:" << (a^=b) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  Vec3D c;

  std::cout << "c=b=a=():" << (c=b=a=Vec3D()) << std::endl;
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl << std::endl;

  std::cout << "c=b=(a=())+(1,0,0)):"
	    << (c=b=(a=Vec3D())+Vec3D(1,0,0)) << std::endl;
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl << std::endl;

  std::cout << "a=(1,0,0):" << (a=Vec3D(1,0,0)) << std::endl;
  std::cout << "b=(1,-1,0):" << (b=Vec3D(1,-1,0)) << std::endl;
  std::cout << "c=a+b:" << (c=a+b) << std::endl;
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl << std::endl;

  std::cout << "c=a-b:" << (c=a-b) << std::endl;
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl << std::endl;

  std::cout << "a*b:" << (a*b) << std::endl;
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl << std::endl;

  std::cout << "(3,4,5)*(2,0,0): " << (Vec3D(3,4,5)*Vec3D(2,0,0)) << std::endl;
  std::cout << "(3,4,5)*(0,3,0): " << (Vec3D(3,4,5)*Vec3D(0,3,0)) << std::endl;
  std::cout << "(3,4,5)*(0,0,4): " << (Vec3D(3,4,5)*Vec3D(0,0,4)) << std::endl;
  std::cout << "(3,4,5)*(2,3,0): " << (Vec3D(3,4,5)*Vec3D(2,3,0)) << std::endl;
  std::cout << "(3,4,5)*(2,0,4): " << (Vec3D(3,4,5)*Vec3D(2,0,4)) << std::endl;
  std::cout << "(3,4,5)*(0,3,4): " << (Vec3D(3,4,5)*Vec3D(0,3,4)) << std::endl;
  std::cout << "(3,4,5)*(2,3,4): " << (Vec3D(3,4,5)*Vec3D(2,3,4)) << std::endl;
  std::cout << std::endl;

  std::cout << "(3,4,5)^(2,0,0): " << (Vec3D(3,4,5)^Vec3D(2,0,0)) << std::endl;
  std::cout << "(3,4,5)^(0,3,0): " << (Vec3D(3,4,5)^Vec3D(0,3,0)) << std::endl;
  std::cout << "(3,4,5)^(0,0,4): " << (Vec3D(3,4,5)^Vec3D(0,0,4)) << std::endl;
  std::cout << "(3,4,5)^(2,3,0): " << (Vec3D(3,4,5)^Vec3D(2,3,0)) << std::endl;
  std::cout << "(3,4,5)^(2,0,4): " << (Vec3D(3,4,5)^Vec3D(2,0,4)) << std::endl;
  std::cout << "(3,4,5)^(0,3,4): " << (Vec3D(3,4,5)^Vec3D(0,3,4)) << std::endl;
  std::cout << "(3,4,5)^(2,3,4): " << (Vec3D(3,4,5)^Vec3D(2,3,4)) << std::endl;
  std::cout << std::endl;

  std::cout << "(3,4,5)+(2,0,0): " << (Vec3D(3,4,5)+Vec3D(2,0,0)) << std::endl;
  std::cout << "(3,4,5)+(0,3,0): " << (Vec3D(3,4,5)+Vec3D(0,3,0)) << std::endl;
  std::cout << "(3,4,5)+(0,0,4): " << (Vec3D(3,4,5)+Vec3D(0,0,4)) << std::endl;
  std::cout << "(3,4,5)+(2,3,0): " << (Vec3D(3,4,5)+Vec3D(2,3,0)) << std::endl;
  std::cout << "(3,4,5)+(2,0,4): " << (Vec3D(3,4,5)+Vec3D(2,0,4)) << std::endl;
  std::cout << "(3,4,5)+(0,3,4): " << (Vec3D(3,4,5)+Vec3D(0,3,4)) << std::endl;
  std::cout << "(3,4,5)+(2,3,4): " << (Vec3D(3,4,5)+Vec3D(2,3,4)) << std::endl;
  std::cout << std::endl;

  std::cout << "(3,4,5)-(2,0,0): " << (Vec3D(3,4,5)-Vec3D(2,0,0)) << std::endl;
  std::cout << "(3,4,5)-(0,3,0): " << (Vec3D(3,4,5)-Vec3D(0,3,0)) << std::endl;
  std::cout << "(3,4,5)-(0,0,4): " << (Vec3D(3,4,5)-Vec3D(0,0,4)) << std::endl;
  std::cout << "(3,4,5)-(2,3,0): " << (Vec3D(3,4,5)-Vec3D(2,3,0)) << std::endl;
  std::cout << "(3,4,5)-(2,0,4): " << (Vec3D(3,4,5)-Vec3D(2,0,4)) << std::endl;
  std::cout << "(3,4,5)-(0,3,4): " << (Vec3D(3,4,5)-Vec3D(0,3,4)) << std::endl;
  std::cout << "(3,4,5)-(2,3,4): " << (Vec3D(3,4,5)-Vec3D(2,3,4)) << std::endl;
  std::cout << std::endl;

  std::cout << "-(3,4,5): " << (-Vec3D(3,4,5)) << std::endl;
  std::cout << "1.5*(3,4,5): " << (1.5*Vec3D(3,4,5)) << std::endl;
  std::cout << "(3,4,5)*1.5: " << (Vec3D(3,4,5)*1.5) << std::endl;
  std::cout << "(3,4,5)/1.5: " << (Vec3D(3,4,5)/1.5) << std::endl << std::endl;

  std::cout << "std::istringstream s(\"(10,11,12)\"); s >> a;" << std::endl;
  std::istringstream s("(10,11,12)");
  s >> a;
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl << std::endl;

  for(unsigned nrot = 1; nrot<=20; nrot++)
    {
      double mean_residual = 0;
      double max_residual = 0;

      for(unsigned n=0; n<1000; n++)
	{
	  double phi = rng.Uniform() * M_PI * 2.0;
	  double costheta = rng.Uniform() * 2.0 - 1.0;
	  double sintheta = sin(acos(costheta));
	  Vec3D a(sintheta*cos(phi),sintheta*sin(phi),costheta);
	  Vec3D b(a);
	  Vec3D rc;

	  std::vector<Vec3D> R;
	  std::vector<Vec3D> Rc;

	  for(unsigned i=0; i<nrot; i++)
	    {
	      phi = rng.Uniform() * M_PI * 2.0;
	      costheta = rng.Uniform() * 2.0 - 1.0;
	      sintheta = sin(acos(costheta));
	      Vec3D r(sintheta*cos(phi),sintheta*sin(phi),costheta);
	      r *= rng.Uniform()*2.0*M_PI;

	      b.Rotate(r);
	      rc &= r;

	      R.push_back(r);
	      Rc.push_back(rc);
	    }

	  a.Rotate(rc);

	  Vec3D c(a-b);
	  double residual = c.Norm();
	  if(residual>max_residual)max_residual=residual;
	  mean_residual+=residual;

#if 0
	  if(residual>1e-14)
	    {
	      for(unsigned i=0; i<R.size(); i++)
		std::cout << R[i].Norm()
			  << ' '  << R[i] << ' ' << Rc[i] << std::endl;
	      std::cout << a << ' ' << b << ' ' << c << std::endl;
	    }
#endif
	}

      std::cout << nrot << " rotations: mean residual: "
		<< mean_residual/1000.0
		<< " max residual: " << max_residual << std::endl;
    }

  Vec3D r1=M_PI/2*Vec3D(1,0,0);
  Vec3D r2=M_PI/2*Vec3D(0,1,0);
  Vec3D rc = r1 & r2;
  std::cout << "r1:" << r1 << " r2:" << r2 << " rc:" << rc << std::endl;

  a=Vec3D(1,0,0); b=a; b.Rotate(r1); b.Rotate(r2); c=a; c.Rotate(rc);
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl;

  a=Vec3D(0,1,0); b=a; b.Rotate(r1); b.Rotate(r2); c=a; c.Rotate(rc);
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl;

  a=Vec3D(0,0,2); b=a; b.Rotate(r1); b.Rotate(r2); c=a; c.Rotate(rc);
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl;

  a=Vec3D(1,1,0); b=a; b.Rotate(r1); b.Rotate(r2); c=a; c.Rotate(rc);
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl << std::endl;

  a=Vec3D(.1,0,0);
  a &= a & a & a;
  std::cout << "a:" << a << std::endl;
}

#endif // TESTMAIN
