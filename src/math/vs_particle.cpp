/*

   calin/math/vs_particle.cpp -- Stephen Fegan -- 2015-11-05

   Class for 3D vector operations. This code is derived from
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

//-*-mode:c++; mode:font-lock;-*-

/*! \file particle.cpp
  relativistic particle class file

  \author   Stephen Fegan             \n
            UCLA                      \n
            sfegan@astro.ucla.edu     \n

  \author   Vladimir Vassiliev        \n
            UCLA                      \n
	    vvv@astro.ucla.edu        \n

  \author   Maciej Nicewicz           \n
            UCLA                      \n
	    nicewicz@physics.ucla.edu \n


  \date     11/09/2004
  \version  1.1
  \note
*/

#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>
#include <cassert>

#include <math/vs_particle.hpp>

using namespace calin::math::vs_physics;

namespace {

std::string DoubleToString( double d )
{
  std::ostringstream ss;
  ss << std::scientific << std::setprecision( 25 ) << d;
  return ss.str();
}

} // anonymous namespace

bool                                                      Particle::s_pdgdata_init = false;
std::map<int, PDGData*>                         Particle::s_pdgdata;
std::map<std::pair<std::string, int>, PDGData*> Particle::s_pdgdata_name_q;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Print vector components and vector norm in square
void Particle::Dump(std::ostream& stream)
{
  stream << std::endl<<"Position:" << std::endl;
  stream << "r.T:   " << r.r0 << std::endl;
  stream << "r.X:   " << r.r.x << std::endl;
  stream << "r.Y:   " << r.r.y << std::endl;
  stream << "r.Z:   " << r.r.z << std::endl;
  stream << "(r,r): " << (r*r)<<std::endl;

  stream<<std::endl<<"Momenta:"<<std::endl;
  stream << "p.T:   " << p.r0 << std::endl;
  stream << "p.X:   " << p.r.x << std::endl;
  stream << "p.Y:   " << p.r.y << std::endl;
  stream << "p.Z:   " << p.r.z << std::endl;
  stream << "(p,p): " << (p*p)<<std::endl;

  stream<<std::endl<<"Charge:"<<charge<<std::endl;
  stream<<std::endl<<"mc2:  "<<sqrt(m2c4)<<std::endl;

  return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*! Propagate particle by positive dist in magnetic field
    \param dist: propagation dist >= 0.      [cm]
    \param    B: magnetic field vector    [gauss]

    \note
     Returns 'true' if particle is propagated, and returns 'false'
     if particle is at rest

     Equation of motion in International System of Units (SI)
     d/dt(p)=q e [v x B]

     q - is particle charge measured in electron charge units
     e - absolute value of electron charge e=1.6021892e-19 Coulomb

     Because magnetic field is involved, the conversion of units
     is commented below in detail for the record

     d/dt(pc)=q e c^2 [v/c x B]
     d/dt(pc)=q e c^2 (1 Tesla) [v/c x B/(1 Tesla)]
     d/dt(pc)=q e c^2 (1 Volt sec m^-2) [v/c x B/(10^4 Gauss)]
     d/dt(pc)=q e c^2 (1 Volt sec cm^-2) 1e-8 [v/c x B/(1 Gauss)]
     d/dt(pc)=(1/sec) q e(1 Volt) (c/1 cm/sec)^2 1e-8 [v/c x B/(1 Gauss)]
     d/dt(pc)=(1/sec) q E0 [v/c x B/(1 Gauss)]
     where E0= e (1 Volt) (c/1 cm/sec)^2 1e-8, and c=2.997924576e+10 cm/s

     e (1 Volt) = 1.6021892e-19 Joule = 1.6021892e-12 erg
     E0 = 1.6021892e-12 erg (2.997924576e+10)^2 1e-8=
        = 1.6021892 (2.997924576)^2 erg = 14.399758 erg

     d/dt(E v/c)=(1/sec) q E0 [v/c x B/(1 Gauss)]
     d/dt(v/c)=(1/sec) q E0/E [v/c x B/(1 Gauss)]

     ---------
     SJF notes
     ---------

     1) The Physics::Constants class has these constants so
        E0 = Constants::cgs_eV * Constants::cgs_c * Constants::cgs_c * 1e-8

     2) For reference the equation of motion in CGS is

        with [e] = e.s.u.
             [B] = gauss
	     [E] = erg

        d/dt(p)                 = q e [v/c x B]
        d/dt(gammma m0 v)       = q e [v/c x B]
        d/dt(gammma m0 c^2 v/c) = q e c [v/c x B]
        d/dt(E v/c)             = q e c [v/c x B]
        d/dt(v/c)               = q e c 1/E [v/c x B]
*/

#define B_EQUATION_CGS

bool Particle::PropagateBField(double dist, const Vec3D& B)
{
  double b=B.Norm();

  // wDc -- frequency of gyration divided by c
#ifdef B_EQUATION_VVV
  double wDc = 1./ 2.997924576e+10   // 1/sec / c
             * charge                // q
             * (14.399758/p.r0)      // E0/(1 erg) * (1 erg)/E
             * b;                    // |B|/ (1 gauss)
#endif

#ifdef B_EQUATION_CONSTANTS
  constexpr double E0 = calin::math::constants::cgs_eV *
                        calin::math::constants::cgs_c *
                        calin::math::constants::cgs_c * 1e-8;
  double wDc = 1./ calin::math::constants::cgs_c  // 1/sec / c
             * charge                             // q
             * (E0/p.r0)                          // E0/(1 erg) * (1 erg)/E
             * b;                                 // |B|/(1 gauss)
#endif

#ifdef B_EQUATION_CGS
  double wDc = charge                             // q     [electic charge]
               * calin::math::constants::cgs_e    // e     [e.s.u.]
             / p.r0                               // 1/E   [1/erg]
             * b;                                 // |B|   [gauss]
#endif

  if ( wDc == 0 )                    // charge=0 or B=0 or energy=inf
    return PropagateFree(dist);

  Vec3D beta=Velocity();

  double s=beta.Norm();
  if( s == 0 )                       // do not propagate; particle is at rest
    return false;

  Vec3D B_hat = B/b;
  double prm=(beta*B_hat);
  Vec3D beta_par=prm*B_hat;
  Vec3D beta_per=beta-beta_par;
  if( beta_per.Norm2() == 0 )        // magnetic field is parallel
    return PropagateFree(dist);      // to particle's velocity vector

  Vec3D beta_rem=beta_per^B_hat;

  double ct=dist/s;
  double swtD2=sin(wDc*ct/2.);
  double cwtD2=cos(wDc*ct/2.);

  r.r += ct*beta_par+2.*(swtD2/wDc)*(cwtD2*beta_per+swtD2*beta_rem);
  r.r0 += ct;

  p.r=p.r0*(beta_par+cos(wDc*ct)*beta_per+sin(wDc*ct)*beta_rem);
  return true;
}

/* This routine actually propagates the particle from its current position
   to the point of intersection of its trajectory with the plane (future and
   past intersections). The plane is specified by a normal (a,b,c) and a real
   number d s.t. for any point (x,y,z) we have ax+by+cz+d = 0; */

bool Particle::PropagateFreeToPlane(const Vec3D& normal, double d, bool time_reversal_ok)
{
  // Direction of particle velocity
  Vec3D v_hat = p.r / p.r.Norm();

  // Rationalized normal (should probably be the case already)
  Vec3D n_hat = normal / normal.Norm();
  double n_hat_dot_v_hat = n_hat*v_hat;

  if( n_hat_dot_v_hat == 0 )     //the particle never hits the plane
    return false;

  // Compute the distance between the plane and one parallel to it
  // which goes through the particle's current position

  double plane_sep = -d - r.r*normal;

  // Distance the particle must travel to reach the plane i.e.
  // n_hat * vhat = cos( theta ) = plane_dist / propagation_dist
  double propagation_dist = plane_sep / n_hat_dot_v_hat;

  // note: propagation_dist < 0 means particle has passed the plane already
  if((propagation_dist<0)&&(!time_reversal_ok))return false;

  return PropagateFree(propagation_dist);
}

//! Propagates free particle to the closest approach with line
bool Particle::PropagateFreeToPointClosestApproach(const Vec3D& r0,
						   bool time_reversal_ok)
{
  // Direction of particle velocity
  Vec3D v_hat = p.r / p.r.Norm();

  // Distance the particle must travel to reach the closest approach
  double propagation_dist = (r.r-r0)*v_hat;

  // note: propagation_dist < 0 means particle has passed the plane already
  if((propagation_dist<0)&&(!time_reversal_ok))return false;

  return PropagateFree(propagation_dist);
}


//! Propagates free particle to the closest approach with line
bool Particle::PropagateFreeToLineClosestApproach(const Vec3D& normal,
						  const Vec3D& r0,
						  bool time_reversal_ok)
{
  // Direction of particle velocity
  Vec3D v_hat = p.r / p.r.Norm();

  // Rationalized normal (should probably be the case already)
  Vec3D n_hat = normal / normal.Norm();
  double n_hat_dot_v_hat = n_hat*v_hat;

  double denom = 1.0-n_hat_dot_v_hat*n_hat_dot_v_hat;
  if(denom<=0)return true;

  // Distance the particle must travel to reach the closest approach
  double propagation_dist = (r0-r.r)*(n_hat*n_hat_dot_v_hat-v_hat)/denom;

  // note: propagation_dist < 0 means particle has passed the plane already
  if((propagation_dist<0)&&(!time_reversal_ok))return false;

  return PropagateFree(propagation_dist);
}

/*! Propagates the particle to either the closer or farther intersection
    with a given sphere
    \param center Vector containing position of center of the sphere
    \param radius Radius of the sphere
    \param ip Desired point of intersection
    \param time_reversal_ok Allow particle to propagate back in time if
           intersection requires it
*/

Particle::IPOut Particle::
PropagateFreeToSphere(const Vec3D& center, double radius,
		      Particle::IntersectionPoint ip,
		      bool time_reversal_ok)
{
  Vec3D v_hat = p.r / p.r.Norm();  // Normalized direction vector
  Vec3D r_i = r.r;                 // Current (initial) position of particle

  // Vector from here to center of circle
  Vec3D to_center = center - r_i;

  // Distance^2 along particle trajectory to closest approach with center
  double prop_dist = v_hat * to_center;
  double prop_dist_2 = prop_dist * prop_dist;

  // Distance^2 to closest approach with center along tangent (negative inside)
  double tangent_dist_2 = to_center.Norm2() - radius*radius;

  if( prop_dist_2 < tangent_dist_2 )   //no intersection
    return IPO_NONE;

  //time of first intersection
  double t1 = prop_dist - sqrt(prop_dist_2 - tangent_dist_2);
  //time of second intersection (t2>t1)
  double t2 = prop_dist + sqrt(prop_dist_2 - tangent_dist_2);

  double time = 0;

  // Select which intersection point to go to depending on what we were
  // asked to do and whether we can move back in time... otherwise do not
  // propagate

  IPOut ipo = IPO_NONE;
  switch(ip)
    {
    case IP_CLOSEST:
      if(time_reversal_ok)
	if(fabs(t1)<fabs(t2))ipo=IPO_FIRST;
	else ipo=IPO_SECOND;
      else
	if(t1>=0)ipo=IPO_FIRST;
	else if(t2>=0)ipo=IPO_SECOND;
	else return IPO_NONE;
      break;

    case IP_FARTHEST:
      if(time_reversal_ok)
	if(fabs(t1)<fabs(t2))ipo=IPO_SECOND;
	else ipo=IPO_FIRST;
      else
	if(t2>=0)ipo=IPO_SECOND;
	else return IPO_NONE;
      break;

    case IP_NEXT:
      if(t1>0)ipo=IPO_SECOND;
      else if(t2>0)ipo=IPO_SECOND;
      else return IPO_NONE;
      break;

    case IP_PREVIOUS:
      if(!time_reversal_ok)return IPO_NONE;
      else if(t2<0)ipo=IPO_SECOND;
      else if(t1<0)ipo=IPO_FIRST;
      else return IPO_NONE;
      break;

    case IP_EARLIEST:
      if((time_reversal_ok)||(t1>=0))ipo=IPO_FIRST;
      else return IPO_NONE;

    case IP_LATEST:
      if((time_reversal_ok)||(t2>=0))ipo=IPO_SECOND;
      else return IPO_NONE;
    };

  switch(ipo)
    {
    case IPO_FIRST:
      time = t1;
      break;
    case IPO_SECOND:
      time = t2;
      break;
    case IPO_NONE:
      assert(0);
    }

  PropagateFree(time);

  // Since v_hat is of unit length, time corresponds to distance to
  // the intersection

  return ipo;
}

/*! Propagates the particle to either the closer or farther intersection
    with a given cylinder
    \param center Vector containing position of point on axis of cylinder
    \param normal Vector describing axis of cylinder
    \param radius Radius of the cylinder
    \param ip Desired point of intersection
    \param time_reversal_ok Allow particle to propagate back in time if
           intersection requires it
*/

Particle::IPOut Particle::
PropagateFreeToCylinder(const Vec3D& center, const Vec3D& normal,
			double radius,
			IntersectionPoint ip,
			bool time_reversal_ok)
{
  Vec3D n1(normal);
  n1 *= (1.0 / n1.Norm());    // Normalized vector along cylinder axis

  Vec3D n2(p.r);
  n2 *= (1.0 / n2.Norm());    // Normalized direction vector

  Vec3D dr = r.r;             // Current (initial) position of particle
  dr -= center;               // relative to point on cylinder axis

  double n1n2 = n1*n2;
  double drn1 = dr*n1;

  double a = 1.0-n1n2*n1n2;
  double b = 2.0*(dr*n2-drn1*n1n2);
  double c = dr*dr-drn1*drn1-radius*radius;

  double disc = b*b-4.0*a*c;
  if(disc < 0)   //no intersection
    return IPO_NONE;

  double q;
  if(b<0.0) q = -0.5*(b-std::sqrt(disc));
  else q = -0.5*(b+std::sqrt(disc));

  //time of first intersection
  double t1 = c/q;
  //time of second intersection (t2>t1)
  double t2 = q/a;
  if(t2<t1)std::swap(t1,t2);

  double time = 0;

  // Select which intersection point to go to depending on what we were
  // asked to do and whether we can move back in time... otherwise do not
  // propagate

  IPOut ipo = IPO_NONE;
  switch(ip)
  {
  case IP_CLOSEST:
    if(time_reversal_ok) {
      if(fabs(t1)<fabs(t2))ipo=IPO_FIRST;
    	else ipo=IPO_SECOND;
    } else {
      if(t1>=0)ipo=IPO_FIRST;
	    else if(t2>=0)ipo=IPO_SECOND;
	    else return IPO_NONE;
    }
    break;

  case IP_FARTHEST:
    if(time_reversal_ok) {
      if(fabs(t1)<fabs(t2))ipo=IPO_SECOND;
      else ipo=IPO_FIRST;
    } else {
      if(t2>=0)ipo=IPO_SECOND;
      else return IPO_NONE;
    }
    break;

  case IP_NEXT:
    if(t1>0)ipo=IPO_FIRST;
    else if(t2>0)ipo=IPO_SECOND;
    else return IPO_NONE;
    break;

  case IP_PREVIOUS:
    if(!time_reversal_ok)return IPO_NONE;
    else if(t2<0)ipo=IPO_SECOND;
    else if(t1<0)ipo=IPO_FIRST;
    else return IPO_NONE;
    break;

  case IP_EARLIEST:
    if((time_reversal_ok)||(t1>=0))ipo=IPO_FIRST;
    else return IPO_NONE;

  case IP_LATEST:
    if((time_reversal_ok)||(t2>=0))ipo=IPO_SECOND;
    else return IPO_NONE;
  };

  switch(ipo)
  {
  case IPO_FIRST:
    time = t1;
    break;
  case IPO_SECOND:
    time = t2;
    break;
  case IPO_NONE:
    assert(0);
  }

  PropagateFree(time);

  // Since v_hat is of unit length, time corresponds to distance to
  // the intersection

  return ipo;
}


/*** Writes all particle data members into a given stream;
     the individual entries are separated by spaces, and each
     particle in the stream is on one line **********************/
void Particle::WriteParticle(std::ostream& stream)
{
  stream << DoubleToString( r.r0 ) << ' '   // output cT
	 << DoubleToString( r.r.x ) << ' '  // output X
	 << DoubleToString( r.r.y ) << ' '  // output Y
	 << DoubleToString( r.r.z ) << ' '  // output Z
	 << DoubleToString( p.r0 ) << ' '   // output Energy
	 << DoubleToString( p.r.x ) << ' '  // output X-Mom
	 << DoubleToString( p.r.y ) << ' '  // output Y-Mom
	 << DoubleToString( p.r.z ) << ' '  // output Z-Mom
	 << DoubleToString( charge ) << ' ' // output Charge
	 << DoubleToString( m2c4 )          // output Rest energy squared (m2c4)
	 << std::endl;
}

bool Particle::ReadParticle(std::istream &stream)
{
  if( stream.eof() )
    return false;

  stream >> r.r0   // input cT
	 >> r.r.x  // input X
	 >> r.r.y  // input Y
	 >> r.r.z  // input Z
	 >> p.r0   // input Energy
	 >> p.r.x  // input X-Mom
	 >> p.r.y  // input Y-Mom
	 >> p.r.z  // input Z-Mom
	 >> charge //input Charge
	 >> m2c4;  //input m2c4

  return stream.good();
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// The following can be used to make the raw PDG file ready for compilation
//
// cat pdg_particles_02.mc |
//    sed -e 's/\\/\\\\/g; s/["]/\\"/g' |
//    awk '{printf "%-89s\n",$0}' |
//    sed -e 's/^/"/; s/$/\\n"/;'
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

static std::string PDGDataTable2002 =
"* MASSES, WIDTHS, AND MC ID NUMBERS FROM 2002 EDITION OF RPP                             \n"
"*                                                                                        \n"
"* The following values were generated on 06-Jun-2002 by the Berkeley Particle            \n"
"* Data Group from the Review of Particle Properties database and were edited             \n"
"* on 27-Jun-2002 and 9-Jul-2002 to remove the second top mass entry (the indirect        \n"
"* top mass), to add the Omega-minus baryon width (which was omitted from the             \n"
"* generated list for unknown reasons), and to correct a typographical error              \n"
"* on the Omega-minus baryon width line.  Our thanks go to Olivier Schneider              \n"
"* for pointing out these errors.                                                         \n"
"*                                                                                        \n"
"* The values on this list are intended for use in Monte Carlo programs.                  \n"
"*                                                                                        \n"
"* For questions regarding distribution or content of this file, contact                  \n"
"* the Particle Data Group at pdg@lbl.gov.                                                \n"
"*                                                                                        \n"
"* To process the images in this file:                                                    \n"
"* 1) ignore documentation lines that begin with an asterisk                              \n"
"* 2) in a FORTRAN program, process data lines with                                       \n"
"*    FORMAT (BN, A1, 4I8, 1X, E15.0, 2(1X, E8.0), 1X, A21)                               \n"
"* 3) column  1   contains either \"M\" or \"W\" indicating mass or width                 \n"
"*       2 -  9 \\ Monte Carlo particle numbers as described in the \"Review of           \n"
"*      10 - 17 | Particle Properties\". Charge states appear, as appropriate,            \n"
"*      18 - 25 | from left-to-right in the order -, 0, +, ++.                            \n"
"*      26 - 33 /                                                                         \n"
"*           34   blank                                                                   \n"
"*      35 - 49   central value of the mass or width (double precision)                   \n"
"*           50   blank                                                                   \n"
"*      51 - 58   positive error                                                          \n"
"*           59   blank                                                                   \n"
"*      60 - 67   negative error                                                          \n"
"*           68   blank                                                                   \n"
"*      69 - 89   particle name left-justified in the field and                           \n"
"*                charge states right-justified in the field.                             \n"
"*                This field is for ease of visual examination of the file and            \n"
"*                should not be taken as a standardized presentation of                   \n"
"*                particle names.                                                         \n"
"*                                                                                        \n"
"* Particle ID(s)                   Value (GeV)    Errors (GeV)      Name          Charges\n"
"M      22                          0.E+00         +0.0E+00 -0.0E+00 gamma               0\n"
"W      22                          0.E+00         +0.0E+00 -0.0E+00 gamma               0\n"
"M      24                          8.042E+01      +4.0E-02 -4.0E-02 W                   +\n"
"W      24                          2.12E+00       +4.0E-02 -4.0E-02 W                   +\n"
"M      23                          9.11876E+01    +2.1E-03 -2.1E-03 Z                   0\n"
"W      23                          2.4952E+00     +2.3E-03 -2.3E-03 Z                   0\n"
"W      11                          0.E+00         +0.0E+00 -0.0E+00 e                   -\n"
"M      11                          5.10998902E-04 +2.1E-11 -2.1E-11 e                   -\n"
"M      13                          1.05658357E-01 +5.0E-09 -5.0E-09 mu                  -\n"
"W      13                          2.99591E-19    +5.0E-24 -5.0E-24 mu                  -\n"
"M      15                          1.77699E+00    +2.9E-04 -2.6E-04 tau                 -\n"
"W      15                          2.265E-12      +9.0E-15 -9.0E-15 tau                 -\n"
"M      12                          0.E+00         +0.0E+00 -0.0E+00 nu(e)               0\n"
"W      12                          0.E+00         +0.0E+00 -0.0E+00 nu(e)               0\n"
"M      14                          0.E+00         +0.0E+00 -0.0E+00 nu(mu)              0\n"
"W      14                          0.E+00         +0.0E+00 -0.0E+00 nu(mu)              0\n"
"M      16                          0.E+00         +0.0E+00 -0.0E+00 nu(tau)             0\n"
"W      16                          0.E+00         +0.0E+00 -0.0E+00 nu(tau)             0\n"
"M       6                          1.74E+02       +5.0E+00 -5.0E+00 t                   0\n"
"M     211                          1.3957018E-01  +3.5E-07 -3.5E-07 pi                  +\n"
"W     211                          2.5284E-17     +5.0E-21 -5.0E-21 pi                  +\n"
"M     111                          1.349766E-01   +6.0E-07 -6.0E-07 pi                  0\n"
"W     111                          7.8E-09        +6.0E-10 -5.0E-10 pi                  0\n"
"M     221                          5.4730E-01     +1.2E-04 -1.2E-04 eta                 0\n"
"W     221                          1.18E-06       +1.1E-07 -1.1E-07 eta                 0\n"
"M     113     213                  7.711E-01      +9.0E-04 -9.0E-04 rho(770)          0,+\n"
"W     113     213                  1.492E-01      +7.0E-04 -7.0E-04 rho(770)          0,+\n"
"M     223                          7.8257E-01     +1.2E-04 -1.2E-04 omega(782)          0\n"
"W     223                          8.44E-03       +9.0E-05 -9.0E-05 omega(782)          0\n"
"M     331                          9.5778E-01     +1.4E-04 -1.4E-04 eta'(958)           0\n"
"W     331                          2.02E-04       +1.6E-05 -1.6E-05 eta'(958)           0\n"
"M 9010221                          9.80E-01       +1.0E-02 -1.0E-02 f(0)(980)           0\n"
"M 9000111 9000211                  9.847E-01      +1.2E-03 -1.2E-03 a(0)(980)         0,+\n"
"M     333                          1.019456E+00   +2.0E-05 -2.0E-05 phi(1020)           0\n"
"W     333                          4.26E-03       +5.0E-05 -5.0E-05 phi(1020)           0\n"
"M   10223                          1.170E+00      +2.0E-02 -2.0E-02 h(1)(1170)          0\n"
"W   10223                          3.6E-01        +4.0E-02 -4.0E-02 h(1)(1170)          0\n"
"M   10113   10213                  1.2295E+00     +3.2E-03 -3.2E-03 b(1)(1235)        0,+\n"
"W   10113   10213                  1.42E-01       +9.0E-03 -9.0E-03 b(1)(1235)        0,+\n"
"M   20113   20213                  1.23E+00       +4.0E-02 -4.0E-02 a(1)(1260)        0,+\n"
"M     225                          1.2754E+00     +1.2E-03 -1.2E-03 f(2)(1270)          0\n"
"W     225                          1.851E-01      +3.4E-03 -2.6E-03 f(2)(1270)          0\n"
"M   20223                          1.2819E+00     +6.0E-04 -6.0E-04 f(1)(1285)          0\n"
"W   20223                          2.40E-02       +1.2E-03 -1.2E-03 f(1)(1285)          0\n"
"M  100221                          1.293E+00      +5.0E-03 -5.0E-03 eta(1295)           0\n"
"W  100221                          5.5E-02        +5.0E-03 -5.0E-03 eta(1295)           0\n"
"M  100111  100211                  1.30E+00       +1.0E-01 -1.0E-01 pi(1300)          0,+\n"
"M     115     215                  1.3180E+00     +6.0E-04 -6.0E-04 a(2)(1320)        0,+\n"
"W     115     215                  1.07E-01       +5.0E-03 -5.0E-03 a(2)(1320)        0,+\n"
"M   20333                          1.4263E+00     +1.1E-03 -1.1E-03 f(1)(1420)          0\n"
"W   20333                          5.55E-02       +2.9E-03 -2.9E-03 f(1)(1420)          0\n"
"M  100223                          1.419E+00      +3.1E-02 -3.1E-02 omega(1420)         0\n"
"W  100223                          1.7E-01        +6.0E-02 -6.0E-02 omega(1420)         0\n"
"M   10111   10211                  1.474E+00      +1.9E-02 -1.9E-02 a(0)(1450)        0,+\n"
"W   10111   10211                  2.65E-01       +1.3E-02 -1.3E-02 a(0)(1450)        0,+\n"
"M  100113  100213                  1.465E+00      +2.5E-02 -2.5E-02 rho(1450)         0,+\n"
"W  100113  100213                  3.1E-01        +6.0E-02 -6.0E-02 rho(1450)         0,+\n"
"M 9020221                          1.507E+00      +5.0E-03 -5.0E-03 f(0)(1500)          0\n"
"W 9020221                          1.09E-01       +7.0E-03 -7.0E-03 f(0)(1500)          0\n"
"M     335                          1.525E+00      +5.0E-03 -5.0E-03 f(2)'(1525)         0\n"
"W     335                          7.6E-02        +1.0E-02 -1.0E-02 f(2)'(1525)         0\n"
"M   30223                          1.649E+00      +2.4E-02 -2.4E-02 omega(1650)         0\n"
"W   30223                          2.20E-01       +3.5E-02 -3.5E-02 omega(1650)         0\n"
"M     227                          1.667E+00      +4.0E-03 -4.0E-03 omega(3)(1670)      0\n"
"W     227                          1.68E-01       +1.0E-02 -1.0E-02 omega(3)(1670)      0\n"
"M   10115   10215                  1.670E+00      +2.0E-02 -2.0E-02 pi(2)(1670)       0,+\n"
"W   10115   10215                  2.59E-01       +1.0E-02 -1.0E-02 pi(2)(1670)       0,+\n"
"M  100333                          1.680E+00      +2.0E-02 -2.0E-02 phi(1680)           0\n"
"W  100333                          1.5E-01        +5.0E-02 -5.0E-02 phi(1680)           0\n"
"M     117     217                  1.691E+00      +5.0E-03 -5.0E-03 rho(3)(1690)      0,+\n"
"W     117     217                  1.61E-01       +1.0E-02 -1.0E-02 rho(3)(1690)      0,+\n"
"M   30113   30213                  1.700E+00      +2.0E-02 -2.0E-02 rho(1700)         0,+\n"
"W   30113   30213                  2.4E-01        +6.0E-02 -6.0E-02 rho(1700)         0,+\n"
"M   10331                          1.713E+00      +6.0E-03 -6.0E-03 f(0)(1710)          0\n"
"W   10331                          1.25E-01       +1.0E-02 -1.0E-02 f(0)(1710)          0\n"
"M  200111  200211                  1.801E+00      +1.3E-02 -1.3E-02 pi(1800)          0,+\n"
"W  200111  200211                  2.10E-01       +1.5E-02 -1.5E-02 pi(1800)          0,+\n"
"M     337                          1.854E+00      +7.0E-03 -7.0E-03 phi(3)(1850)        0\n"
"W     337                          8.7E-02        +2.8E-02 -2.3E-02 phi(3)(1850)        0\n"
"M  100335                          2.01E+00       +6.0E-02 -8.0E-02 f(2)(2010)          0\n"
"W  100335                          2.0E-01        +6.0E-02 -6.0E-02 f(2)(2010)          0\n"
"M     119     219                  2.011E+00      +1.3E-02 -1.3E-02 a(4)(2040)        0,+\n"
"W     119     219                  3.6E-01        +4.0E-02 -4.0E-02 a(4)(2040)        0,+\n"
"M     229                          2.025E+00      +8.0E-03 -8.0E-03 f(4)(2050)          0\n"
"W     229                          1.94E-01       +1.3E-02 -1.3E-02 f(4)(2050)          0\n"
"M 9060225                          2.297E+00      +2.8E-02 -2.8E-02 f(2)(2300)          0\n"
"W 9060225                          1.5E-01        +4.0E-02 -4.0E-02 f(2)(2300)          0\n"
"M 9070225                          2.34E+00       +6.0E-02 -6.0E-02 f(2)(2340)          0\n"
"W 9070225                          3.2E-01        +8.0E-02 -7.0E-02 f(2)(2340)          0\n"
"M     321                          4.93677E-01    +1.6E-05 -1.6E-05 K                   +\n"
"W     321                          5.315E-17      +1.0E-19 -1.0E-19 K                   +\n"
"M     311                          4.97672E-01    +3.1E-05 -3.1E-05 K                   0\n"
"M     310                          4.97672E-01    +3.1E-05 -3.1E-05 K(S)                0\n"
"M     130                          4.97672E-01    +3.1E-05 -3.1E-05 K(L)                0\n"
"W     310                          7.367E-15      +7.0E-18 -7.0E-18 K(S)                0\n"
"W     130                          1.272E-17      +1.0E-19 -1.0E-19 K(L)                0\n"
"M     323                          8.9166E-01     +2.6E-04 -2.6E-04 K*(892)             +\n"
"M     313                          8.9610E-01     +2.7E-04 -2.7E-04 K*(892)             0\n"
"W     323                          5.08E-02       +9.0E-04 -9.0E-04 K*(892)             +\n"
"W     313                          5.07E-02       +6.0E-04 -6.0E-04 K*(892)             0\n"
"M   10313   10323                  1.272E+00      +7.0E-03 -7.0E-03 K(1)(1270)        0,+\n"
"W   10313   10323                  9.0E-02        +2.0E-02 -2.0E-02 K(1)(1270)        0,+\n"
"M   20313   20323                  1.402E+00      +7.0E-03 -7.0E-03 K(1)(1400)        0,+\n"
"W   20313   20323                  1.74E-01       +1.3E-02 -1.3E-02 K(1)(1400)        0,+\n"
"M  100313  100323                  1.414E+00      +1.5E-02 -1.5E-02 K*(1410)          0,+\n"
"W  100313  100323                  2.32E-01       +2.1E-02 -2.1E-02 K*(1410)          0,+\n"
"M   10311   10321                  1.412E+00      +6.0E-03 -6.0E-03 K(0)*(1430)       0,+\n"
"W   10311   10321                  2.94E-01       +2.3E-02 -2.3E-02 K(0)*(1430)       0,+\n"
"M     325                          1.4256E+00     +1.5E-03 -1.5E-03 K(2)*(1430)         +\n"
"M     315                          1.4324E+00     +1.3E-03 -1.3E-03 K(2)*(1430)         0\n"
"W     325                          9.85E-02       +2.7E-03 -2.7E-03 K(2)*(1430)         +\n"
"W     315                          1.09E-01       +5.0E-03 -5.0E-03 K(2)*(1430)         0\n"
"M   30313   30323                  1.717E+00      +2.7E-02 -2.7E-02 K*(1680)          0,+\n"
"W   30313   30323                  3.2E-01        +1.1E-01 -1.1E-01 K*(1680)          0,+\n"
"M   10315   10325                  1.773E+00      +8.0E-03 -8.0E-03 K(2)(1770)        0,+\n"
"W   10315   10325                  1.86E-01       +1.4E-02 -1.4E-02 K(2)(1770)        0,+\n"
"M     317     327                  1.776E+00      +7.0E-03 -7.0E-03 K(3)*(1780)       0,+\n"
"W     317     327                  1.59E-01       +2.1E-02 -2.1E-02 K(3)*(1780)       0,+\n"
"M   20315   20325                  1.816E+00      +1.3E-02 -1.3E-02 K(2)(1820)        0,+\n"
"W   20315   20325                  2.76E-01       +3.5E-02 -3.5E-02 K(2)(1820)        0,+\n"
"M     319     329                  2.045E+00      +9.0E-03 -9.0E-03 K(4)*(2045)       0,+\n"
"W     319     329                  1.98E-01       +3.0E-02 -3.0E-02 K(4)*(2045)       0,+\n"
"M     411                          1.8693E+00     +5.0E-04 -5.0E-04 D                   +\n"
"W     411                          6.26E-13       +8.0E-15 -8.0E-15 D                   +\n"
"M     421                          1.8645E+00     +5.0E-04 -5.0E-04 D                   0\n"
"W     421                          1.599E-12      +1.0E-14 -1.0E-14 D                   0\n"
"M     423                          2.0067E+00     +5.0E-04 -5.0E-04 D*(2007)            0\n"
"M     413                          2.0100E+00     +5.0E-04 -5.0E-04 D*(2010)            +\n"
"W     413                          9.6E-05        +2.2E-05 -2.2E-05 D*(2010)            +\n"
"M   10423                          2.4222E+00     +1.8E-03 -1.8E-03 D(1)(2420)          0\n"
"W   10423                          1.89E-02       +4.6E-03 -3.5E-03 D(1)(2420)          0\n"
"M     425                          2.4589E+00     +2.0E-03 -2.0E-03 D(2)*(2460)         0\n"
"W     425                          2.3E-02        +5.0E-03 -5.0E-03 D(2)*(2460)         0\n"
"M     415                          2.459E+00      +4.0E-03 -4.0E-03 D(2)*(2460)         +\n"
"W     415                          2.5E-02        +8.0E-03 -7.0E-03 D(2)*(2460)         +\n"
"M     431                          1.9685E+00     +6.0E-04 -6.0E-04 D(s)                +\n"
"W     431                          1.342E-12      +2.6E-14 -2.6E-14 D(s)                +\n"
"M     433                          2.1124E+00     +7.0E-04 -7.0E-04 D(s)*               +\n"
"M   10433                          2.5353E+00     +6.0E-04 -6.0E-04 D(s1)(2536)         +\n"
"M     521                          5.2790E+00     +5.0E-04 -5.0E-04 B                   +\n"
"W     521                          3.93E-13       +4.0E-15 -4.0E-15 B                   +\n"
"M     511                          5.2794E+00     +5.0E-04 -5.0E-04 B                   0\n"
"W     511                          4.27E-13       +4.0E-15 -4.0E-15 B                   0\n"
"M     513     523                  5.3250E+00     +6.0E-04 -6.0E-04 B*                0,+\n"
"M     531                          5.3696E+00     +2.4E-03 -2.4E-03 B(s)                0\n"
"W     531                          4.51E-13       +1.8E-14 -1.8E-14 B(s)                0\n"
"M     541                          6.4E+00        +4.0E-01 -4.0E-01 B(c)                +\n"
"W     541                          1.4E-12        +8.0E-13 -4.0E-13 B(c)                +\n"
"M     441                          2.9797E+00     +1.5E-03 -1.5E-03 eta(c)(1S)          0\n"
"W     441                          1.60E-02       +3.6E-03 -3.2E-03 eta(c)(1S)          0\n"
"M     443                          3.09687E+00    +4.0E-05 -4.0E-05 J/psi(1S)           0\n"
"W     443                          8.7E-05        +5.0E-06 -5.0E-06 J/psi(1S)           0\n"
"M   10441                          3.4151E+00     +8.0E-04 -8.0E-04 chi(c0)(1P)         0\n"
"W   10441                          1.62E-02       +2.3E-03 -2.3E-03 chi(c0)(1P)         0\n"
"M   20443                          3.51051E+00    +1.2E-04 -1.2E-04 chi(c1)(1P)         0\n"
"W   20443                          9.2E-04        +1.3E-04 -1.3E-04 chi(c1)(1P)         0\n"
"M     445                          3.55618E+00    +1.3E-04 -1.3E-04 chi(c2)(1P)         0\n"
"W     445                          2.08E-03       +1.7E-04 -1.7E-04 chi(c2)(1P)         0\n"
"M  100443                          3.68596E+00    +9.0E-05 -9.0E-05 psi(2S)             0\n"
"W  100443                          3.00E-04       +2.5E-05 -2.5E-05 psi(2S)             0\n"
"M   30443                          3.7699E+00     +2.5E-03 -2.5E-03 psi(3770)           0\n"
"W   30443                          2.36E-02       +2.7E-03 -2.7E-03 psi(3770)           0\n"
"M 9000443                          4.040E+00      +1.0E-02 -1.0E-02 psi(4040)           0\n"
"W 9000443                          5.2E-02        +1.0E-02 -1.0E-02 psi(4040)           0\n"
"M 9010443                          4.159E+00      +2.0E-02 -2.0E-02 psi(4160)           0\n"
"W 9010443                          7.8E-02        +2.0E-02 -2.0E-02 psi(4160)           0\n"
"M 9020443                          4.415E+00      +6.0E-03 -6.0E-03 psi(4415)           0\n"
"W 9020443                          4.3E-02        +1.5E-02 -1.5E-02 psi(4415)           0\n"
"M     553                          9.46030E+00    +2.6E-04 -2.6E-04 Upsilon(1S)         0\n"
"W     553                          5.25E-05       +1.8E-06 -1.8E-06 Upsilon(1S)         0\n"
"M   10551                          9.8599E+00     +1.0E-03 -1.0E-03 chi(b0)(1P)         0\n"
"M   20553                          9.8927E+00     +6.0E-04 -6.0E-04 chi(b1)(1P)         0\n"
"M     555                          9.9126E+00     +5.0E-04 -5.0E-04 chi(b2)(1P)         0\n"
"M  100553                          1.002326E+01   +3.1E-04 -3.1E-04 Upsilon(2S)         0\n"
"W  100553                          4.4E-05        +7.0E-06 -7.0E-06 Upsilon(2S)         0\n"
"M  110551                          1.02321E+01    +6.0E-04 -6.0E-04 chi(b0)(2P)         0\n"
"M  120553                          1.02552E+01    +5.0E-04 -5.0E-04 chi(b1)(2P)         0\n"
"M  100555                          1.02685E+01    +4.0E-04 -4.0E-04 chi(b2)(2P)         0\n"
"M  200553                          1.03552E+01    +5.0E-04 -5.0E-04 Upsilon(3S)         0\n"
"W  200553                          2.63E-05       +3.5E-06 -3.5E-06 Upsilon(3S)         0\n"
"M  300553                          1.05800E+01    +3.5E-03 -3.5E-03 Upsilon(4S)         0\n"
"W  300553                          1.4E-02        +5.0E-03 -5.0E-03 Upsilon(4S)         0\n"
"M 9000553                          1.0865E+01     +8.0E-03 -8.0E-03 Upsilon(10860)      0\n"
"W 9000553                          1.10E-01       +1.3E-02 -1.3E-02 Upsilon(10860)      0\n"
"M 9010553                          1.1019E+01     +8.0E-03 -8.0E-03 Upsilon(11020)      0\n"
"W 9010553                          7.9E-02        +1.6E-02 -1.6E-02 Upsilon(11020)      0\n"
"W    2212                          0.E+00         +0.0E+00 -0.0E+00 p                   +\n"
"M    2212                          9.3827200E-01  +4.0E-08 -4.0E-08 p                   +\n"
"M    2112                          9.3956533E-01  +4.0E-08 -4.0E-08 n                   0\n"
"W    2112                          7.432E-28      +7.0E-31 -7.0E-31 n                   0\n"
"M   12112   12212                  1.440E+00      +3.0E-02 -1.0E-02 N(1440)           0,+\n"
"W   12112   12212                  3.5E-01        +1.0E-01 -1.0E-01 N(1440)           0,+\n"
"M    1214    2124                  1.520E+00      +1.0E-02 -5.0E-03 N(1520)           0,+\n"
"W    1214    2124                  1.20E-01       +1.5E-02 -1.0E-02 N(1520)           0,+\n"
"M   22112   22212                  1.535E+00      +2.0E-02 -1.5E-02 N(1535)           0,+\n"
"W   22112   22212                  1.5E-01        +5.0E-02 -5.0E-02 N(1535)           0,+\n"
"M   32112   32212                  1.650E+00      +3.0E-02 -1.0E-02 N(1650)           0,+\n"
"W   32112   32212                  1.50E-01       +4.0E-02 -5.0E-03 N(1650)           0,+\n"
"M    2116    2216                  1.675E+00      +1.0E-02 -5.0E-03 N(1675)           0,+\n"
"W    2116    2216                  1.50E-01       +3.0E-02 -1.0E-02 N(1675)           0,+\n"
"M   12116   12216                  1.680E+00      +1.0E-02 -5.0E-03 N(1680)           0,+\n"
"W   12116   12216                  1.30E-01       +1.0E-02 -1.0E-02 N(1680)           0,+\n"
"M   21214   22124                  1.70E+00       +5.0E-02 -5.0E-02 N(1700)           0,+\n"
"W   21214   22124                  1.0E-01        +5.0E-02 -5.0E-02 N(1700)           0,+\n"
"M   42112   42212                  1.710E+00      +3.0E-02 -3.0E-02 N(1710)           0,+\n"
"W   42112   42212                  1.0E-01        +1.5E-01 -5.0E-02 N(1710)           0,+\n"
"M   31214   32124                  1.720E+00      +3.0E-02 -7.0E-02 N(1720)           0,+\n"
"W   31214   32124                  1.5E-01        +5.0E-02 -5.0E-02 N(1720)           0,+\n"
"M    1218    2128                  2.190E+00      +1.0E-02 -9.0E-02 N(2190)           0,+\n"
"W    1218    2128                  4.5E-01        +1.0E-01 -1.0E-01 N(2190)           0,+\n"
"M    1114    2114    2214    2224  1.2320E+00     +2.0E-03 -2.0E-03 Delta(1232)  -,0,+,++\n"
"W    1114    2114    2214    2224  1.20E-01       +5.0E-03 -5.0E-03 Delta(1232)  -,0,+,++\n"
"M   31114   32114   32214   32224  1.60E+00       +1.0E-01 -5.0E-02 Delta(1600)  -,0,+,++\n"
"W   31114   32114   32214   32224  3.5E-01        +1.0E-01 -1.0E-01 Delta(1600)  -,0,+,++\n"
"M    1112    1212    2122    2222  1.620E+00      +6.0E-02 -5.0E-03 Delta(1620)  -,0,+,++\n"
"W    1112    1212    2122    2222  1.50E-01       +3.0E-02 -3.0E-02 Delta(1620)  -,0,+,++\n"
"M   11114   12114   12214   12224  1.700E+00      +7.0E-02 -3.0E-02 Delta(1700)  -,0,+,++\n"
"W   11114   12114   12214   12224  3.0E-01        +1.0E-01 -1.0E-01 Delta(1700)  -,0,+,++\n"
"M    1116    1216    2126    2226  1.905E+00      +1.5E-02 -3.5E-02 Delta(1905)  -,0,+,++\n"
"W    1116    1216    2126    2226  3.5E-01        +9.0E-02 -7.0E-02 Delta(1905)  -,0,+,++\n"
"M   21112   21212   22122   22222  1.910E+00      +1.0E-02 -4.0E-02 Delta(1910)  -,0,+,++\n"
"W   21112   21212   22122   22222  2.50E-01       +2.0E-02 -6.0E-02 Delta(1910)  -,0,+,++\n"
"M   21114   22114   22214   22224  1.920E+00      +5.0E-02 -2.0E-02 Delta(1920)  -,0,+,++\n"
"W   21114   22114   22214   22224  2.0E-01        +1.0E-01 -5.0E-02 Delta(1920)  -,0,+,++\n"
"M   11116   11216   12126   12226  1.930E+00      +4.0E-02 -1.0E-02 Delta(1930)  -,0,+,++\n"
"W   11116   11216   12126   12226  3.5E-01        +1.0E-01 -1.0E-01 Delta(1930)  -,0,+,++\n"
"M    1118    2118    2218    2228  1.950E+00      +1.0E-02 -1.0E-02 Delta(1950)  -,0,+,++\n"
"W    1118    2118    2218    2228  3.00E-01       +5.0E-02 -1.0E-02 Delta(1950)  -,0,+,++\n"
"M    3122                          1.115683E+00   +6.0E-06 -6.0E-06 Lambda              0\n"
"W    3122                          2.501E-15      +1.9E-17 -1.9E-17 Lambda              0\n"
"M   13122                          1.407E+00      +4.0E-03 -4.0E-03 Lambda(1405)        0\n"
"W   13122                          5.00E-02       +2.0E-03 -2.0E-03 Lambda(1405)        0\n"
"M    3124                          1.5195E+00     +1.0E-03 -1.0E-03 Lambda(1520)        0\n"
"W    3124                          1.56E-02       +1.0E-03 -1.0E-03 Lambda(1520)        0\n"
"M   23122                          1.60E+00       +1.0E-01 -4.0E-02 Lambda(1600)        0\n"
"W   23122                          1.5E-01        +1.0E-01 -1.0E-01 Lambda(1600)        0\n"
"M   33122                          1.670E+00      +1.0E-02 -1.0E-02 Lambda(1670)        0\n"
"W   33122                          3.5E-02        +1.5E-02 -1.0E-02 Lambda(1670)        0\n"
"M   13124                          1.690E+00      +5.0E-03 -5.0E-03 Lambda(1690)        0\n"
"W   13124                          6.0E-02        +1.0E-02 -1.0E-02 Lambda(1690)        0\n"
"M   43122                          1.80E+00       +5.0E-02 -8.0E-02 Lambda(1800)        0\n"
"W   43122                          3.0E-01        +1.0E-01 -1.0E-01 Lambda(1800)        0\n"
"M   53122                          1.81E+00       +4.0E-02 -6.0E-02 Lambda(1810)        0\n"
"W   53122                          1.5E-01        +1.0E-01 -1.0E-01 Lambda(1810)        0\n"
"M    3126                          1.820E+00      +5.0E-03 -5.0E-03 Lambda(1820)        0\n"
"W    3126                          8.0E-02        +1.0E-02 -1.0E-02 Lambda(1820)        0\n"
"M   13126                          1.830E+00      +0.0E+00 -2.0E-02 Lambda(1830)        0\n"
"W   13126                          9.5E-02        +1.5E-02 -3.5E-02 Lambda(1830)        0\n"
"M   23124                          1.890E+00      +2.0E-02 -4.0E-02 Lambda(1890)        0\n"
"W   23124                          1.0E-01        +1.0E-01 -4.0E-02 Lambda(1890)        0\n"
"M    3128                          2.100E+00      +1.0E-02 -1.0E-02 Lambda(2100)        0\n"
"W    3128                          2.0E-01        +5.0E-02 -1.0E-01 Lambda(2100)        0\n"
"M   23126                          2.110E+00      +3.0E-02 -2.0E-02 Lambda(2110)        0\n"
"W   23126                          2.0E-01        +5.0E-02 -5.0E-02 Lambda(2110)        0\n"
"M    3222                          1.18937E+00    +7.0E-05 -7.0E-05 Sigma               +\n"
"W    3222                          8.209E-15      +2.7E-17 -2.7E-17 Sigma               +\n"
"M    3212                          1.192642E+00   +2.4E-05 -2.4E-05 Sigma               0\n"
"W    3212                          8.9E-06        +9.0E-07 -8.0E-07 Sigma               0\n"
"M    3112                          1.197449E+00   +3.0E-05 -3.0E-05 Sigma               -\n"
"W    3112                          4.450E-15      +3.2E-17 -3.2E-17 Sigma               -\n"
"M    3224                          1.3828E+00     +4.0E-04 -4.0E-04 Sigma(1385)         +\n"
"M    3214                          1.3837E+00     +1.0E-03 -1.0E-03 Sigma(1385)         0\n"
"M    3114                          1.3872E+00     +5.0E-04 -5.0E-04 Sigma(1385)         -\n"
"W    3224                          3.58E-02       +8.0E-04 -8.0E-04 Sigma(1385)         +\n"
"W    3214                          3.6E-02        +5.0E-03 -5.0E-03 Sigma(1385)         0\n"
"W    3114                          3.94E-02       +2.1E-03 -2.1E-03 Sigma(1385)         -\n"
"M   13112   13212   13222          1.660E+00      +3.0E-02 -3.0E-02 Sigma(1660)     -,0,+\n"
"W   13112   13212   13222          1.0E-01        +1.0E-01 -6.0E-02 Sigma(1660)     -,0,+\n"
"M   13114   13214   13224          1.670E+00      +1.5E-02 -5.0E-03 Sigma(1670)     -,0,+\n"
"W   13114   13214   13224          6.0E-02        +2.0E-02 -2.0E-02 Sigma(1670)     -,0,+\n"
"M   23112   23212   23222          1.750E+00      +5.0E-02 -2.0E-02 Sigma(1750)     -,0,+\n"
"W   23112   23212   23222          9.0E-02        +7.0E-02 -3.0E-02 Sigma(1750)     -,0,+\n"
"M    3116    3216    3226          1.775E+00      +5.0E-03 -5.0E-03 Sigma(1775)     -,0,+\n"
"W    3116    3216    3226          1.20E-01       +1.5E-02 -1.5E-02 Sigma(1775)     -,0,+\n"
"M   13116   13216   13226          1.915E+00      +2.0E-02 -1.5E-02 Sigma(1915)     -,0,+\n"
"W   13116   13216   13226          1.2E-01        +4.0E-02 -4.0E-02 Sigma(1915)     -,0,+\n"
"M   23114   23214   23224          1.940E+00      +1.0E-02 -4.0E-02 Sigma(1940)     -,0,+\n"
"W   23114   23214   23224          2.2E-01        +8.0E-02 -7.0E-02 Sigma(1940)     -,0,+\n"
"M    3118    3218    3228          2.030E+00      +1.0E-02 -5.0E-03 Sigma(2030)     -,0,+\n"
"W    3118    3218    3228          1.80E-01       +2.0E-02 -3.0E-02 Sigma(2030)     -,0,+\n"
"M    3322                          1.31483E+00    +2.0E-04 -2.0E-04 Xi                  0\n"
"W    3322                          2.27E-15       +7.0E-17 -7.0E-17 Xi                  0\n"
"M    3312                          1.32131E+00    +1.3E-04 -1.3E-04 Xi                  -\n"
"W    3312                          4.02E-15       +4.0E-17 -4.0E-17 Xi                  -\n"
"M    3324                          1.53180E+00    +3.2E-04 -3.2E-04 Xi(1530)            0\n"
"M    3314                          1.5350E+00     +6.0E-04 -6.0E-04 Xi(1530)            -\n"
"W    3324                          9.1E-03        +5.0E-04 -5.0E-04 Xi(1530)            0\n"
"W    3314                          9.9E-03        +1.7E-03 -1.9E-03 Xi(1530)            -\n"
"M   13314   13324                  1.823E+00      +5.0E-03 -5.0E-03 Xi(1820)          -,0\n"
"W   13314   13324                  2.4E-02        +1.5E-02 -1.0E-02 Xi(1820)          -,0\n"
"M    3334                          1.67245E+00    +2.9E-04 -2.9E-04 Omega               -\n"
"W    3334                          8.02E-15       +1.1E-16 -1.1E-16 Omega               -\n"
"M    4122                          2.2849E+00     +6.0E-04 -6.0E-04 Lambda(c)           +\n"
"W    4122                          3.30E-12       +9.0E-14 -9.0E-14 Lambda(c)           +\n"
"M   14122                          2.5939E+00     +8.0E-04 -8.0E-04 Lambda(c)(2593)     +\n"
"W   14122                          3.6E-03        +2.0E-03 -1.3E-03 Lambda(c)(2593)     +\n"
"M    4222                          2.4526E+00     +6.0E-04 -6.0E-04 Sigma(c)(2455)     ++\n"
"M    4212                          2.4513E+00     +7.0E-04 -7.0E-04 Sigma(c)(2455)      +\n"
"M    4112                          2.4522E+00     +6.0E-04 -6.0E-04 Sigma(c)(2455)      0\n"
"W    4222                          2.0E-03        +5.0E-04 -5.0E-04 Sigma(c)(2455)     ++\n"
"W    4112                          1.5E-03        +5.0E-04 -5.0E-04 Sigma(c)(2455)      0\n"
"M    4232                          2.4663E+00     +1.4E-03 -1.4E-03 Xi(c)               +\n"
"W    4232                          1.49E-12       +9.0E-14 -9.0E-14 Xi(c)               +\n"
"M    4132                          2.4718E+00     +1.4E-03 -1.4E-03 Xi(c)               0\n"
"W    4132                          6.7E-12        +1.3E-12 -1.3E-12 Xi(c)               0\n"
"M    4322                          2.5741E+00     +3.3E-03 -3.3E-03 Xi(c)'              +\n"
"M    4312                          2.5788E+00     +3.2E-03 -3.2E-03 Xi(c)'              0\n"
"M    4332                          2.6975E+00     +2.6E-03 -2.6E-03 Omega(c)            0\n"
"W    4332                          1.02E-11       +4.8E-12 -2.4E-12 Omega(c)            0\n"
"M    5122                          5.624E+00      +9.0E-03 -9.0E-03 Lambda(b)           0\n"
"W    5122                          5.36E-13       +3.7E-14 -3.3E-14 Lambda(b)           0\n";

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*! Initialize the particle information list from a file. The file should
  be in the format of the computer readable data tables from the Particle
  Data Group (PDG). At time of writing the 2002 table was available from
  http://pdg.lbl.gov/rpp/mcdata/mass_width_02.mc.

  \param filename: File to read the particle information from. If the
  filename is not supplied then the 2002 list, which is compiled into
  the library, is used.

  \param reinitialize: Force reinitialization of the particle list if
  true. If false and the list is initialized already then this fuction
  returns immediately without overwriting the existing list.
*/

bool Particle::InitializePDGData(const std::string& filename, bool reinitialize)
{
  // Only initialize one time unless user demands more
  if((s_pdgdata_init)&&(!reinitialize))return true;

  // Delete old data
  for(std::map<int, PDGData*>::iterator i=s_pdgdata.begin();
      i!=s_pdgdata.end();i++)
    delete i->second;

  s_pdgdata_init = false;
  s_pdgdata.clear();
  s_pdgdata_name_q.clear();

  // Setup input stream -- either istringstream or ifstream
  std::unique_ptr<std::istream> stream;

  if(filename.length()==0)
    {
      stream.reset(new std::istringstream(PDGDataTable2002));
    }
  else
    {
      stream.reset(new std::ifstream(filename.c_str()));
      if(!(*stream))
	return false;
    }

  // Vector of M-Lines from PDG file and Map of W-Lines for fast lookup
  std::vector<std::string> m_lines;
  std::map<std::string, std::string> w_lines;

  std::string line;
  while(getline(*stream,line))
    {
      switch(line[0])
	{
	case 'M':
	  // Push mass lines onto the back of the vector
	  m_lines.push_back(line);
	  break;
	case 'W':
	  // Store the width lines away keyed by the MC ID entry
	  w_lines[line.substr(1,32)]=line;
	  break;
	case '*':
	  // Skip the comment lines
	  continue;
	default:
	  std::cerr << "Particle::InitialisePDGData: Parse error" << std::endl
		    << line;
	  break;
	};
    }

  // Loop through M-Lines, matching corresponding W-Line and extracting data
  for(std::vector<std::string>::iterator i=m_lines.begin(); i!=m_lines.end(); i++)
    {
      std::string m_line = *i;
      std::string w_line = w_lines[m_line.substr(1,32)];

      std::istringstream str_mcids(m_line.substr(1,32));
      std::istringstream str_mass(m_line.substr(34,33));
      std::istringstream str_name_charge(m_line.substr(68,21));

      double mass;
      double mass_err_lo;
      double mass_err_hi;
      if(!(str_mass >> mass >> mass_err_hi >> mass_err_lo))
	{
	  std::cerr << "Particle::InitialisePDGData: Parse error" << std::endl
		    << "Cannot extract mass terms: \"" << m_line.substr(34,33) << '"'
		    << std::endl << m_line << std::endl;
	  continue;
	}

      double width = -1;
      double width_err_lo = 0;
      double width_err_hi = 0;

      if(w_line.length()==0)
	{
	  if(filename.length() != 0) // only print warning when reading external file
	    std::cerr << "No width for particle: " << m_line << std::endl;
	}
      else
	{
	  std::istringstream str_width(w_line.substr(34,33));
	  if(!(str_width >> width >> width_err_hi >> width_err_lo))
	    {
	      std::cerr << "Particle::InitialisePDGData: Parse error" << std::endl
			<< "Cannot extract width terms: \"" << w_line.substr(34,33) << '"'
			<< std::endl << m_line << std::endl << w_line << std::endl;
	      continue;
	    }
	}

      std::string name;
      std::string charges;
      if(!(str_name_charge >> name >> charges))
	{
	  std::cerr << "Particle::InitialisePDGData: Parse error" << std::endl
		    << "Cannot extract anme and charge terms" << std::endl
		    << std::endl << m_line << std::endl << w_line << std::endl;
	  continue;
	}

      // The PDG table allows multiple (four) particles per table
      // entry, each differing only by charge, (e.g. Delta -,0,+,++ on
      // one line). Each has a different MC ID number and charge but
      // the same mass and width. So loop through the IDs on each
      // line, extrating the charge and make a new PDGData entry for
      // each

      unsigned charge_index=0;
      int id;
      while(str_mcids >> id)
	{
	  int charge=0;
	  bool charge_stop=false;
	  bool charge_good=false;
	  while((charge_index<charges.length())&&(!charge_stop))
	    {
	      switch(charges[charge_index++])
		{
		case '0':
		  charge=0;
		  charge_good=true;
		  break;
		case '+':
		  charge++;
		  charge_good=true;
		  break;
		case '-':
		  charge--;
		  charge_good=true;
		  break;
		case ',':
		  charge_stop=true;
		  break;
		default:
		  charge_good=false;
		  charge_stop=true;
		  break;
		}
	    }

	  if(!charge_good)
	    {
	      std::cerr << "Particle::InitialisePDGData: Parse error" << std::endl
			<< "Cannot extract charge for paricle id:" << id << std::endl
			<< std::endl << m_line << std::endl << w_line << std::endl;
	      m_line.clear();
	      w_line.clear();
	      continue;
	    }

	  PDGData* particle = new PDGData;
	  constexpr double cgs_eV = calin::math::constants::cgs_eV;

	  particle->mcid         = id;
	  particle->name         = name;
	  particle->mass         = mass         * cgs_eV*1e9;
	  particle->mass_err_lo  = mass_err_lo  * cgs_eV*1e9;
	  particle->mass_err_hi  = mass_err_hi  * cgs_eV*1e9;
	  particle->width        = width        * cgs_eV*1e9;
	  particle->width_err_lo = width_err_lo * cgs_eV*1e9;
	  particle->width_err_hi = width_err_hi * cgs_eV*1e9;
	  particle->charge       = charge;

	  s_pdgdata[id]=particle;
          std::pair<std::string, int> key(name,charge);
	  s_pdgdata_name_q[key]  = particle;
	}
    }

  s_pdgdata_init=true;
  return true;
}

std::ostream& operator << (std::ostream& stream,
                           const calin::math::vs_physics::PDGData& p)
{
  std::ostringstream s;
  s << std::setw(7) << std::right << p.mcid << ' '
    << std::setw(15) << std::left << p.name << ' ' << std::showpos
    << std::scientific << std::setprecision(8) << p.mass << ' '
    << std::scientific << std::setprecision(4) << p.mass_err_lo << ' '
    << std::scientific << std::setprecision(4) << p.mass_err_hi << ' '
    // << std::scientific << std::setprecision(8) << p.width << ' '
    // << std::scientific << std::setprecision(4) << p.width_err_lo << ' '
    // << std::scientific << std::setprecision(4) << p.width_err_hi << ' '
    << std::scientific << std::setw(15) << std::setprecision(8) << p.Lifetime() << ' '
    << std::setw(2) << p.charge;
  stream << s.str();
  return stream;
}

std::ostream& operator << (std::ostream& stream,
                           const calin::math::vs_physics::Particle& p)
{
  std::ostringstream s;
  s << p.Position() << ' '
    << std::scientific << std::setprecision(8) << p.Momenta() << ' '
    << std::setprecision(8) << p.p4_p4() << ' '
    << p.MassSquared() << ' '
    << std::fixed << std::setprecision(0) << std::showpos << p.Charge();
  stream << s.str();
  return stream;
}

#ifdef TEST_MAIN

int main()
{
  Particle::InitializePDGData();
  const std::map<int, PDGData*>& table = Particle::GetPDGDataList();
  for(std::map<int, PDGData*>::const_iterator i = table.begin();
      i!=table.end(); i++)std::cout << *(i->second) << std::endl;

  const PDGData* e_data(Particle::GetPDGDataByID(11));
  const PDGData* p_data(Particle::GetPDGDataByNameCharge("p",+1));
  std::cout << std::endl << *e_data << std::endl << *p_data << std::endl << std::endl;

  Particle g(Vec4D(0,0,0,0), Vec3D(1,0,0), .0016021765, 0, 0);   // .001 TeV ph / +x-dir
  Particle e(Vec4D(0,0,0,0), Vec3D(1,0,0), .0016021765, e_data); // .001 TeV e- / +x-dir
  Particle p(Vec4D(0,0,0,0), Vec3D(1,0,0), .0016021765, p_data); // .001 TeV p+ / +x-dir

  std::cout << g << std::endl
	    << e << std::endl
	    << p << std::endl << std::endl;

  g.PropagateFree(1);
  e.PropagateFree(1);
  p.PropagateFree(1);

  std::cout << g << std::endl
	    << e << std::endl
	    << p << std::endl << std::endl;

  g.SetPosition(Vec4D());
  e.SetPosition(Vec4D());
  p.SetPosition(Vec4D());

  g.SetDirection(Vec3D(1,0,0));
  e.SetDirection(Vec3D(1,0,0));
  p.SetDirection(Vec3D(1,0,0));

  for(unsigned i=0;i<10;i++)
    {
      g.PropagateFree(1);
      e.PropagateFree(1);
      p.PropagateFree(1);

      std::ostringstream s;
      s << std::fixed << std::setprecision(6)
	<< g.Position() << ' ' << e.Position() << ' ' << p.Position();
      std::cout << s.str() << std::endl;
    }

  std::cout << std::endl;

  g.SetPosition(Vec4D());
  e.SetPosition(Vec4D());
  p.SetPosition(Vec4D());

  g.SetDirection(Vec3D(1,0,0));
  e.SetDirection(Vec3D(1,0,0));
  p.SetDirection(Vec3D(1,0,0));

  for(unsigned i=0;i<1000;i++)
    {
      g.PropagateBField(1,Vec3D(0,0,1)*100000);
      e.PropagateBField(1,Vec3D(0,0,1)*100000);
      p.PropagateBField(1,Vec3D(0,0,1)*100000);

      std::ostringstream s;
      s << std::fixed << std::setprecision(6)
	<< g.Position() << ' ' << e.Position() << ' ' << p.Position();
      std::cout << s.str() << std::endl;
    }

  std::cout << std::endl;

  Vec3D normal(0,0,1);
  normal.Rotate(Vec3D(0,1,0)*M_PI/10);
  normal.Rotate(Vec3D(1,0,0)*M_PI/10);

  for(int y=0;y<11;y++)
    {
      for(int x=0;x<11;x++)
	{
	  Particle ray(Vec4D(0,double(x-5),double(y-5),0), Vec4D(1,0,0,1), 0);

	  ray.PropagateFreeToPlane(normal,100,true);
	  std::ostringstream s;
	  s << std::fixed << std::setprecision(1) << std::setw(5) << ray.Position().r.z;
	  if(x!=0)std::cout << ' ';
	  std::cout << s.str();
      }
      std::cout << std::endl;
    }
  std::cout << std::endl;

  Vec3D center(0,0,9);

  for(int y=0;y<101;y++)
    {
      for(int x=0;x<11;x++)
	{
	  Vec3D dir(0,0,1); dir.Rotate(Vec3D(1,0,0)*0.1);
	  Particle ray(Vec4D(0,double(x-5),double(y-50)/10,0), Vec4D(1,dir), 0);
	  ray.PropagateFreeToSphere(center,4,false);
	  std::ostringstream s;
	  s << std::fixed << std::setprecision(2) << std::setw(5) << ray.Position().r.z;
	  if(x!=0)std::cout << ' ';
	  std::cout << s.str();
      }
      std::cout << std::endl;
    }
  std::cout << std::endl;

  return 1;
}

#endif // TESTMAIN
