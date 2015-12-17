/* 

   calin/math/vs_particle.hpp -- Stephen Fegan -- 2015-11-05

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

/*! \file Particle.hpp
  relativistic particle class header file

  \author   Stephen Fegan             \n
            UCLA                      \n
	    sfegan@astro.ucla.edu     \n
  
  \author   Vladimir Vassiliev        \n
            UCLA                      \n
	    vvv@astro.ucla.edu        \n

  \author   Maciej Nicewicz           \n
            UCLA                      \n
	    nicewicz@physics.ucla.edu \n


  \date     09/13/2004
  \version  1.1
  \note
*/

#pragma once

#include <fstream>
#include <string>
#include <utility>
#include <map>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <math/constants.hpp>
#include <math/vs_vec4d.hpp>

namespace calin { namespace math { namespace vs_physics {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Particle Data Group list of particle mass, width and charges
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*! \class PDGData
  \brief Container for Particle Data Group mass and width information
*/
class PDGData
{
 public:
  PDGData(): 
      mcid(), name(), 
      mass(), mass_err_lo(), mass_err_hi(),
      width(), width_err_lo(), width_err_hi(), charge()
  { /* nothing to see here */ }
  
  int mcid;                // PDG Monte Carlo ID (see PDG review)
  std::string name;        // PDG (non-authoritative) name
  double mass;             // ergs
  double mass_err_lo;      // ergs
  double mass_err_hi;      // ergs
  double width;            // ergs
  double width_err_lo;     // ergs
  double width_err_hi;     // ergs
  int charge;              // e  --  we do not handle fractional charges

  inline double Lifetime() const; // resonance lifetime [s]
};

/*! \class Particle
  \brief relativistic particle class
*/
  
class Particle
{
 public:
  inline Particle();                                 //!<Construct massless particle at rest
  inline Particle(const Particle& o);                //!<Copy construct
  
  inline Particle(double mc2, double q);             //!<Construct massive particle at rest
  inline Particle(const PDGData* data);              //!<Construct massive particle at rest
  
  inline Particle(const Vec4D& _r, const Vec4D& _p, 
                  double q);                         //!<Construct from 4-pos, 4-mom and charge
  
  inline Particle(const Vec4D& _r, const Vec3D& v_hat, double E,
                  double mc2, double q);             //!<Construct from 4-pos, 3-dir, energy, rest mass and charge
  inline Particle(const Vec4D& _r, const Vec3D& v_hat, double E,
                  const PDGData* data);              //!<Construct from 4-pos, 3-dir, energy and PDG data
  
  inline const Vec4D& Position() const;   //!<accessor
  inline const Vec4D& Momenta() const;    //!<accessor
  inline const Vec3D Velocity() const;   //!<particle's velocity  v/c  [1]
  inline double MassSquared() const;      //!<accessor
  inline double Charge() const;           //!<accessor
    
  inline void SetPosition( const Vec4D& _r ); //!<accessor 
  inline void SetMomenta( const Vec4D& _p );  //!<accessor 
  inline void SetDirection( const Vec3D& v_hat); //!<accesor
  inline void SetCharge( double _q );         //!<accessor
    
  inline void Rotate(const Vec3D& axis); //!<rotates particle vectors around axis
  inline bool Boost(const Vec3D& axis);  //!<boosts 4-vectors
  inline void TranslateOrigin(const Vec4D& new_origin);  //!<translates position 4-vector
  inline void P();                       //!<parity transformation
  inline void T();                       //!<time inversion
  inline void Reflect(const Vec3D& norm);  //!< reflect vector

  inline void Repair();                  //!<restores p4 on mass shell
    
  inline double r4_r4() const;           //!<calculates 4-position square
  inline double p4_p4() const;           //!<calculates 4-momenta square
  
  //! Propagate free particle fixed distance
  inline bool PropagateFree(double dist);
  
  //! Propagate particle fixed distance in magnetic field
  bool PropagateBField(double dist, const Vec3D& B);
  
  //! Propagates free particle to the given plane
  bool PropagateFreeToPlane(const Vec3D& normal, double d, 
                            bool time_reversal_ok=true); 

  //! Propagates free particle to the closest approach with line
  bool PropagateFreeToPointClosestApproach(const Vec3D& r0, 
                                           bool time_reversal_ok=true); 

  //! Propagates free particle to the closest approach with line
  bool PropagateFreeToLineClosestApproach(const Vec3D& normal,
                                          const Vec3D& r0, 
                                          bool time_reversal_ok=true); 

  //! Propagates particle to a given sphere
  enum IntersectionPoint { IP_CLOSEST, IP_FARTHEST, 
                           IP_NEXT, IP_PREVIOUS, 
                           IP_EARLIEST, IP_LATEST };
  
  enum IPOut { IPO_NONE, IPO_FIRST, IPO_SECOND };
  
  IPOut PropagateFreeToSphere(const Vec3D& center, double radius, 
                              IntersectionPoint ip = IP_CLOSEST,
                              bool time_reversal_ok = true);
  
  IPOut PropagateFreeToCylinder(const Vec3D& center, const Vec3D& normal, 
                                double radius, 
                                IntersectionPoint ip = IP_CLOSEST,
                                bool time_reversal_ok = true);
    
  void Dump(std::ostream& stream=std::cout);
    
  void WriteParticle(std::ostream &stream);  //!<writes particle data to a stream
  bool ReadParticle(std::istream &stream);  //!<initializes the particle from a stream  

  inline static bool InitializePDGData(const std::string& filename=""); //!< Initialize PDG data
  static bool InitializePDGData(const std::string& filename,bool reinitialize); //!< Initialize or reinitialize PDG data
  inline static const std::map<int, PDGData*>& GetPDGDataList() { return s_pdgdata; } //!< Return all particle information
  inline static const PDGData* GetPDGDataByID(int id); //!< Get one particle by its ID
  inline static const PDGData* GetPDGDataByNameCharge(const std::string& name, int charge); //!< Get one particle by its name and charge
    
 private:
  Vec4D  r;        //!<4-position: (c*t, x, y, z)                   [cm]
  Vec4D  p;        //!<4-momenta:  (e, px*c, py*c, pz*c)           [erg]
  double m2c4;     //!<rest (m*c^2)^2 of particle                [erg^2]
  double charge;   //!<particle's charge in units of electron charge [e]    

  static bool                                              s_pdgdata_init;
  static std::map<int, PDGData*>                           s_pdgdata;
  static std::map<std::pair<std::string, int>, PDGData*>   s_pdgdata_name_q;
};

//  std::string DoubleToString( double ); //!<converts double to string

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** Default class constructor. Construct a massless particle at rest
    at the origin
*/  
inline Particle::Particle(): r(), p(), m2c4(p*p), charge()
{
  // nothing to see here
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//! Copy constructor
inline Particle::Particle(const Particle& o):
    r(o.r), p(o.p), m2c4(o.m2c4), charge(o.charge)
{
  // nothing to see here
}
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** Overloaded class constructor. Construct a particle at rest
    at the origin.
    
    \param mc2: rest mass                          [erg]
    \param   q: charge in units of electron charge   [e]
*/
inline Particle::Particle (double mc2, double q): 
    r(), p(mc2,Vec3D()), m2c4(mc2*mc2), charge(q)
{
  // nothing to see here
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** Overloaded class constructor. Construct a particle at rest
    at the origin.

    \param data: PDG data for particle
*/
inline Particle::Particle (const PDGData* data):
    r(), p(data->mass,Vec3D()), m2c4(data->mass*data->mass), charge(data->charge)
{
  // nothing to see here
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** Overloaded class constructor. Construct a particle with
    specified 4-position, 4-momentum and charge.
    
    \param _r: 4-position: (c*t, x, y, z)          [cm]
    \param _p: 4-momenta:  (e, px*c, py*c, pz*c)  [erg]
    \param  q: charge in units of electron charge   [e] 
*/
inline Particle::Particle(const Vec4D& _r, const Vec4D& _p, double q):
    r(_r), p(_p), m2c4(p*p), charge(q)
{
  // nothing to see here
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** Overloaded class constructor. Construct a particle with
    specified 4-position, 3-velocity direction, energy, 
    rest mass and charge
    
    \param  _r: 4-position: (c*t, x, y, z)            [cm]
    \param   v: direction of velocity           [unitless]
    \param   E: energy                               [erg]
    \param mc2: rest mass                            [erg]
    \param   q: charge in units of electron charge     [e] 
*/
inline Particle::Particle(const Vec4D& _r, const Vec3D& v_hat, double E, 
                          double mc2, double q):
    r(_r), p(E,v_hat*(sqrt(E*E-mc2*mc2)/v_hat.Norm())), 
    m2c4(mc2*mc2), charge(q)
{
  // nothing to see here
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** Overloaded class constructor. Construct a particle with
    specified 4-position, 3-velocity direction, energy, 
    rest mass and charge
    
    \param   _r: 4-position: (c*t, x, y, z)            [cm]
    \param    v: direction of velocity           [unitless]
    \param    E: energy                               [erg]
    \param data: PDG data for particle
*/
inline Particle::Particle(const Vec4D& _r, const Vec3D& v_hat, double E, 
                          const PDGData* data):
    r(_r), p(E,v_hat*(sqrt(E*E-data->mass*data->mass)/v_hat.Norm())), 
    m2c4(data->mass*data->mass), charge(double(data->charge))
{
  // nothing to see here
}
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Accessor
inline const Vec4D& Particle::Position() const
{
  return r;
}
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Accessor
inline const Vec4D& Particle::Momenta() const
{
  return p;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Returns 3D velocity vector  v/c
inline const Vec3D Particle::Velocity() const
{
  return p.r/p.r0;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Accessor
inline double Particle::MassSquared() const
{
  return m2c4;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Accessor
inline double Particle::Charge() const
{
  return charge;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Accessor
/// \param r4: 4-position: (c*t, x, y, z)          [cm]
inline void Particle::SetPosition(const Vec4D& _r)
{
  r=_r;
  return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Accessor
/// \param p4: 4-momenta:  (e, px*c, py*c, pz*c)  [erg]
inline void Particle::SetMomenta(const Vec4D& _p)
{
  p=_p;
  return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Reset the direction of propagation, keeping energy fixed
/// \param v_hat: Direction of proagation, i.e. (vx/v, vy/v, vz/v)
inline void Particle::SetDirection( const Vec3D& v_hat)
{
  p.r = v_hat*(sqrt(p.r0*p.r0-m2c4)/v_hat.Norm());
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Accessor
/// \param q: charge in units of electron charge   [e]
inline void Particle::SetCharge(double q)
{
  charge=q;
  return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Rotation around vector "axis".
/// \param axis: axis of rotation vector [rad] 
/*! \note 
  The modulus of rotation angle is equal to the axis 
  norm and is given in radians. The rotation of e_x 
  around e_z by PI/2 is equal to e_y.
*/

inline void Particle::Rotate(const Vec3D& axis)
{
  r.Rotate(axis);
  p.Rotate(axis);
  return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Boost by vector v/c
/// \param boost: boost vector (boost,boost) < 1.0 
/*! \note 
  
 */
inline bool Particle::Boost(const Vec3D& axis)
{
  if(!r.Boost(axis))return false;   // boost r4
  if(!p.Boost(axis))return false;   // boost p4
  return true;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Translate position 4-vector
/// \param trst: translation 4-vector 
inline void Particle::TranslateOrigin(const Vec4D& new_origin)
{
  r-=new_origin;
  return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Parity transformation
inline void Particle::P()
{
  r.P();
  p.P();
  return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Time inversion
inline void Particle::T()
{
  r.T();
  p.T();
  return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for reflection in mirror (defined by norm)
inline void Particle::Reflect(const Vec3D& norm)
{
  p.Reflect(norm);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Restores p4 on mass shell.  
/*! \note 
  During calculations with large relativistic factors an 
  accuracy of (p4,p4) can be lost. This method recalculates 
  energy to restore particle's mass. 
*/
inline void Particle::Repair()
{
  p.r0=sqrt(m2c4+(p.r*p.r));
  return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// 4-position square
inline double Particle::r4_r4() const
{
  return r*r;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// 4-momenta square
inline double Particle::p4_p4() const
{
  return p*p;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*! Propagate particle by positive dist
  \param dist: propagation distance        [cm]
  
  \note 
  Returns 'true' if particle is propagated, and returns 'false'
  if particle is at rest
*/
inline bool Particle::PropagateFree(double dist)
{
  Vec3D beta=Velocity();
  double s=beta.Norm();
  
  if( s == 0 )      // do not propagate; particle is at rest       
    return false;
  
  double ct=dist/s;
  r.r+=ct*beta;
  r.r0+=ct;
  
  return true;
}
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*! Initialize the particle information list from a file. The file should
  be in the format of the computer readable data tables from the Particle
  Data Group (PDG). At time of writing the 2002 table was available from
  http://pdg.lbl.gov/rpp/mcdata/mass_width_02.mc. This function is just
  an inline call to InitializePDGData(const std::string&,bool) which
  tests whether or not the particle list is initialized first. This may
  save some time.
  
  \param filename: File to read the particle information from. If it is
  blank then the 2002 list, which is compiled into the library, is used.
*/
inline bool Particle::InitializePDGData(const std::string& filename)
{
  if(!s_pdgdata_init)return InitializePDGData(filename,true);     
  else return true;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*! Get information on one particle from the list. The particle is specified
  by the PDG Monte Carlo ID number. 
      
  \param id: PDG MC number, see "Monte Carlo Particle Numbering Scheme"
  review in the PDG Review of Particle Physics.
*/
inline const PDGData* Particle::GetPDGDataByID(int id)
{
  InitializePDGData();
  return s_pdgdata[id];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*! Get information on one particle from the list. The particle is specified
  by the PDG Monte Carlo ID number. 
  
  \param id: PDG MC number, see "Monte Carlo Particle Numbering Scheme"
  review in the PDG Review of Particle Physics.
*/
inline const PDGData* 
Particle::GetPDGDataByNameCharge(const std::string& name, int charge)
{
  InitializePDGData();
  std::pair<std::string, int> key(name,charge);
  return s_pdgdata_name_q[key];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*! Return the lifetime of the particle (resonance) in seconds
 */
inline double PDGData::Lifetime() const
{ 
  return calin::math::constants::cgs_hbar/width; 
}

} } } // namespace calin::math::vs_physics

#ifndef SWIG
std::ostream& operator << (std::ostream& stream,
                           const calin::math::vs_physics::PDGData& p);
std::ostream& operator << (std::ostream& stream,
                           const calin::math::vs_physics::Particle& p);
#endif
