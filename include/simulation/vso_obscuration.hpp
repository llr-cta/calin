/* 

   calin/simulation/vso_obscuration.hpp -- Stephen Fegan -- 2015-11-09

   Class for telescope data and functionality. This code is derived
   from raytracing code largely implemented by the author at UCLA in
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

#include <math/vs_vec3d.hpp>
#include <math/vs_particle.hpp>
#include <simulation/vs_optics.pb.h>

namespace calin { namespace simulation { namespace vs_optics {

class VSOObscuration
{
 public:
  virtual ~VSOObscuration();
  virtual bool doesObscure(const math::vs_physics::Particle& p_in)
  {
    math::vs_physics::Particle p_out;
    return doesObscure(p_in, p_out);
  }
  virtual bool doesObscure(const math::vs_physics::Particle& p_in,
                           math::vs_physics::Particle& p_out) const = 0;

  virtual VSOObscuration* clone() const = 0;

  virtual void
  dump_to_proto(ix::simulation::vs_optics::VSOObscurationData* d) const = 0;
  
  static VSOObscuration*
  create_from_proto(const ix::simulation::vs_optics::VSOObscurationData& d);
  
#if 0
  virtual std::string dumpToString() const = 0;
  static std::vector<VSOObscuration*> 
  createObsVecFromString(const std::string &str);
  static std::string 
  dumpObsVecToString(const std::vector<VSOObscuration*>& vo);
  
  static bool tokenize(std::string& str, const std::string& name,
                       std::vector<std::string>& tokens);
#endif
};

class VSODiskObscuration: public VSOObscuration
{
 public:
  VSODiskObscuration(const math::vs_physics::Vec3D& center,
                     const math::vs_physics::Vec3D& normal,
                     double radius, bool incoming_only):
      VSOObscuration(), fX0(center), fN(normal), fR(radius), fD0(),
      fICO(incoming_only)
  {
    fD0 = center*normal;
  }
  virtual ~VSODiskObscuration();
  bool doesObscure(const math::vs_physics::Particle& p_in,
                   math::vs_physics::Particle& p_out) const override;
  VSOObscuration* clone() const override;

  void
  dump_to_proto(ix::simulation::vs_optics::VSOObscurationData* d)
      const override;
  
  static VSOObscuration*
  create_from_proto(const ix::simulation::vs_optics::VSODiskObscurationData& d);

#if 0
  virtual std::string dumpToString() const;
  static VSODiskObscuration* createFromString(std::string& str);
#endif

  const math::vs_physics::Vec3D& center_pos() const { return fX0; }
  const math::vs_physics::Vec3D& normal() const { return fN; }
  double diameter() const { return 2.0*fR; }
  bool incoming_only() const { return fICO; }

 private:
  math::vs_physics::Vec3D fX0;
  math::vs_physics::Vec3D fN;
  double                  fR;
  double                  fD0;
  bool                    fICO;
};

class VSOTubeObscuration: public VSOObscuration
{
 public:
  VSOTubeObscuration(const math::vs_physics::Vec3D& x1,
                         const math::vs_physics::Vec3D& x2,
                         double radius, bool incoming_only):
      VSOObscuration(), fX1(x1), fX2(x2), fR(radius), fN(), fD1(), fD2(),
      fICO(incoming_only)
  { 
    fN = x2-x1;
    fN /= fN.Norm();
    fD1 = fN*x1;
    fD2 = fN*x2;
    fD  = std::fabs(fD2-fD1);
  }

  virtual ~VSOTubeObscuration();
  bool doesObscure(const math::vs_physics::Particle& p_in,
                   math::vs_physics::Particle& p_out) const override;
  VSOObscuration* clone() const override;

  void dump_to_proto(ix::simulation::vs_optics::VSOObscurationData* d)
      const override;
  
  static VSOObscuration*
  create_from_proto(const ix::simulation::vs_optics::VSOTubeObscurationData& d);

#if 0
  virtual std::string dumpToString() const;
  static VSOTubeObscuration* createFromString(std::string& str);
#endif

  const math::vs_physics::Vec3D& end1_pos() const { return fX1; }
  const math::vs_physics::Vec3D& end2_pos() const { return fX2; }
  double diameter() const { return 2.0*fR; }
  bool incoming_only() const { return fICO; }
  
 private:
  math::vs_physics::Vec3D fX1;
  math::vs_physics::Vec3D fX2;
  double                  fR;
  math::vs_physics::Vec3D fN;
  double                  fD1;
  double                  fD2;
  double                  fD;
  bool                    fICO;
};

} } } // namespace calin::simulation::vs_optics

