/* 

   calin/simulation/vso_obscuration.cpp -- Stephen Fegan -- 2015-11-09

   Classes for obscurations (such as the camera and telescope arms)

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

#include <math/vs_particle.hpp>
#include <simulation/vso_obscuration.hpp>

using namespace calin::math::vs_physics;
using namespace calin::simulation::vs_optics;

VSOObscuration::~VSOObscuration()
{
  // nothing to see here
}

#if 0
std::vector<VSOObscuration*> 
VSOObscuration::createObsVecFromString(const std::string &str)
{
  std::vector<VSOObscuration*> v;
  std::string s(str);
  while(!s.empty())
    {
      size_t slen = s.length();

      if(s[0] == ',' || s[0] == ' ')
	{
	  s = s.substr(1);
	  continue;
	}

      VSOObscuration* obs = 0;

      obs = VSODiskObscuration::createFromString(s);
      if(obs)v.push_back(obs);
      if(s.length() != slen)continue;

      obs = VSOCylinderObscuration::createFromString(s);
      if(obs)v.push_back(obs);
      if(s.length() != slen)continue;
      
      std::string err("Unknown obscuration type: ");
      err += s;
      std::cerr << err << '\n';
      throw(err);
    }
  return v;
}

std::string 
VSOObscuration::dumpObsVecToString(const std::vector<VSOObscuration*>& vo)
{  
  std::string s;
  for(unsigned i=0;i<vo.size();i++)
    {
      if(i)s += ",";
      s += vo[i]->dumpToString();
    }
  return s;
}

bool VSOObscuration::tokenize(std::string& str, const std::string& name,
			      std::vector<std::string>& tokens)
{
  if(str.substr(0,name.length()) != name)return false;
  if(str[name.length()] != '(')return false;
  size_t iend = str.find(')', name.length()+1);
  if(iend == std::string::npos)
    {
      std::cerr << "Closing parenthesis not found: " << str << '\n';
      str = std::string();
      return false;
    }
  std::string data_str = str.substr(name.length()+1,iend-name.length()-1);
  str = str.substr(iend+1);

  tokens.clear();
  while(!data_str.empty())
    {
      iend = data_str.find(',');
      if(iend == std::string::npos)
	{
	  tokens.push_back(data_str);
	  break;
	}
      else
	{
	  tokens.push_back(data_str.substr(0,iend));
	  data_str = data_str.substr(iend+1);
	}
    }
  return true;
}
#endif

VSODiskObscuration::~VSODiskObscuration()
{
  // nothing to see here
}

bool VSODiskObscuration::doesObscure(const Particle& p_in,
				    Particle& p_out) const
{
  if(fICO && p_in.Velocity().y>0)return false;
  p_out = p_in;
  if(p_out.PropagateFreeToPlane(fN, -fD0, false)
     && (p_out.Position().r-fX0).Norm()<=fR) return true;
   return false;
}

#if 0
VSODiskObscuration* VSODiskObscuration::
createFromString(std::string& str)
{
  std::vector<std::string> tokens;
  if(!tokenize(str,"DISK",tokens))return 0;
  if(tokens.size() != 8)
    {
      std::cerr << "DISK: need exactly 8 arguments (" << tokens.size()
		<< "found) - skipping\n";
      return 0;
    }

  Vec3D X0;
  Vec3D N;
  double r;
  bool ico;
  VSDataConverter::fromString(X0.x, tokens[0]);
  VSDataConverter::fromString(X0.y, tokens[1]);
  VSDataConverter::fromString(X0.z, tokens[2]);
  VSDataConverter::fromString(N.x,  tokens[3]);
  VSDataConverter::fromString(N.y,  tokens[4]);
  VSDataConverter::fromString(N.z,  tokens[5]);
  VSDataConverter::fromString(r,    tokens[6]);
  VSDataConverter::fromString(ico,  tokens[7]);
  return new VSODiskObscuration(X0, N, r, ico);
}

std::string VSODiskObscuration::dumpToString() const
{
  std::string s = "DISK(";
  s += VSDataConverter::toString(fX0.x); s += ",";
  s += VSDataConverter::toString(fX0.y); s += ",";
  s += VSDataConverter::toString(fX0.z); s += ",";
  s += VSDataConverter::toString(fN.x); s += ",";
  s += VSDataConverter::toString(fN.y); s += ",";
  s += VSDataConverter::toString(fN.z); s += ",";
  s += VSDataConverter::toString(fR); s += ",";
  s += VSDataConverter::toString(fICO); s += ")";
  return s;
}
#endif

VSOObscuration* VSODiskObscuration::clone() const
{
  return new VSODiskObscuration(*this);
}


VSOCylinderObscuration::~VSOCylinderObscuration()
{
  // nothing to see here
}

bool VSOCylinderObscuration::doesObscure(const Particle& p_in,
					 Particle& p_out) const
{
  if(fICO && p_in.Velocity().y>0)return false;

  p_out = p_in;
  Particle::IPOut ipo;
  ipo = p_out.PropagateFreeToCylinder(fX1, fN, fR, Particle::IP_NEXT, false);

#if 0
  static unsigned iprint=0;
  static std::ofstream stream;
  if(iprint==0)stream.open("test.dat",std::ofstream::out|std::ofstream::true);
  if(ipo!=Particle::IPO_NONE && iprint<1000)
    { 
      if(iprint==0)stream.open("test.dat",std::ofstream::out|std::ofstream::trunc);
      stream << ipo << ' ' << p_in.Position() << ' '
	     << p_out.Position().r << ' ' << p_out.Position().r*fN << ' '
	     << fD1 << ' ' << fD2 << '\n';
      iprint++;
      if(iprint==1000)stream.close();
    }
#endif

  if(ipo != Particle::IPO_NONE)
    {
      double Dc = p_out.Position().r*fN;
      if((std::fabs(Dc-fD1)<=fD)&&(std::fabs(Dc-fD2)<=fD))return true;
      if(ipo == Particle::IPO_SECOND)return false;
      ipo = 
	p_out.PropagateFreeToCylinder(fX1, fN, fR, Particle::IP_LATEST, false);
      if(ipo != Particle::IPO_SECOND)return false;
      Dc = p_out.Position().r*fN;
      if((std::fabs(Dc-fD1)<=fD)&&(std::fabs(Dc-fD2)<=fD))return true;
    }
  
  return false;
}

#if 0
VSOCylinderObscuration* VSOCylinderObscuration::
createFromString(std::string& str)
{
  std::vector<std::string> tokens;
  if(!tokenize(str,"TUBE",tokens))return 0;
  if(tokens.size() != 8)
    {
      std::cerr << "TUBE: need exactly 8 arguments (" << tokens.size()
		<< "found)\n";
      return 0;
    }
  Vec3D X0;
  Vec3D X1;
  double r;
  bool ico;
  VSDataConverter::fromString(X0.x, tokens[0]);
  VSDataConverter::fromString(X0.y, tokens[1]);
  VSDataConverter::fromString(X0.z, tokens[2]);
  VSDataConverter::fromString(X1.x, tokens[3]);
  VSDataConverter::fromString(X1.y, tokens[4]);
  VSDataConverter::fromString(X1.z, tokens[5]);
  VSDataConverter::fromString(r,    tokens[6]);
  VSDataConverter::fromString(ico,  tokens[7]);
  return new VSOCylinderObscuration(X0, X1, r, ico);
}

std::string VSOCylinderObscuration::dumpToString() const
{
  std::string s("TUBE(");
  s += VSDataConverter::toString(fX1.x); s += ",";
  s += VSDataConverter::toString(fX1.y); s += ",";
  s += VSDataConverter::toString(fX1.z); s += ",";
  s += VSDataConverter::toString(fX2.x); s += ",";
  s += VSDataConverter::toString(fX2.y); s += ",";
  s += VSDataConverter::toString(fX2.z); s += ",";
  s += VSDataConverter::toString(fR); s += ",";
  s += VSDataConverter::toString(fICO); s += ")";
  return s;
}
#endif

VSOObscuration* VSOCylinderObscuration::clone() const
{
  return new VSOCylinderObscuration(*this);
}
