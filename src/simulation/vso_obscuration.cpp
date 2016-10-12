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

VSOObscuration* VSOObscuration::
create_from_proto(const ix::simulation::vs_optics::VSOObscurationData& d)
{
  if(d.has_disk())return VSODiskObscuration::create_from_proto(d.disk());
  else if(d.has_tube())return VSOTubeObscuration::create_from_proto(d.tube());
  else if(d.has_aligned_box())return VSOAlignedBoxObscuration::create_from_proto(d.aligned_box());
  else {
    throw std::runtime_error("VSOObscuration::create_from_proto: unknown obscuration type");
    return 0;
  }
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

      obs = VSOTubeObscuration::createFromString(s);
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

calin::ix::simulation::vs_optics::VSOObscurationData* VSODiskObscuration::
dump_as_proto(ix::simulation::vs_optics::VSOObscurationData* d) const
{
  if(d == nullptr)d = new ix::simulation::vs_optics::VSOObscurationData;
  auto* dd = d->mutable_disk();
  fX0.dump_as_proto(dd->mutable_center_pos());
  fN.dump_as_proto(dd->mutable_normal());
  dd->set_diameter(2.0*fR);
  dd->set_incoming_only(fICO);
  return d;
}

VSODiskObscuration* VSODiskObscuration::
create_from_proto(const ix::simulation::vs_optics::VSODiskObscurationData& d)
{
  return new VSODiskObscuration(d.center_pos(), d.normal(),
                                d.diameter()/2.0, d.incoming_only());
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

VSODiskObscuration* VSODiskObscuration::clone() const
{
  return new VSODiskObscuration(*this);
}


VSOTubeObscuration::~VSOTubeObscuration()
{
  // nothing to see here
}

bool VSOTubeObscuration::doesObscure(const Particle& p_in,
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

calin::ix::simulation::vs_optics::VSOObscurationData* VSOTubeObscuration::
dump_as_proto(ix::simulation::vs_optics::VSOObscurationData* d) const
{
  if(d == nullptr)d = new calin::ix::simulation::vs_optics::VSOObscurationData;
  auto* dd = d->mutable_tube();
  fX1.dump_as_proto(dd->mutable_end1_pos());
  fX2.dump_as_proto(dd->mutable_end2_pos());
  dd->set_diameter(2.0*fR);
  dd->set_incoming_only(fICO);
  return d;
}

VSOTubeObscuration* VSOTubeObscuration::
create_from_proto(const ix::simulation::vs_optics::VSOTubeObscurationData& d)
{
  return new VSOTubeObscuration(d.end1_pos(), d.end2_pos(),
                                d.diameter()/2.0, d.incoming_only());
}

#if 0
VSOTubeObscuration* VSOTubeObscuration::
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
  return new VSOTubeObscuration(X0, X1, r, ico);
}

std::string VSOTubeObscuration::dumpToString() const
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

VSOTubeObscuration* VSOTubeObscuration::clone() const
{
  return new VSOTubeObscuration(*this);
}


VSOAlignedBoxObscuration::~VSOAlignedBoxObscuration()
{
  // nothing to see here
}

bool VSOAlignedBoxObscuration::doesObscure(
  const math::vs_physics::Particle& p_in,
  math::vs_physics::Particle& p_out) const
{
  // See: https://tavianator.com/fast-branchless-raybounding-box-intersections/
  // and: https://tavianator.com/fast-branchless-raybounding-box-intersections-part-2-nans/

  p_out = p_in;

  // Normalized direction vector
  Vec3D v_hat = p_out.Velocity() / p_out.Velocity().Norm();

  const double tx1 = (min_corner_.x - p_out.Position().r.x)/v_hat.x;
  const double tx2 = (max_corner_.x - p_out.Position().r.x)/v_hat.x;
  double tmin = std::min(tx1, tx2);
  double tmax = std::max(tx1, tx2);

  const double ty1 = (min_corner_.y - p_out.Position().r.y)/v_hat.y;
  const double ty2 = (max_corner_.y - p_out.Position().r.y)/v_hat.y;
  tmin = std::max(tmin, std::min(std::min(ty1, ty2), tmax));
  tmax = std::min(tmax, std::max(std::max(ty1, ty2), tmin));

  const double tz1 = (min_corner_.z - p_out.Position().r.z)/v_hat.z;
  const double tz2 = (max_corner_.z - p_out.Position().r.z)/v_hat.z;
  tmin = std::max(tmin, std::min(std::min(tz1, tz2), tmax));
  tmax = std::min(tmax, std::max(std::max(tz1, tz2), tmin));

  if(tmax > std::max(tmin, 0.0)) {
    if(tmin > 0)p_out.PropagateFree(tmin);
    return true;
  }

  return false;
}

VSOAlignedBoxObscuration* VSOAlignedBoxObscuration::clone() const
{
  return new VSOAlignedBoxObscuration(*this);
}

calin::ix::simulation::vs_optics::VSOObscurationData*
VSOAlignedBoxObscuration::dump_as_proto(
  calin::ix::simulation::vs_optics::VSOObscurationData* d) const
{
  if(d == nullptr)d = new calin::ix::simulation::vs_optics::VSOObscurationData;
  auto* dd = d->mutable_aligned_box();
  max_corner_.dump_as_proto(dd->mutable_max_corner());
  min_corner_.dump_as_proto(dd->mutable_min_corner());
  dd->set_incoming_only(incoming_only_);
  return d;
}

VSOAlignedBoxObscuration* VSOAlignedBoxObscuration::create_from_proto(
  const ix::simulation::vs_optics::VSOAlignedBoxObscurationData& d)
{
  return new VSOAlignedBoxObscuration(d.max_corner(), d.min_corner(),
    d.incoming_only());
}
