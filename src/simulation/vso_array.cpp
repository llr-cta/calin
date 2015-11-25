/* 

   calin/simulation/vso_array.cpp -- Stephen Fegan -- 2015-11-25

   Class for array of telescopes

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

/*! \file ArrayParameters.cpp
  Array code file

  \author  Stephen Fegan               \n
           UCLA                        \n
	   sfegan@astro.ucla.edu       \n

  \date    12/05/2004
  \version 0.2
  \note
*/

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <set>

#include <simulation/vso_array.hpp>

using namespace calin::simulation::vs_optics;
using namespace calin::ix::simulation::vs_optics;

VSOArray::VSOArray():
    fLatitude(), fLongitude(), fAltitude(), fSpacing(), fArrayParity(), 
    fTelescopes(), fTelescopesByHexID()
{
  // nothing to see here
}

VSOArray::~VSOArray()
{
  for(std::vector<VSOTelescope*>::iterator i = fTelescopes.begin();
      i!=fTelescopes.end(); i++)delete *i;
}

// ****************************************************************************
// General functions
// ****************************************************************************

bool VSOArray::pointTelescopes(const math::vs_physics::Vec3D& v)
{
  if(v.Norm2() == 0 )
    return false;

  bool good = true;
  for(std::vector<VSOTelescope*>::iterator i = fTelescopes.begin();
      i!=fTelescopes.end(); i++)good &= (*i)->pointTelescope(v);

  return good;
}

bool VSOArray::
pointTelescopesAzEl(const double az_rad, const double el_rad)
{
  bool good = true;
  for(std::vector<VSOTelescope*>::iterator i = fTelescopes.begin();
      i!=fTelescopes.end(); i++)
    good &= (*i)->pointTelescopeAzEl(az_rad,el_rad);
  return good;
}

// ****************************************************************************
// Array creation 
// ****************************************************************************

#if 0
void VSOArray::
generateFromArrayParameters(const IsotropicDCArrayParameters& param,
                            math::rng::RNG& rng)
{
  // Array
  fLatitude    = param.array_origin().latitude()/180.0*M_PI;
  fLongitude   = param.array_origin().longitude()/180.0*M_PI;
  fAltitude    = param.array_origin().elevation();

  if(param.array_layout_case() ==
     IsotropicDCArrayParameters::kHexArrayLayout)
  {
    const auto& layout = param.hex_array_layout();
    fSpacing     = layout.scope_spacing();
    fArrayParity = layout.scope_labeling_parity();

  // Telescopes
  unsigned num_telescopes = 
    3*param.ArrayNumTelRings*(param.ArrayNumTelRings+1)+1;
  std::set<unsigned> scopes_missing;
  Physics::tokenize(param.ScopeMissingList,scopes_missing);

  }
  else if(param.array_layout_case() ==
          IsotropicDCArrayParameters::kPrescribedArrayLayout)
  {

  }
  
  
  std::vector<double> pos_north, pos_east, pos_asl;
  VSDatumConverter< std::vector<double> >::
    fromString(pos_north,param.ScopePosNorth.c_str());
  VSDatumConverter< std::vector<double> >::
    fromString(pos_east,param.ScopePosEast.c_str());
  VSDatumConverter< std::vector<double> >::
    fromString(pos_asl,param.ScopePosASL.c_str());

  // Mirrors
  unsigned num_hex_mirror_rings = 
    unsigned(floor(param.ReflectorApert/2.0/param.MirrorFacetSpacing/cos(30.0/180.0*M_PI)))+2;

  // Camera
  Vec3D CameraFPTrans(param.CameraFPTransX, param.CameraFPTransY,
		      param.CameraFPTransZ);

  double FoV = 
    2.0*atan(param.CameraDiameter/2.0/CameraFPTrans.Norm())*180.0/M_PI;

  unsigned id = 0;
  for(unsigned i=0; i<num_telescopes; i++)
    {
      int hexid(i+1);

      if(scopes_missing.find(hexid) != scopes_missing.end())
	{
	  fTelescopesByHexID.push_back(0);
	  continue;
	}

      // Position
      Vec3D pos(0,0,param.ArrayAltitude+param.ReflectorIP);

      if(pos_north.size() > i && pos_east.size() > i && pos_asl.size() > i)
	{
	  pos.x = pos_east[i]  + rng.Normal() * param.ScopePosXYDisp;
	  pos.y = pos_north[i] + rng.Normal() * param.ScopePosXYDisp;
	  pos.z += pos_asl[i]  + rng.Normal() * param.ScopePosZDisp;
	}
      else
	{
	  nh_to_xy(&hexid,&pos.x,&pos.y);
	  if(fArrayParity)pos.x = -pos.x;
	  pos.x = pos.x * fSpacing + rng.Normal() * param.ScopePosXYDisp;
	  pos.y = pos.y * fSpacing + rng.Normal() * param.ScopePosXYDisp;
	  pos.z += rng.Normal() * param.ScopePosZDisp;
	}

      std::vector<VSOObscuration*> obsvec;      
      obsvec = 
	VSOObscuration::createObsVecFromString(param.ReflectorObscuration);

      VSOTelescope* telescope =
	new VSOTelescope(id, hexid, pos,
			 param.ScopeDeltaY, 
			 param.ScopeAlphaX, param.ScopeAlphaY,
			 param.ScopeElevation, param.ScopeAzimuth,
			 Vec3D(param.ScopeTranslationX, 
			       param.ScopeTranslationY, 
			       param.ScopeTranslationZ),
			 param.ReflectorCurvR, param.ReflectorApert,
			 param.MirrorFacetSpacing, param.MirrorFacetSize, 
			 param.ReflectorRot,
			 num_hex_mirror_rings,
			 param.ReflectorIP, param.MirrorLabelingParity,
			 CameraFPTrans, param.CameraDiameter, FoV, 
			 param.PixelDiameter, param.PixelSpacing, 
			 param.PixelConcSurvProb,
			 Vec3D(param.CameraFPRotX,param.CameraFPRotY,
			       param.CameraFPRotZ), param.CameraIP, 
			 param.PixelLabelingParity,
#if 0
			 param.SecondaryPresent,
			 param.SecondaryRefInd, param.SecondaryRefInd1,
			 param.SecondaryRefInd2, param.SecondaryCE10,
			 param.SecondaryCE12, param.SecondaryCE13,
			 param.SecondaryCE14, param.SecondaryCE15,
			 param.SecondaryCE20, param.SecondaryCE22,
			 param.SecondaryCE23, param.SecondaryCE24,
			 param.SecondaryCE25,
#endif			     
			 obsvec
			 );
      
      telescope->populateMirrorsAndPixelsRandom(param,rng);
      
      fTelescopes.push_back(telescope);      
      fTelescopesByHexID.push_back(telescope);
      id++;
    }

  return;
}
#endif

#if 0
void VSOArray::dumpShort(const std::string& filename) const
{
  std::ofstream stream(filename.c_str());
  if(stream.good())dumpShort(stream);
}

void VSOArray::dumpShort(std::ostream& stream) const
{
  stream 
    << "ARRAY "
    << VSDataConverter::toString(fTelescopes.size()) << ' '
    << VSDataConverter::toString(fSpacing) << ' '
    << VSDataConverter::toString(fLatitude) << ' '
    << VSDataConverter::toString(fLongitude) << ' '
    << VSDataConverter::toString(fAltitude) << ' '
    << VSDataConverter::toString(fArrayParity) << std::endl;
  
  for(std::vector<VSOTelescope*> ::const_iterator i = fTelescopes.begin();
      i!=fTelescopes.end(); i++)
    (*i)->dumpShort(stream);
}

void VSOArray::dump(std::ostream& stream, unsigned l) const
{
  stream << FDAV("Num telescopes", fTelescopes.size(), "", 30, l) << std::endl
	 << FDAV("Telescope spacing", fSpacing, "cm", 30, l) << std::endl
	 << FDAV("Latitude", fLatitude, "rad", 30, l) << std::endl
	 << FDAV("Longitude", fLongitude, "rad", 30, l) << std::endl
	 << FDAV("Altitude", fAltitude, "cm", 30, l) << std::endl
	 << FDAV("Array Parity", fArrayParity, "", 30, l) << std::endl;

  for(std::vector<VSOTelescope*> ::const_iterator i = fTelescopes.begin();
      i!=fTelescopes.end(); i++)
    {
      stream << std::endl;
      (*i)->dump(stream,l+1);
    }
}

bool VSOArray::readFromShortDump(const std::string& filename)
{
  std::ifstream stream(filename.c_str());
  if(stream.good())return readFromShortDump(stream);
  else return false;
}

bool VSOArray::readFromShortDump(std::istream& stream)
{
  std::string line;
  std::getline(stream,line);
  if(line.empty())return false;

  std::istringstream linestream(line);

  std::string keyword;
  linestream >> keyword;
  if(keyword!=std::string("ARRAY"))return false;

  unsigned telescopes_size;

  linestream
    >> telescopes_size
    >> fSpacing
    >> fLatitude
    >> fLongitude
    >> fAltitude
    >> fArrayParity;
  
  for(unsigned i=0; i<telescopes_size; i++)
    {
      VSOTelescope* telescope = VSOTelescope::createFromShortDump(stream);
      if(telescope==0)
	{
	  for(std::vector<VSOTelescope*>::iterator itel = fTelescopes.begin();
	      itel!=fTelescopes.end(); itel++)delete *itel;
	  delete telescope;
	  return false;
	}
      
      if(telescope->id() >= fTelescopes.size())
	fTelescopes.resize(telescope->id()+1);
      fTelescopes[telescope->id()]=telescope;
      
      if(telescope->hexID() > fTelescopesByHexID.size())
	fTelescopesByHexID.resize(telescope->hexID());
      fTelescopesByHexID[telescope->hexID()-1]=telescope;
    }

  return true;
}
#endif

#ifdef TEST_MAIN
  
#include <fstream>

int main(int argc, char** argv)
{
  RandomNumbers rng("random.seeds");

  ArrayParameters param;
  param.readFromArrayINIFile("array.ini");
  param.writeToArrayINIFile("array1.ini");

  VSOArray array;
  array.generateFromArrayParameters(param, rng);

  std::ofstream f1("array1.txt");
  array.dump(f1);

  Database db("array",true,false,false);
  ArrayParameters::createSimulationParametersTable(&db);
  VSOArray::createArrayTable(&db);
  VSOTelescope::createTelescopeTable(&db);
  VSOMirror::createMirrorTable(&db);
  VSOPixel::createPixelTable(&db);

  param.writeToDatabase(&db);
  array.writeToDatabase(&db);

  ArrayParameters param2;
  param2.readFromDatabase(&db);  
  param2.writeToArrayINIFile("array2.ini");

  VSOArray array2;
  array2.readFromDatabase(&db);

  std::ofstream f2("array2.txt");
  array2.dump(f2);
}

#endif



#ifdef TEST_MAIN_2
  
#include <fstream>

int main(int argc, char** argv)
{
  RandomNumbers rng("random.seeds");

  ArrayParameters param;
  param.readFromArrayINIFile("array.ini");

  VSOArray array;
  array.generateFromArrayParameters(param, rng);

  for(unsigned t=0; t< array.numTelescopes(); t++)
    for(unsigned m=0; m<array.telescope(t)->numMirrors(); m++)
      {
	Vec3D a(array.telescope(t)->mirror(m)->pos()+
		array.telescope(t)->mirror(m)->align());
	std::cout << a << ' ';
	array.telescope(t)->mirror(m)->reflectorToMirror(a);
	std::cout << a << std::endl;
      }
}

#endif
