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

#include <math/hex_array.hpp>
#include <simulation/vso_array.hpp>

using namespace calin::simulation::vs_optics;
using namespace calin::ix::simulation::vs_optics;
using namespace calin::math::vs_physics;

VSOArray::VSOArray():
    fLatitude(), fLongitude(), fAltitude(), /* fSpacing(), fArrayParity(), */
    fTelescopes() /* , fTelescopesByHexID() */
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

bool VSOArray::pointTelescopes(const Vec3D& v)
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

void VSOArray::
generateFromArrayParameters(const IsotropicDCArrayParameters& param,
                            math::rng::RNG& rng)
{
  // Array
  fLatitude    = param.array_origin().latitude()/180.0*M_PI;
  fLongitude   = param.array_origin().longitude()/180.0*M_PI;
  fAltitude    = param.array_origin().elevation();

  std::vector<Vec3D> scope_pos;

  if(param.array_layout_case() ==
     IsotropicDCArrayParameters::kHexArrayLayout)
  {
    const auto& layout = param.hex_array_layout();
    double spacing     = layout.scope_spacing();
    bool array_parity  = layout.scope_labeling_parity();

    unsigned num_telescopes =
        math::hex_array::ringid_to_nsites_contained(layout.num_scope_rings());
    std::set<unsigned> scopes_missing;
    for(auto hexid : layout.scope_missing_list())
      scopes_missing.insert(hexid);

    for(unsigned hexid=0; hexid<num_telescopes; hexid++)
      if(scopes_missing.find(hexid) == scopes_missing.end())
      {
        Vec3D pos;
        math::hex_array::hexid_to_xy(hexid, pos.x, pos.y);
        if(array_parity)pos.x = -pos.x;
        pos.x  = pos.x * spacing +
                 rng.normal() * layout.scope_position_dispersion_xy();
        pos.y  = pos.y * spacing +
                 rng.normal() * layout.scope_position_dispersion_xy();
        pos.z += rng.normal() * layout.scope_position_dispersion_z();
        scope_pos.push_back(pos);
      }
  }
  else if(param.array_layout_case() ==
          IsotropicDCArrayParameters::kPrescribedArrayLayout)
  {
    for(auto pos : param.prescribed_array_layout().scope_positions())
      scope_pos.push_back(pos);
  }
  else
  {
    assert(0);
  }

  // Mirrors
  unsigned num_hex_mirror_rings = param.reflector().facet_num_hex_rings();
  if(param.reflector().aperture() > 0)
  {
    unsigned aperture_num_hex_mirror_rings =
        std::floor(param.reflector().aperture()/
           (2.0*param.reflector().facet_spacing()*CALIN_HEX_ARRAY_SQRT3*0.5))+2;
    if(aperture_num_hex_mirror_rings<num_hex_mirror_rings)
      num_hex_mirror_rings = aperture_num_hex_mirror_rings;
  }

  // Camera
  Vec3D camera_fp_trans(param.focal_plane().translation());

  double FoV =
      2.0*atan(param.focal_plane().camera_diameter()/
               (2.0*camera_fp_trans.Norm()))*180.0/M_PI;

  for(unsigned i=0; i<scope_pos.size(); i++)
    {
      // Position
      Vec3D pos(scope_pos[i]);

      std::vector<VSOObscuration*> obsvec;
      for(const auto& obs : param.obscurations())
        obsvec.push_back(VSOObscuration::create_from_proto(obs));

      VSOTelescope* telescope =
	new VSOTelescope(i, pos,
			 param.reflector_frame().delta_y()*M_PI/180.0,
                         param.reflector_frame().alpha_x()*M_PI/180.0,
                         param.reflector_frame().alpha_y()*M_PI/180.0,
                         param.reflector_frame().altaz().altitude()*M_PI/180.0,
                         param.reflector_frame().altaz().azimuth()*M_PI/180.0,
			 param.reflector_frame().translation(),
			 param.reflector().curvature_radius(),
                         param.reflector().aperture(),
			 param.reflector().facet_spacing(),
                         param.reflector().facet_size(),
                        param.reflector_frame().optic_axis_rotation()*M_PI/180.0,
			 num_hex_mirror_rings,
			 0.0, /*param.reflector().reflector_ip(),*/
                         param.reflector().facet_labeling_parity(),
			 camera_fp_trans,
                         param.focal_plane().camera_diameter(),
                         FoV,
			 param.pixel().cone_inner_diameter(),
                         param.pixel().spacing(),
			 param.pixel().cone_survival_prob(),
			 Vec3D(param.focal_plane().rotation(), M_PI/180.0),
                         0.0,
			 param.pixel().pixel_labeling_parity(),
			 obsvec
			 );

      telescope->populateMirrorsAndPixelsRandom(param,rng);

      fTelescopes.push_back(telescope);
    }

  return;
}

calin::ix::simulation::vs_optics::VSOArrayData* VSOArray::
dump_as_proto(calin::ix::simulation::vs_optics::VSOArrayData* d) const
{
  if(d == nullptr)d = new calin::ix::simulation::vs_optics::VSOArrayData;
  d->mutable_array_origin()->set_latitude(fLatitude*180.0/M_PI);
  d->mutable_array_origin()->set_longitude(fLongitude*180.0/M_PI);
  d->mutable_array_origin()->set_elevation(fAltitude);
  for(const auto* scope : fTelescopes)
    scope->dump_as_proto(d->add_telescope());
  return d;
}

VSOArray* VSOArray::
create_from_proto(const ix::simulation::vs_optics::VSOArrayData& d)
{
  VSOArray* array = new VSOArray;
  array->fLatitude  = d.array_origin().latitude()*M_PI/180.0;
  array->fLongitude = d.array_origin().longitude()*M_PI/180.0;
  array->fAltitude  = d.array_origin().elevation();
  for(const auto& scope : d.telescope())
    array->fTelescopes.push_back(VSOTelescope::create_from_proto(scope));
  return array;
}

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
