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
#include <math/vector3d_util.hpp>
#include <simulation/vso_array.hpp>

using namespace calin::simulation::vs_optics;
using namespace calin::ix::simulation::vs_optics;

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

bool VSOArray::pointTelescopes(const Eigen::Vector3d& v)
{
  if(v.squaredNorm() == 0)
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

  std::vector<Eigen::Vector3d> scope_pos;

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
        Eigen::Vector3d pos;
        math::hex_array::hexid_to_xy(hexid, pos.x(), pos.y());
        if(array_parity)pos.x() = -pos.x();
        pos.x()  = pos.x() * spacing +
                   rng.normal() * layout.scope_position_dispersion_xy();
        pos.y()  = pos.y() * spacing +
                   rng.normal() * layout.scope_position_dispersion_xy();
        pos.z() += rng.normal() * layout.scope_position_dispersion_z();
        scope_pos.push_back(pos);
      }
  }
  else if(param.array_layout_case() ==
          IsotropicDCArrayParameters::kPrescribedArrayLayout)
  {
    for(auto pos : param.prescribed_array_layout().scope_positions())
      scope_pos.emplace_back(calin::math::vector3d_util::from_proto(pos));
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
  Eigen::Vector3d camera_fp_trans(
    calin::math::vector3d_util::from_proto(param.focal_plane().translation()));

  double FoV =
      2.0*atan(param.focal_plane().camera_diameter()/
               (2.0*camera_fp_trans.norm()))*180.0/M_PI;

  for(unsigned i=0; i<scope_pos.size(); i++)
    {
      // Position
      Eigen::Vector3d pos(scope_pos[i]);

      std::vector<VSOObscuration*> obsvec;
      for(const auto& obs : param.obscurations())
        obsvec.push_back(VSOObscuration::create_from_proto(obs));

      VSOTelescope* telescope =
      	new VSOTelescope(fTelescopes.size(), pos,
      			 param.reflector_frame().delta_y()*M_PI/180.0,
             param.reflector_frame().alpha_x()*M_PI/180.0,
             param.reflector_frame().alpha_y()*M_PI/180.0,
             param.reflector_frame().altaz().altitude()*M_PI/180.0,
             param.reflector_frame().altaz().azimuth()*M_PI/180.0,
      			 calin::math::vector3d_util::from_proto(param.reflector_frame().translation()),
      			 param.reflector().curvature_radius(),
             param.reflector().aperture(),
             param.reflector().facet_spacing(),
             param.reflector().facet_size(),
             param.reflector_frame().optic_axis_rotation()*M_PI/180.0,
             param.reflector_frame().facet_grid_shift_x(),
             param.reflector_frame().facet_grid_shift_z(),
             num_hex_mirror_rings,
             0.0, Eigen::Vector3d::Zero(), /*param.reflector().reflector_ip(),*/
             param.reflector().facet_labeling_parity(),
	           camera_fp_trans,
             param.focal_plane().camera_diameter(),
             FoV,
             param.pixel().cone_inner_diameter(),
             param.pixel().spacing(),
             param.pixel().grid_rotation()*M_PI/180.0,
             param.pixel().cone_survival_prob(),
             calin::math::vector3d_util::from_scaled_proto(param.focal_plane().rotation(), M_PI/180.0),
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

calin::ix::iact_data::instrument_layout::ArrayLayout*
calin::simulation::vs_optics::dc_parameters_to_array_layout(
  const ix::simulation::vs_optics::IsotropicDCArrayParameters& param,
  calin::ix::iact_data::instrument_layout::ArrayLayout* d)
{

}
