/*

   calin/simulation/world_magnetic_model.cpp -- Stephen Fegan -- 2016-11-22

   Interface to the NOAA World Magnetic Model (WMM)

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <cstring>
#include <cerrno>

#include <provenance/system_info.hpp>
#include <simulation/GeomagnetismHeader.h>
#include <simulation/EGM9615.h>
#include <simulation/world_magnetic_model.hpp>

using namespace calin::simulation::world_magnetic_model;

FieldVsElevation::~FieldVsElevation()
{
  // nothing to see here
}

WMM_FieldVsElevation::WMM_FieldVsElevation(WMM* wmm, double latitude_deg,
    double longitude_deg_pos_east): FieldVsElevation(), wmm_(wmm),
  latitude_deg_(latitude_deg), longitude_deg_pos_east_(longitude_deg_pos_east)
{
  // nothing to see here
}

WMM_FieldVsElevation::~WMM_FieldVsElevation()
{
  // nothing to see here
}

Eigen::Vector3d WMM_FieldVsElevation::field_nT(double ellipsoid_elevation_cm)
{
  return wmm_->field_nT(latitude_deg_, longitude_deg_pos_east_,
    ellipsoid_elevation_cm);
}

namespace calin { namespace simulation { namespace world_magnetic_model { namespace internal {

struct WMM_Internals
{
  MAGtype_MagneticModel* magnetic_models[1] = { nullptr };
  double date;
  MAGtype_MagneticModel* timed_magnetic_model = nullptr;
  MAGtype_Ellipsoid ellipsoid;
  MAGtype_Geoid geoid;
};

} } } } // namespace calin::simulation::world_magnetic_model::internal


WMM::WMM(const std::string& cof_file, double date):
  mag_(new calin::simulation::world_magnetic_model::internal::WMM_Internals)
{
  if(!MAG_robustReadMagModels(cof_file.c_str(),mag_->magnetic_models, 1))
    throw std::runtime_error("Cannot open WMM COF file: " + cof_file
      + "\n" + strerror(errno));
  if(mag_->magnetic_models[0] == nullptr)
    throw std::runtime_error("Magnetic model is NULL");

  if(date == 0)mag_->date = mag_->magnetic_models[0]->epoch;
  else mag_->date = date;

  int nmax = mag_->magnetic_models[0]->nMax;
  int nterms = (nmax + 1) * (nmax + 2) / 2;
  mag_->timed_magnetic_model = MAG_AllocateModelMemory(nterms); /* For storing the time modified WMM Model parameters */
  if(mag_->timed_magnetic_model == nullptr)
    throw std::runtime_error("Could not allocate timed magnetic model");

  MAGtype_Date mag_date { 0, 0, 0, mag_->date };
  MAG_TimelyModifyMagneticModel(mag_date, mag_->magnetic_models[0],
    mag_->timed_magnetic_model); /* Time adjust the coefficients, Equation 19, WMM Technical report */

  MAG_SetDefaults(&mag_->ellipsoid, &mag_->geoid); /* Set default values and constants */
  /* Check for Geographic Poles */

  /* Set EGM96 Geoid parameters */
  mag_->geoid.GeoidHeightBuffer = GeoidHeightBuffer;
  mag_->geoid.Geoid_Initialized = 1;
  /* Set EGM96 Geoid parameters END */
}

WMM::~WMM()
{
  MAG_FreeMagneticModelMemory(mag_->timed_magnetic_model);
  MAG_FreeMagneticModelMemory(mag_->magnetic_models[0]);
}

Eigen::Vector3d WMM::field_nT(double latitude_deg, double longitude_deg_pos_east,
  double ellipsoid_elevation_cm)
{
  MAGtype_CoordGeodetic coord_geodetic;
  coord_geodetic.lambda               = longitude_deg_pos_east;
  coord_geodetic.phi                  = latitude_deg;
  coord_geodetic.HeightAboveEllipsoid = ellipsoid_elevation_cm*1.0e-5;
  coord_geodetic.HeightAboveGeoid     = 0.0;
  coord_geodetic.UseGeoid             = 0;

  MAGtype_CoordSpherical coord_spherical;
  MAG_GeodeticToSpherical(mag_->ellipsoid, coord_geodetic, &coord_spherical); /*Convert from geodetic to Spherical Equations: 17-18, WMM Technical report*/

  MAGtype_GeoMagneticElements geo_magnetic_elements;
  MAG_Geomag(mag_->ellipsoid, coord_spherical, coord_geodetic,
    mag_->timed_magnetic_model, &geo_magnetic_elements); /* Computes the geoMagnetic field elements and their time change*/

  return {
    geo_magnetic_elements.Y, /*6. Eastern component of the magnetic field vector*/
    geo_magnetic_elements.X, /*5. Northern component of the magnetic field vector*/
    -geo_magnetic_elements.Z /*7. Downward component of the magnetic field vector*/
  };
}

WMM_FieldVsElevation WMM::field_vs_elevation(double latitude_deg,
  double longitude_deg_pos_east)
{
  return WMM_FieldVsElevation(this, latitude_deg, longitude_deg_pos_east);
}

std::string WMM::default_cof_file()
{
  return calin::provenance::system_info::build_info()->data_install_dir()
    + "/simulation/WMM.COF";
}
