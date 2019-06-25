/*

   calin/simulation/world_magnetic_model.hpp -- Stephen Fegan -- 2016-11-22

   Interface to the NOAA World Magnetic Model (WMM)

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#pragma once

#include <string>
#include <Eigen/Dense>

namespace calin { namespace simulation { namespace world_magnetic_model {

namespace internal {
struct WMM_Internals;
}

#ifdef SWIG
//%apply Eigen::Vector3d &OUTPUT { Eigen::Vector3d& field_nT };
#endif

class FieldVsElevation
{
public:
  virtual ~FieldVsElevation();
  virtual Eigen::Vector3d field_nT(double ellipsoid_elevation_cm) = 0;
};

class WMM;

class WMM_FieldVsElevation: public FieldVsElevation
{
public:
  WMM_FieldVsElevation(WMM* wmm, double latitude_deg,
    double longitude_deg_pos_east);
  virtual ~WMM_FieldVsElevation();
  Eigen::Vector3d field_nT(double ellipsoid_elevation_cm) override;
private:
  WMM* wmm_;
  double latitude_deg_;
  double longitude_deg_pos_east_;
};

class WMM
{
public:
  WMM(const std::string& cof_file = default_cof_file(), double date = 0);
  ~WMM();
  Eigen::Vector3d field_nT(double latitude_deg, double longitude_deg_pos_east,
    double ellipsoid_elevation_cm);
  static std::string default_cof_file();
  WMM_FieldVsElevation field_vs_elevation(double latitude_deg,
    double longitude_deg_pos_east);
protected:
  internal::WMM_Internals* mag_ = nullptr;
};

} } } // namespace calin::simulation::world_magnetic_model
