/*

   calin/simulation/sct_ray_tracer.hpp -- Stephen Fegan -- 2021-04-27

   Class for SCT ray tracing

   Copyright 2021, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <cmath>
#include <Eigen/Dense>

#include <math/rng.hpp>
#include <math/ray.hpp>
#include <simulation/vso_obscuration.hpp>
#include <simulation/sct_optics.pb.h>
#include <simulation/sct_facet_scheme.hpp>

namespace calin { namespace simulation { namespace sct_optics {

enum SCTRayTracerStatus
{
  RTS_MISSED_TELESCOPE_SPHERE   = 0,
  RTS_OBSCURED_BEFORE_PRIMARY   = 1,
  RTS_MISSED_PRIMARY            = 2,
  RTS_NO_PRIMARY_FACET          = 3,

  RTS_OBSCURED_BEFORE_SECONDARY = 4,
  RTS_MISSED_SECONDARY          = 5,
  RTS_NO_SECONDARY_FACET        = 6,

  RTS_OBSCURED_BEFORE_CAMERA    = 7,
  RTS_MISSED_CAMERA             = 8,

  RTS_COMPLETE                  = 9
};

struct SCTRayTracerResults
{
  SCTRayTracerStatus status;
  Eigen::Vector3d primary_position;
  int primary_facet;

  Eigen::Vector3d secondary_position;
  int secondary_facet;

  Eigen::Vector3d camera_position;
  Eigen::Vector3d final_position;
};

class SCTRayTracer
{
public:
  SCTRayTracer(const calin::ix::simulation::sct_optics::SCTArray* array,
    calin::math::rng::RNG* rng = nullptr, bool adopt_array = false, bool adopt_rng = false);
  ~SCTRayTracer();
  bool trace_ray_in_reflector_frame(unsigned iscope, calin::math::ray::Ray& ray,
    SCTRayTracerResults& results);
private:
  struct Telescope
  {
    ~Telescope() {
      delete p_scheme;
      delete s_scheme;
      for(auto* obs : primary_obscuration)delete obs;
      for(auto* obs : secondary_obscuration)delete obs;
      for(auto* obs : camera_obscuration)delete obs;
    }
    const calin::ix::simulation::sct_optics::SCTTelescope* param;

    SCTPrimaryFacetScheme* p_scheme;
    bool p_has_frame_change;
    Eigen::Vector3d p_offset;
    Eigen::Matrix3d p_rotation;
    double p_rho_min;
    double p_rho_max;
    const double* p_surface;
    unsigned p_surface_n;

    SCTSecondaryFacetScheme* s_scheme;
    bool s_has_frame_change;
    Eigen::Vector3d s_offset;
    Eigen::Matrix3d s_rotation;
    double s_rho_min;
    double s_rho_max;
    const double* s_surface;
    unsigned s_surface_n;

    bool c_has_frame_change;
    Eigen::Vector3d c_offset;
    Eigen::Matrix3d c_rotation;
    double c_rho_max;
    const double* c_surface;
    unsigned c_surface_n;

    std::vector<calin::simulation::vs_optics::VSOObscuration*> primary_obscuration;
    std::vector<calin::simulation::vs_optics::VSOObscuration*> secondary_obscuration;
    std::vector<calin::simulation::vs_optics::VSOObscuration*> camera_obscuration;
  };

  const calin::ix::simulation::sct_optics::SCTArray* array_ = nullptr;
  calin::math::rng::RNG* rng_ = nullptr;
  bool adopt_array_ = false;
  bool adopt_rng_ = false;
  std::vector<Telescope*> scopes_;
  double n_ = 1.0;
};

#ifndef SWIG
calin::ix::simulation::sct_optics::SCTTelescope* make_sct_telescope(
  calin::ix::simulation::sct_optics::SCTRandomArrayParameters& param,
  double x, double y, double z,
  calin::math::rng::RNG* rng = nullptr,
  calin::ix::simulation::sct_optics::SCTTelescope* telescope = nullptr);
calin::ix::simulation::sct_optics::SCTArray* make_sct_array(
  calin::ix::simulation::sct_optics::SCTRandomArrayParameters& param,
  calin::math::rng::RNG* rng = nullptr,
  calin::ix::simulation::sct_optics::SCTArray* array = nullptr);
#else
calin::ix::simulation::sct_optics::SCTTelescope* make_sct_telescope(
  calin::ix::simulation::sct_optics::SCTRandomArrayParameters& param,
  double x, double y, double z, calin::math::rng::RNG* rng = nullptr);
void make_sct_telescope(
  calin::ix::simulation::sct_optics::SCTRandomArrayParameters& param,
  double x, double y, double z, calin::math::rng::RNG* rng,
  calin::ix::simulation::sct_optics::SCTTelescope* telescope);
calin::ix::simulation::sct_optics::SCTArray* make_sct_array(
  calin::ix::simulation::sct_optics::SCTRandomArrayParameters& param,
  calin::math::rng::RNG* rng = nullptr);
void make_sct_array(
  calin::ix::simulation::sct_optics::SCTRandomArrayParameters& param,
  calin::math::rng::RNG* rng,
  calin::ix::simulation::sct_optics::SCTArray* array);
#endif

} } } // namespace calin::simulations::sct_optics
