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
  double primary_reflection_cosine;

  Eigen::Vector3d secondary_position;
  int secondary_facet;
  double secondary_reflection_cosine;

  int obscuration_id;

  Eigen::Vector3d fp_position;

  Eigen::Vector3d camera_position;
  double camera_time;
  int camera_module_id;
  int camera_pixel_id;

  Eigen::Vector3d final_position;
};

class SCTRayTracer
{
public:
  SCTRayTracer(const calin::ix::simulation::sct_optics::SCTArray* array,
    calin::math::rng::RNG* rng = nullptr, bool adopt_array = false, bool adopt_rng = false);
  ~SCTRayTracer();

  unsigned num_telescopes() const { return scopes_.size(); }
  const Eigen::Vector3d& detector_sphere_center(unsigned iscope) {
    return scopes_.at(iscope)->detector_sphere_center; }
  double detector_sphere_radius(unsigned iscope) {
    return scopes_.at(iscope)->p_r_max; }

  const Eigen::Vector3d& reflector_position(unsigned iscope) const {
    return scopes_.at(iscope)->reflector_position; }
  const Eigen::Matrix3d& reflector_rotation(unsigned iscope) const {
    return scopes_.at(iscope)->reflector_rotation; }

  void point_telescope(unsigned iscope, double el_deg, double az_deg, double phi_deg = 0);
  void point_all_telescopes(double el_deg, double az_deg, double phi_deg = 0);

  bool trace_ray_in_global_frame(unsigned iscope, calin::math::ray::Ray& ray,
    SCTRayTracerResults& results, bool translate_back_to_global_frame = true) const;
  bool trace_ray_in_reflector_frame(unsigned iscope, calin::math::ray::Ray& ray,
    SCTRayTracerResults& results, bool skip_primary = false) const;


private:
  struct Facet
  {
    bool removed;
    double roughness;
    bool has_frame_change;
    Eigen::Vector3d offset;
    Eigen::Matrix3d rotation;
  };

  struct Telescope
  {
    ~Telescope() {
      delete p_scheme;
      delete p_scheme_loose;
      delete s_scheme;
      delete s_scheme_loose;
      for(auto* obs : primary_obscuration)delete obs;
      for(auto* obs : secondary_obscuration)delete obs;
      for(auto* obs : camera_obscuration)delete obs;
    }

    const calin::ix::simulation::sct_optics::SCTTelescope* param;

    SCTPrimaryFacetScheme* p_scheme;
    SCTPrimaryFacetScheme* p_scheme_loose;
    bool p_has_frame_change;
    bool p_facets_have_frame_change;
    Eigen::Vector3d p_offset;
    Eigen::Matrix3d p_rotation;
    double p_rho_min;
    double p_rho_max;
    const double* p_surface;
    unsigned p_surface_n;
    std::vector<Facet> p_facets;
    double p_r_max;
    double p_y_max;

    SCTSecondaryFacetScheme* s_scheme;
    SCTSecondaryFacetScheme* s_scheme_loose;
    bool s_has_frame_change;
    bool s_facets_have_frame_change;
    Eigen::Vector3d s_offset;
    Eigen::Matrix3d s_rotation;
    double s_rho_min;
    double s_rho_max;
    const double* s_surface;
    unsigned s_surface_n;
    std::vector<Facet> s_facets;

    bool c_has_frame_change;
    Eigen::Vector3d c_offset;
    Eigen::Matrix3d c_rotation;
    double c_rho_max;
    const double* c_surface;
    unsigned c_surface_n;

    std::vector<calin::simulation::vs_optics::VSOObscuration*> primary_obscuration;
    std::vector<calin::simulation::vs_optics::VSOObscuration*> secondary_obscuration;
    std::vector<calin::simulation::vs_optics::VSOObscuration*> camera_obscuration;

    Eigen::Vector3d reflector_position;
    Eigen::Matrix3d reflector_rotation;

    Eigen::Vector3d detector_sphere_center;
  };

  bool trace_ray_to_primary_in_reflector_frame(const Telescope* scope,
    calin::math::ray::Ray& ray, SCTRayTracerResults& results) const;
  bool trace_ray_to_secondary_in_reflector_frame(const Telescope* scope,
    calin::math::ray::Ray& ray, SCTRayTracerResults& results) const;

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
