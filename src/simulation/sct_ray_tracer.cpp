/*

   calin/simulation/sct_ray_tracer.cpp -- Stephen Fegan -- 2021-04-27

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

#include <cmath>
#include <algorithm>
#include <stdexcept>

#include <util/log.hpp>
#include <math/least_squares.hpp>
#include <math/geometry.hpp>
#include <math/vector3d_util.hpp>
#include <math/special.hpp>
#include <simulation/sct_facet_scheme.hpp>
#include <simulation/sct_ray_tracer.hpp>

using namespace calin::simulation::sct_optics;
using namespace calin::util::log;
using calin::math::special::SQR;

SCTRayTracer::SCTRayTracer(const calin::ix::simulation::sct_optics::SCTArray* array,
    calin::math::rng::RNG* rng, bool adopt_array, bool adopt_rng):
  array_(array), rng_(rng ? rng : new calin::math::rng::RNG(__PRETTY_FUNCTION__,
    "SCT ray tracer")), adopt_array_(adopt_array), adopt_rng_(rng ? adopt_rng : true)
{
  for(const auto& scope_params : array->telescope())
  {
    auto* scope = new Telescope;
    scope->param = &scope_params;

    // *************************************************************************
    // PRIMARY
    // *************************************************************************

    scope->p_scheme = new SCTPrimaryFacetScheme(scope_params.primary_facet_scheme());
    scope->p_scheme_loose = nullptr;
    scope->p_surface = scope_params.primary_surface_polynomial().data();
    scope->p_surface_n = scope_params.primary_surface_polynomial_size();
    if(scope->p_surface_n == 0) {
      throw std::runtime_error("Primary surface polynomial must have at least one coefficient");
    }
    scope->p_rotation = calin::math::geometry::euler_to_matrix(scope_params.primary_rotation()).transpose();
    scope->p_offset = calin::math::vector3d_util::from_proto(scope_params.primary_offset());
    scope->p_has_frame_change = (scope->p_offset.squaredNorm()>0) or
      (not calin::math::geometry::euler_is_zero(scope_params.primary_rotation()));
    scope->p_facets_have_frame_change = false;
    scope->p_rho_max = SQR(scope->p_scheme->outer_radius());
    scope->p_rho_min = SQR(scope->p_scheme->inner_radius());
    scope->p_facets.resize(scope->p_scheme->num_facets());
    for(unsigned i=0;i<scope->p_scheme->num_facets();++i) {
      auto& facet = scope->p_facets[i];
      if(i<scope_params.primary_facets_size()) {
        Eigen::Vector3d nominal_position = scope->p_scheme->facet_centroid_3d(i,
          scope->p_surface,scope->p_surface_n);
        Eigen::Vector3d facet_normal = calin::math::geometry::norm_of_polynomial_surface(
          nominal_position.x(), nominal_position.z(), scope->p_surface, scope->p_surface_n);
        Eigen::Matrix3d facet_frame_matrix = calin::math::geometry::rotation_y_to_vec_Ryx(
          facet_normal);
        facet.removed = scope_params.primary_facets(i).removed();
        facet.roughness = 0.5 * scope_params.primary_facets(i).spot_size() * (M_PI/180.0);
        facet.rotation = facet_frame_matrix * calin::math::geometry::euler_to_matrix(
          scope_params.primary_facets(i).rotation()).transpose() * facet_frame_matrix.transpose();
        Eigen::Vector3d actual_position = nominal_position;
        if(scope_params.primary_facets(i).has_position()) {
          actual_position = calin::math::vector3d_util::from_proto(
            scope_params.primary_facets(i).position());
        }
        facet.offset = actual_position - facet.rotation.transpose()*nominal_position;
        facet.has_frame_change = (facet.offset.squaredNorm()>0) or
          (not calin::math::geometry::euler_is_zero(scope_params.primary_facets(i).rotation()));
        scope->p_facets_have_frame_change |= facet.has_frame_change;
      } else {
        facet.removed = true;
      }
    }
    if(scope->p_facets_have_frame_change and scope_params.has_primary_facet_scheme_loose()) {
      scope->p_scheme_loose = new SCTPrimaryFacetScheme(scope_params.primary_facet_scheme_loose());
      scope->p_rho_max = SQR(scope->p_scheme_loose->outer_radius());
      scope->p_rho_min = SQR(scope->p_scheme_loose->inner_radius());
      if(scope->p_scheme->num_facets() != scope->p_scheme_loose->num_facets()) {
        throw std::runtime_error("Primary surface facet schemes must have same number of facets.");
      }
    }
    scope->p_r_max = std::sqrt(scope->p_rho_max);
    scope->p_y_max = calin::math::least_squares::polyval(
      scope->p_surface, scope->p_surface_n, scope->p_rho_max);

    // *************************************************************************
    // SECONDARY
    // *************************************************************************

    scope->s_scheme = new SCTSecondaryFacetScheme(scope_params.secondary_facet_scheme());
    scope->s_scheme_loose = nullptr;

    scope->s_surface = scope_params.secondary_surface_polynomial().data();
    scope->s_surface_n = scope_params.secondary_surface_polynomial_size();
    if(scope->s_surface_n == 0) {
      throw std::runtime_error("Secondary surface polynomial must have at least one coefficient");
    }
    scope->s_rotation = calin::math::geometry::euler_to_matrix(scope_params.secondary_rotation()).transpose();
    scope->s_offset = Eigen::Vector3d(0, scope->s_surface[0], 0)
      + calin::math::vector3d_util::from_proto(scope_params.secondary_offset())
      - scope->s_rotation.transpose() * Eigen::Vector3d(0, scope->s_surface[0], 0);
    scope->s_has_frame_change = (scope->s_offset.squaredNorm()>0) or
      (not calin::math::geometry::euler_is_zero(scope_params.secondary_rotation()));
    scope->s_facets_have_frame_change = false;
    scope->s_rho_max = SQR(scope->s_scheme->outer_radius());
    scope->s_rho_min = SQR(scope->s_scheme->inner_radius());
    scope->s_facets.resize(scope->s_scheme->num_facets());
    for(unsigned i=0;i<scope->s_scheme->num_facets();++i) {
      auto& facet = scope->s_facets[i];
      if(i<scope_params.secondary_facets_size()) {
        Eigen::Vector3d nominal_position = scope->s_scheme->facet_centroid_3d(i,
          scope->s_surface,scope->s_surface_n);
        Eigen::Vector3d facet_normal = calin::math::geometry::norm_of_polynomial_surface(
          nominal_position.x(), nominal_position.z(), scope->s_surface, scope->s_surface_n);
        Eigen::Matrix3d facet_frame_matrix = calin::math::geometry::rotation_y_to_vec_Ryx(
          facet_normal);
        facet.removed = scope_params.secondary_facets(i).removed();
        facet.roughness = 0.5 * scope_params.secondary_facets(i).spot_size() * (M_PI/180.0);
        facet.rotation = facet_frame_matrix * calin::math::geometry::euler_to_matrix(
          scope_params.secondary_facets(i).rotation()).transpose() * facet_frame_matrix.transpose();
        Eigen::Vector3d actual_position = nominal_position;
        if(scope_params.secondary_facets(i).has_position()) {
          actual_position = calin::math::vector3d_util::from_proto(
            scope_params.secondary_facets(i).position());
        }
        facet.offset = actual_position - facet.rotation.transpose()*nominal_position;
        facet.has_frame_change = (facet.offset.squaredNorm()>0) or
          (not calin::math::geometry::euler_is_zero(scope_params.secondary_facets(i).rotation()));
        scope->s_facets_have_frame_change |= facet.has_frame_change;
      } else {
        facet.removed = true;
      }
    }

    if(scope->s_facets_have_frame_change and scope_params.has_secondary_facet_scheme_loose()) {
      scope->s_scheme_loose = new SCTSecondaryFacetScheme(scope_params.secondary_facet_scheme_loose());
      scope->s_rho_max = SQR(scope->s_scheme_loose->outer_radius());
      scope->s_rho_min = SQR(scope->s_scheme_loose->inner_radius());
      if(scope->s_scheme->num_facets() != scope->s_scheme_loose->num_facets()) {
        throw std::runtime_error("Secondary surface facet schemes must have same number of facets.");
      }
    }

    // *************************************************************************
    // CAMERA
    // *************************************************************************

    scope->c_surface = scope_params.camera_surface_polynomial().data();
    scope->c_surface_n = scope_params.camera_surface_polynomial_size();
    if(scope->c_surface_n == 0) {
      throw std::runtime_error("Camera surface polynomial must have at least one coefficient");
    }
    scope->c_rotation = calin::math::geometry::euler_to_matrix(scope_params.camera_rotation());
    scope->c_offset = Eigen::Vector3d(0, scope->c_surface[0], 0)
      + calin::math::vector3d_util::from_proto(scope_params.camera_offset())
      - scope->c_rotation.transpose() * Eigen::Vector3d(0, scope->c_surface[0], 0);
    scope->c_has_frame_change = (scope->c_offset.squaredNorm()>0) or
      (not calin::math::geometry::euler_is_zero(scope_params.camera_rotation()));
    scope->c_rho_max = SQR(scope_params.camera_radius());

    scope->c_module_pitch_inv = 1.0/scope_params.camera_module_grid().spacing();
    scope->c_module_nside = scope_params.camera_module_grid().num_side();
    scope->c_module_xc = scope_params.camera_module_grid().x_center();
    scope->c_module_yc = scope_params.camera_module_grid().y_center();
    scope->c_modules.resize(SQR(scope->c_module_nside), nullptr);
    for(const auto& m : scope_params.camera_modules()) {
      if(m.grid_id() >= scope->c_modules.size()) {
        throw std::runtime_error("Module grid_id out of range:\n"+m.DebugString());
      }
      if(scope->c_modules[m.grid_id()] != nullptr) {
        throw std::runtime_error("Duplicate module grid_id:\n"+m.DebugString());
      }
      auto* cm = scope->c_modules[m.grid_id()] = new CameraModule;
      cm->id = m.id();
      cm->first_pixel_id = m.first_pixel_id();
      calin::math::vector3d_util::set_from_proto(cm->position, m.center_position());
      cm->position.x() += scope_params.camera_module_pixel_grid().x_center();
      cm->position.z() += scope_params.camera_module_pixel_grid().y_center();
    }

    scope->c_pixel_pitch_inv = 1.0/scope_params.camera_module_pixel_grid().spacing();
    scope->c_pixel_nside = scope_params.camera_module_pixel_grid().num_side();
    scope->c_pixel_dead_space_fraction =
      scope_params.camera_module_pixel_grid().dead_space_width() * scope->c_pixel_pitch_inv;

    // *************************************************************************
    // OBSCURATIONS
    // *************************************************************************

    for(const auto& obs_param : scope_params.primary_obscuration()) {
      scope->primary_obscuration.push_back(
        calin::simulation::vs_optics::VSOObscuration::create_from_proto(obs_param));
    }
    for(const auto& obs_param : scope_params.secondary_obscuration()) {
      scope->secondary_obscuration.push_back(
        calin::simulation::vs_optics::VSOObscuration::create_from_proto(obs_param));
    }
    for(const auto& obs_param : scope_params.camera_obscuration()) {
      scope->camera_obscuration.push_back(
        calin::simulation::vs_optics::VSOObscuration::create_from_proto(obs_param));
    }

    scopes_.push_back(scope);

    point_telescope(scopes_.size()-1, 0, 0);
  }
}

SCTRayTracer::~SCTRayTracer()
{
  if(adopt_array_)delete array_;
  if(adopt_rng_)delete rng_;
  for(auto* scope : scopes_) {
    delete scope;
  }
}

bool SCTRayTracer::trace_ray_in_global_frame(unsigned iscope, calin::math::ray::Ray& ray,
  SCTRayTracerResults& results, bool translate_back_to_global_frame) const
{
  if(iscope >= scopes_.size()) {
    throw std::runtime_error("SCTRayTracer::trace_ray_in_global_frame: iscope out of range");
  }

  const Telescope* scope = scopes_[iscope];

  ray.translate_origin(scope->reflector_position);
  ray.rotate(scope->reflector_rotation);

  bool good = this->trace_ray_in_reflector_frame(iscope, ray, results);

  if(translate_back_to_global_frame) {
    ray.derotate(scope->reflector_rotation);
    ray.untranslate_origin(scope->reflector_position);
  }

  return good;
}

bool SCTRayTracer::trace_ray_to_primary_in_reflector_frame(const Telescope* scope,
  calin::math::ray::Ray& ray, SCTRayTracerResults& results) const
{
  // ***************************************************************************
  // Ray starts (and ends) in nominal telescope reflector frame
  // ***************************************************************************

  // Test for obscuration of incoming ray
  unsigned nobs = scope->primary_obscuration.size();
  unsigned obs_ihit  = nobs;
  double   obs_time  = ray.ct();
  for(unsigned iobs=0;iobs<nobs;iobs++)
  {
    math::ray::Ray ray_out;
    if(scope->primary_obscuration[iobs]->doesObscure(ray, ray_out, n_))
    {
      if((obs_ihit==nobs)||(ray_out.ct()<obs_time)) {
        obs_ihit = iobs;
        obs_time = ray_out.ct();
      }
    }
  }

  if(obs_ihit!=nobs)
  {
    results.status = RTS_OBSCURED_BEFORE_PRIMARY;
    results.obscuration_id = obs_ihit;
    return false;
  }

  // ***************************************************************************
  // TRANSFORM RAY INTO PRIMARY FRAME
  // ***************************************************************************

  if(scope->p_has_frame_change) {
    ray.translate_origin(scope->p_offset);
    ray.rotate(scope->p_rotation);
  }

  bool good = ray.propagate_to_polynomial_surface(scope->p_surface, scope->p_surface_n,
    scope->p_rho_min, scope->p_rho_max, /* time_reversal_ok= */ false, /* n= */ n_,
    /* ray_is_close_to_surface= */ false, /* tol= */ 1e-8);

  if(!good) {
    if(scope->p_has_frame_change) {
      ray.derotate(scope->p_rotation);
      ray.untranslate_origin(scope->p_offset);
    }
    results.status = RTS_MISSED_PRIMARY;
    return false;
  }

  if(scope->p_scheme_loose) {
    results.primary_facet = scope->p_scheme_loose->find_facet(ray.x(), ray.z());
  } else {
    results.primary_facet = scope->p_scheme->find_facet(ray.x(), ray.z());
  }
  const auto& primary_facet = scope->p_facets[results.primary_facet];
  if((results.primary_facet<0) or (primary_facet.removed))
  {
    if(scope->p_has_frame_change) {
      ray.derotate(scope->p_rotation);
      ray.untranslate_origin(scope->p_offset);
    }
    results.primary_position = ray.position();
    results.status = RTS_NO_PRIMARY_FACET;
    return false;
  }

  // ***************************************************************************
  // TRANSFORM INTO PRIMARY FACET FRAME IF THERE IS ONE
  // ***************************************************************************

  if(primary_facet.has_frame_change) {

    ray.translate_origin(primary_facet.offset);
    ray.rotate(primary_facet.rotation);

    good = ray.propagate_to_polynomial_surface(scope->p_surface, scope->p_surface_n,
      scope->p_rho_min, scope->p_rho_max, /* time_reversal_ok= */ true, /* n= */ n_,
      /* ray_is_close_to_surface= */ true, /* tol= */ 1e-8);
    if(not good) {
      ray.derotate(primary_facet.rotation);
      ray.untranslate_origin(primary_facet.offset);
      if(scope->p_has_frame_change) {
        ray.derotate(scope->p_rotation);
        ray.untranslate_origin(scope->p_offset);
      }
      results.status = RTS_MISSED_PRIMARY;
      return false;
    }
  }

  if(scope->p_scheme_loose) {
    if(scope->p_scheme->find_facet(ray.x(), ray.z()) != results.primary_facet) {
      if(primary_facet.has_frame_change) {
        ray.derotate(primary_facet.rotation);
        ray.untranslate_origin(primary_facet.offset);
      }
      if(scope->p_has_frame_change) {
        ray.derotate(scope->p_rotation);
        ray.untranslate_origin(scope->p_offset);
      }
      results.primary_position = ray.position();
      results.status = RTS_NO_PRIMARY_FACET;
      return false;
    }
  }

  if(primary_facet.roughness > 0) {
    results.primary_reflection_cosine =
      ray.reflect_from_rough_polynomial_surface(
        scope->p_surface, scope->p_surface_n, primary_facet.roughness, *rng_);
  } else {
    results.primary_reflection_cosine =
      ray.reflect_from_polynomial_surface(
        scope->p_surface, scope->p_surface_n);
  }

  // ***************************************************************************
  // TRANSFORM BACK INTO PRIMARY FRAME
  // ***************************************************************************

  if(primary_facet.has_frame_change) {
    ray.derotate(primary_facet.rotation);
    ray.untranslate_origin(primary_facet.offset);
  }

  // ***************************************************************************
  // TRANSFORM RAY BACK INTO REFLECTOR FRAME
  // ***************************************************************************

  if(scope->p_has_frame_change) {
    ray.derotate(scope->p_rotation);
    ray.untranslate_origin(scope->p_offset);
  }

  results.primary_position = ray.position();
  return true;
}

bool SCTRayTracer::trace_ray_to_secondary_in_reflector_frame(const Telescope* scope,
  calin::math::ray::Ray& ray, SCTRayTracerResults& results) const
{
  // ***************************************************************************
  // Ray starts (and ends) in nominal telescope reflector frame
  // ***************************************************************************

  // Test for obscuration of ray between primary and secondary
  unsigned nobs = scope->secondary_obscuration.size();
  unsigned obs_ihit  = nobs;
  double   obs_time  = ray.ct();
  for(unsigned iobs=0;iobs<nobs;iobs++)
  {
    math::ray::Ray ray_out;
    if(scope->secondary_obscuration[iobs]->doesObscure(ray, ray_out, n_))
    {
      if((obs_ihit==nobs)||(ray_out.ct()<obs_time)) {
        obs_ihit = iobs;
        obs_time = ray_out.ct();
      }
    }
  }

  if(obs_ihit!=nobs)
  {
    results.status = RTS_OBSCURED_BEFORE_SECONDARY;
    results.obscuration_id = obs_ihit + scope->primary_obscuration.size();
    return false;
  }

  // ***************************************************************************
  // TRANSFORM RAY INTO SECONDARY FRAME
  // ***************************************************************************

  if(scope->s_has_frame_change) {
    ray.translate_origin(scope->s_offset);
    ray.rotate(scope->s_rotation);
  }

  bool good = ray.propagate_to_polynomial_surface(scope->s_surface, scope->s_surface_n,
    scope->s_rho_min, scope->s_rho_max, /* time_reversal_ok= */ false, /* n= */ n_,
    /* ray_is_close_to_surface= */ false, /* tol= */ 1e-8);

  if(!good) {
    if(scope->s_has_frame_change) {
      ray.derotate(scope->s_rotation);
      ray.untranslate_origin(scope->s_offset);
    }
    results.status = RTS_MISSED_SECONDARY;
    return false;
  }

  if(scope->s_scheme_loose) {
    results.secondary_facet = scope->s_scheme_loose->find_facet(ray.x(), ray.z());
  } else {
    results.secondary_facet = scope->s_scheme->find_facet(ray.x(), ray.z());
  }
  const auto& secondary_facet = scope->s_facets[results.secondary_facet];
  if((results.secondary_facet<0) or (secondary_facet.removed))
  {
    if(scope->s_has_frame_change) {
      ray.derotate(scope->s_rotation);
      ray.untranslate_origin(scope->s_offset);
    }
    results.secondary_position = ray.position();
    results.status = RTS_NO_SECONDARY_FACET;
    return false;
  }

  // ***************************************************************************
  // TRANSFORM INTO SECONDARY FACET FRAME IF THERE IS ONE
  // ***************************************************************************

  if(secondary_facet.has_frame_change) {
    ray.translate_origin(secondary_facet.offset);
    ray.rotate(secondary_facet.rotation);

    good = ray.propagate_to_polynomial_surface(scope->s_surface, scope->s_surface_n,
      scope->s_rho_min, scope->s_rho_max, /* time_reversal_ok= */ true, /* n= */ n_,
      /* ray_is_close_to_surface= */ true, /* tol= */ 1e-8);
    if(not good) {
      ray.derotate(secondary_facet.rotation);
      ray.untranslate_origin(secondary_facet.offset);
      if(scope->s_has_frame_change) {
        ray.derotate(scope->s_rotation);
        ray.untranslate_origin(scope->s_offset);
      }
      results.status = RTS_MISSED_SECONDARY;
      return false;
    }
  }

  if(scope->s_scheme_loose) {
    if(scope->s_scheme->find_facet(ray.x(), ray.z()) != results.secondary_facet) {
      if(secondary_facet.has_frame_change) {
        ray.derotate(secondary_facet.rotation);
        ray.untranslate_origin(secondary_facet.offset);
      }
      if(scope->s_has_frame_change) {
        ray.derotate(scope->s_rotation);
        ray.untranslate_origin(scope->s_offset);
      }
      results.secondary_position = ray.position();
      results.status = RTS_NO_SECONDARY_FACET;
      return false;
    }
  }

  if(secondary_facet.roughness > 0) {
    results.secondary_reflection_cosine =
      ray.reflect_from_rough_polynomial_surface(
        scope->s_surface, scope->s_surface_n, secondary_facet.roughness, *rng_);
  } else {
    results.secondary_reflection_cosine =
      ray.reflect_from_polynomial_surface(
        scope->s_surface, scope->s_surface_n);
  }

  // ***************************************************************************
  // TRANSFORM BACK INTO SECONDARY FRAME
  // ***************************************************************************

  if(secondary_facet.has_frame_change) {
    ray.derotate(secondary_facet.rotation);
    ray.untranslate_origin(secondary_facet.offset);
  }

  // ***************************************************************************
  // TRANSFORM RAY BACK INTO REFLECTOR FRAME
  // ***************************************************************************

  if(scope->s_has_frame_change) {
    ray.derotate(scope->s_rotation);
    ray.untranslate_origin(scope->s_offset);
  }

  results.secondary_position = ray.position();
  return true;
}

bool SCTRayTracer::trace_ray_in_reflector_frame(unsigned iscope, calin::math::ray::Ray& ray,
  SCTRayTracerResults& results, bool skip_primary) const
{
  if(iscope >= scopes_.size()) {
    throw std::runtime_error("SCTRayTracer::trace_ray_in_reflector_frame: iscope out of range");
  }

  const Telescope* scope = scopes_[iscope];

  if((not skip_primary) and (not trace_ray_to_primary_in_reflector_frame(scope, ray, results))) {
    return false;
  }

  if(not trace_ray_to_secondary_in_reflector_frame(scope, ray, results)) {
    return false;
  }

  // ***************************************************************************
  // Ray starts (and ends) in nominal telescope reflector frame
  // ***************************************************************************

  // Test for obscuration of ray between secondary and camera
  unsigned nobs = scope->camera_obscuration.size();
  unsigned obs_ihit  = nobs;
  double   obs_time  = ray.ct();
  for(unsigned iobs=0;iobs<nobs;iobs++)
  {
    math::ray::Ray ray_out;
    if(scope->camera_obscuration[iobs]->doesObscure(ray, ray_out, n_))
    {
      if((obs_ihit==nobs)||(ray_out.ct()<obs_time)) {
        obs_ihit = iobs;
        obs_time = ray_out.ct();
      }
    }
  }

  if(obs_ihit!=nobs)
  {
    results.status = RTS_OBSCURED_BEFORE_CAMERA;
    results.obscuration_id = obs_ihit
      + scope->secondary_obscuration.size() + scope->primary_obscuration.size();
    return false;
  }

  // ***************************************************************************
  // TRANSFORM RAY INTO CAMERA FRAME
  // ***************************************************************************

  if(scope->c_has_frame_change) {
    ray.translate_origin(scope->c_offset);
    ray.rotate(scope->c_rotation);
  }

  bool good = ray.propagate_to_polynomial_surface(scope->c_surface, scope->c_surface_n,
    0, scope->c_rho_max, /* time_reversal_ok= */ false, /* n= */ n_,
    /* ray_is_close_to_surface= */ false, /* tol= */ 1e-8);

  if(!good) {
    if(scope->c_has_frame_change) {
      ray.derotate(scope->c_rotation);
      ray.untranslate_origin(scope->c_offset);
    }
    results.status = RTS_MISSED_CAMERA;
    return false;
  }

  results.fp_position = ray.position();

  // ***************************************************************************
  // Find camera module
  // ***************************************************************************

  results.camera_module_grid_id = calin::math::geometry::find_square_grid_site(
    ray.x(), ray.z(), scope->c_module_pitch_inv, scope->c_module_nside,
    scope->c_module_xc, scope->c_module_yc);

  if(results.camera_module_grid_id < 0) {
    if(scope->c_has_frame_change) {
      ray.derotate(scope->c_rotation);
      ray.untranslate_origin(scope->c_offset);
    }
    results.status = RTS_NO_MODULE;
    return false;
  }

  const auto* module = scope->c_modules[results.camera_module_grid_id];
  if(module == nullptr) {
    if(scope->c_has_frame_change) {
      ray.derotate(scope->c_rotation);
      ray.untranslate_origin(scope->c_offset);
    }
    results.status = RTS_NO_MODULE;
    return false;
  }

  results.camera_module_id = module->id;

  // ***************************************************************************
  // Find pixel
  // ***************************************************************************

  good = ray.propagate_to_y_plane(-module->position.y(), /* time_reversal_ok= */ true,
    scope->c_surface_n);
  if(!good) {
    if(scope->c_has_frame_change) {
      ray.derotate(scope->c_rotation);
      ray.untranslate_origin(scope->c_offset);
    }
    results.status = RTS_MISSED_MODULE;
    return false;
  }

  results.camera_module_pixel_id = calin::math::geometry::find_square_grid_site(
    ray.x(), ray.z(), scope->c_pixel_pitch_inv, scope->c_pixel_nside,
    module->position.x(), module->position.z(), scope->c_pixel_dead_space_fraction);

  if(results.camera_module_pixel_id < 0) {
    if(scope->c_has_frame_change) {
      ray.derotate(scope->c_rotation);
      ray.untranslate_origin(scope->c_offset);
    }
    results.status = RTS_NO_PIXEL;
    return false;
  }

  results.camera_pixel_id = module->first_pixel_id + results.camera_module_pixel_id;
  results.camera_position = ray.position();
  results.camera_time = ray.time();

  // ***************************************************************************
  // TRANSFORM RAY BACK INTO REFLECTOR FRAME
  // ***************************************************************************

  if(scope->c_has_frame_change) {
    ray.derotate(scope->c_rotation);
    ray.untranslate_origin(scope->c_offset);
  }

  results.final_position = ray.position();
  results.status = RTS_COMPLETE;
  return true;
}

void SCTRayTracer::point_telescope(unsigned iscope, double el_deg, double az_deg, double phi_deg)
{
  if(iscope >= scopes_.size()) {
    throw std::runtime_error("SCTRayTracer::point_telescope: iscope out of range");
  }

  Telescope* scope = scopes_[iscope];

  Eigen::Quaterniond rot;
  Eigen::Vector3d pos;

  pos = calin::math::vector3d_util::from_proto(scope->param->position());
  rot = Eigen::AngleAxisd(az_deg*M_PI/180.0, Eigen::Vector3d::UnitZ());

  pos += rot.inverse() * Eigen::Vector3d(0, scope->param->azimuth_elevation_axes_separation(), 0);
  rot = Eigen::AngleAxisd(-el_deg*M_PI/180.0, Eigen::Vector3d::UnitX()) * rot;

  pos += rot.inverse() * calin::math::vector3d_util::from_proto(scope->param->reflector_origin());
  rot = Eigen::AngleAxisd(phi_deg*M_PI/180.0, Eigen::Vector3d::UnitZ()) * rot;

  scope->reflector_position = pos;
  scope->reflector_rotation = rot;

  scope->detector_sphere_center = pos
    + scope->reflector_rotation.transpose()*Eigen::Vector3d(0,scope->p_y_max,0);
}

void SCTRayTracer::point_all_telescopes(double el_deg, double az_deg, double phi_deg)
{
  for(unsigned iscope=0; iscope<scopes_.size(); ++iscope) {
    point_telescope(iscope, el_deg, az_deg, phi_deg);
  }
}

void SCTRayTracer::pixel_centers(unsigned iscope, Eigen::VectorXd& x_out, Eigen::VectorXd& z_out) const
{
  if(iscope >= scopes_.size()) {
    throw std::runtime_error("SCTRayTracer::point_telescope: iscope out of range");
  }

  calin::simulation::sct_optics::pixel_centers(*scopes_[iscope]->param, x_out, z_out);
}


calin::ix::simulation::sct_optics::SCTTelescope*
calin::simulation::sct_optics::make_sct_telescope(
  calin::ix::simulation::sct_optics::SCTRandomArrayParameters& param,
  double telescope_x, double telescope_y, double telescope_z,
  calin::math::rng::RNG* rng,
  calin::ix::simulation::sct_optics::SCTTelescope* telescope)
{
  if(rng == nullptr) {
    rng = new calin::math::rng::RNG(__PRETTY_FUNCTION__, "Random SCT telescope generation");
  }
  if(telescope == nullptr) {
    telescope = new calin::ix::simulation::sct_optics::SCTTelescope();
  } else {
    telescope->Clear();
  }

  if(param.telescope_position_xy_dispersion() > 0) {
    telescope_x += param.telescope_position_xy_dispersion()*rng->normal();
    telescope_y += param.telescope_position_xy_dispersion()*rng->normal();
  }
  if(param.telescope_position_z_dispersion() > 0) {
    telescope_z += param.telescope_position_z_dispersion()*rng->normal();
  }

  telescope->mutable_position()->set_x(telescope_x);
  telescope->mutable_position()->set_y(telescope_y);
  telescope->mutable_position()->set_z(telescope_z);

  telescope->set_azimuth_elevation_axes_separation(
    param.azimuth_elevation_axes_separation());
  telescope->mutable_reflector_origin()->CopyFrom(
    param.reflector_origin());

  // ***************************************************************************
  // ***************************************************************************
  //
  // PRIMARY MIRROR SURFACE AND FACETS
  //
  // ***************************************************************************
  // ***************************************************************************

  for(unsigned i=0;i<param.primary_sag_polynomial_size();++i) {
    double p_i = param.primary_sag_polynomial(i);
    telescope->add_primary_surface_polynomial(p_i);
  }

  Eigen::Vector3d primary_offset =
    calin::math::vector3d_util::from_proto(param.primary_offset());
  if(param.primary_offset_xz_dispersion() > 0) {
    primary_offset.x() += param.primary_offset_xz_dispersion()*rng->normal();
    primary_offset.z() += param.primary_offset_xz_dispersion()*rng->normal();
  }
  if(param.primary_offset_y_dispersion() > 0) {
    primary_offset.y() += param.primary_offset_y_dispersion()*rng->normal();
  }
  if(primary_offset.squaredNorm() > 0) {
    calin::math::vector3d_util::dump_as_proto(primary_offset,
      telescope->mutable_primary_offset());
  }

  Eigen::Quaterniond primary_rotation =
    calin::math::geometry::euler_to_quaternion(param.primary_rotation());
  if(param.primary_rotation_dispersion() >= 0) {
    primary_rotation =
      calin::math::geometry::euler_to_quaternion(
        calin::math::geometry::scattering_euler(param.primary_rotation_dispersion(),
          *rng, calin::ix::common_types::EulerAngles3D::YXY)) * primary_rotation;
  }
  if(primary_rotation.vec().squaredNorm() > 0) {
    telescope->mutable_primary_rotation()->set_rotation_order(
      calin::ix::common_types::EulerAngles3D::YXY);
    calin::math::geometry::quaternion_to_euler(
      telescope->mutable_primary_rotation(), primary_rotation);
  }
  telescope->mutable_primary_facet_scheme()->CopyFrom(param.primary_facet_scheme());
  if(param.has_primary_facet_scheme_loose()) {
    telescope->mutable_primary_facet_scheme_loose()->CopyFrom(param.primary_facet_scheme_loose());
  }

  auto primary_facet_scheme = calin::simulation::sct_optics::SCTPrimaryFacetScheme(
    param.primary_facet_scheme());
  for(unsigned ifacet=0; ifacet<primary_facet_scheme.num_facets(); ++ifacet) {
    auto* facet = telescope->add_primary_facets();
    facet->set_id(ifacet);
    Eigen::Vector3d facet_position(0,0,0);
    primary_facet_scheme.facet_centroid(ifacet, facet_position.x(), facet_position.z());
    double facet_off_axis_dist = facet_position.norm();
    facet_position.y() = calin::math::least_squares::polyval(
      telescope->primary_surface_polynomial().data(),
      telescope->primary_surface_polynomial_size(), facet_position.squaredNorm());
    Eigen::Vector3d facet_normal = calin::math::geometry::norm_of_polynomial_surface(
      facet_position.x(), facet_position.z(), telescope->primary_surface_polynomial().data(),
      telescope->primary_surface_polynomial_size());
    Eigen::Matrix3d facet_frame_matrix = calin::math::geometry::rotation_y_to_vec_Ryx(
      facet_normal);
    Eigen::Vector3d facet_offset(0,0,0);
    if(param.primary_facet_offset_xz_dispersion() > 0) {
      facet_offset.x() = param.primary_facet_offset_xz_dispersion()*rng->normal();
      facet_offset.z() = param.primary_facet_offset_xz_dispersion()*rng->normal();
    }
    if(param.primary_facet_offset_y_dispersion() > 0) {
      facet_offset.y() = param.primary_facet_offset_y_dispersion()*rng->normal();
    }
    facet_position += facet_frame_matrix * facet_offset;
    calin::math::vector3d_util::dump_as_proto(facet_position, facet->mutable_position());

    Eigen::Quaterniond facet_rotation = Eigen::Quaterniond::Identity();
    if(param.primary_facet_rho_rotation_dispersion() > 0) {
      facet_rotation = Eigen::AngleAxisd(
        param.primary_facet_rho_rotation_dispersion()*rng->normal()*(M_PI/180.0),
        Eigen::Vector3d::UnitY()) * facet_rotation;
    }
    if(param.primary_facet_phi_rotation_dispersion() > 0) {
      facet_rotation = Eigen::AngleAxisd(
        param.primary_facet_phi_rotation_dispersion()*rng->normal()*(M_PI/180.0),
        Eigen::Vector3d::UnitZ()) * facet_rotation;
    }
    if(param.primary_facet_theta_rotation_dispersion() > 0) {
      facet_rotation = Eigen::AngleAxisd(
        (param.primary_facet_theta_rotation_dispersion()*rng->normal()
          + facet_off_axis_dist*param.primary_facet_theta_canting())*(M_PI/180.0),
        Eigen::Vector3d::UnitX()) * facet_rotation;
    } else if(param.primary_facet_theta_canting() != 0) {
      facet_rotation = Eigen::AngleAxisd(
        facet_off_axis_dist*param.primary_facet_theta_canting()*(M_PI/180.0),
        Eigen::Vector3d::UnitX()) * facet_rotation;
    }
    if(facet_rotation.vec().squaredNorm() > 0) {
      auto* facet_euler = facet->mutable_rotation();
      facet_euler->set_rotation_order(calin::ix::common_types::EulerAngles3D::YXY);
      calin::math::geometry::quaternion_to_euler(facet_euler, facet_rotation);
        // facet_frame_matrix * facet_rotation * facet_frame_matrix.transpose());
    }

    if(param.primary_facet_spot_size_mean()>0
        and param.primary_facet_spot_size_dispersion()>0) {
      facet->set_spot_size(rng->gamma_by_mean_and_sigma(
        param.primary_facet_spot_size_mean(), param.primary_facet_spot_size_dispersion()));
    } else if(param.primary_facet_spot_size_mean()>0) {
      facet->set_spot_size(param.primary_facet_spot_size_mean());
    }
  }
  for(auto id : param.primary_facets_removed()) {
    if(id < primary_facet_scheme.num_facets()) {
      telescope->mutable_primary_facets(id)->set_removed(true);
    }
  }

  // ***************************************************************************
  // ***************************************************************************
  //
  // SECONDARY MIRROR SURFACE AND FACETS
  //
  // ***************************************************************************
  // ***************************************************************************

  for(unsigned i=0;i<param.secondary_sag_polynomial_size();++i) {
    double p_i = param.secondary_sag_polynomial(i);
    if(i == 0) {
      p_i += param.secondary_distance();
    }
    telescope->add_secondary_surface_polynomial(p_i);
  }

  Eigen::Vector3d secondary_offset =
    calin::math::vector3d_util::from_proto(param.secondary_offset());
  if(param.secondary_offset_xz_dispersion() > 0) {
    secondary_offset.x() += param.secondary_offset_xz_dispersion()*rng->normal();
    secondary_offset.z() += param.secondary_offset_xz_dispersion()*rng->normal();
  }
  if(param.secondary_offset_y_dispersion() > 0) {
    secondary_offset.y() += param.secondary_offset_y_dispersion()*rng->normal();
  }
  if(secondary_offset.squaredNorm() > 0) {
    calin::math::vector3d_util::dump_as_proto(secondary_offset,
      telescope->mutable_secondary_offset());
  }

  Eigen::Quaterniond secondary_rotation =
    calin::math::geometry::euler_to_quaternion(param.secondary_rotation());
  if(param.secondary_rotation_dispersion() >= 0) {
    secondary_rotation =
      calin::math::geometry::euler_to_quaternion(
        calin::math::geometry::scattering_euler(param.secondary_rotation_dispersion(),
          *rng, calin::ix::common_types::EulerAngles3D::YXY)) * secondary_rotation;
  }
  if(secondary_rotation.vec().squaredNorm() > 0) {
    telescope->mutable_secondary_rotation()->set_rotation_order(
      calin::ix::common_types::EulerAngles3D::YXY);
    calin::math::geometry::quaternion_to_euler(
      telescope->mutable_secondary_rotation(), secondary_rotation);
  }
  telescope->mutable_secondary_facet_scheme()->CopyFrom(param.secondary_facet_scheme());
  if(param.has_secondary_facet_scheme_loose()) {
    telescope->mutable_secondary_facet_scheme_loose()->CopyFrom(param.secondary_facet_scheme_loose());
  }

  auto secondary_facet_scheme = calin::simulation::sct_optics::SCTSecondaryFacetScheme(
    param.secondary_facet_scheme());
  for(unsigned ifacet=0; ifacet<secondary_facet_scheme.num_facets(); ++ifacet) {
    auto* facet = telescope->add_secondary_facets();
    facet->set_id(ifacet);
    Eigen::Vector3d facet_position(0,0,0);
    secondary_facet_scheme.facet_centroid(ifacet, facet_position.x(), facet_position.z());
    double facet_off_axis_dist = facet_position.norm();
    facet_position.y() = calin::math::least_squares::polyval(
      telescope->secondary_surface_polynomial().data(),
      telescope->secondary_surface_polynomial_size(), facet_position.squaredNorm());
    Eigen::Vector3d facet_normal = calin::math::geometry::norm_of_polynomial_surface(
      facet_position.x(), facet_position.z(), telescope->secondary_surface_polynomial().data(),
      telescope->secondary_surface_polynomial_size());
    Eigen::Matrix3d facet_frame_matrix = calin::math::geometry::rotation_y_to_vec_Ryx(
      facet_normal);
    Eigen::Vector3d facet_offset(0,0,0);
    if(param.secondary_facet_offset_xz_dispersion() > 0) {
      facet_offset.x() = param.secondary_facet_offset_xz_dispersion()*rng->normal();
      facet_offset.z() = param.secondary_facet_offset_xz_dispersion()*rng->normal();
    }
    if(param.secondary_facet_offset_y_dispersion() > 0) {
      facet_offset.y() = param.secondary_facet_offset_y_dispersion()*rng->normal();
    }
    facet_position += facet_frame_matrix * facet_offset;
    calin::math::vector3d_util::dump_as_proto(facet_position, facet->mutable_position());

    Eigen::Quaterniond facet_rotation = Eigen::Quaterniond::Identity();
    if(param.secondary_facet_rho_rotation_dispersion() > 0) {
      facet_rotation = Eigen::AngleAxisd(
        param.secondary_facet_rho_rotation_dispersion()*rng->normal()*(M_PI/180.0),
        Eigen::Vector3d::UnitY()) * facet_rotation;
    }
    if(param.secondary_facet_phi_rotation_dispersion() > 0) {
      facet_rotation = Eigen::AngleAxisd(
        param.secondary_facet_phi_rotation_dispersion()*rng->normal()*(M_PI/180.0),
        Eigen::Vector3d::UnitZ()) * facet_rotation;
    }
    if(param.secondary_facet_theta_rotation_dispersion() > 0) {
      facet_rotation = Eigen::AngleAxisd(
        (param.secondary_facet_theta_rotation_dispersion()*rng->normal()
          + facet_off_axis_dist*param.secondary_facet_theta_canting())*(M_PI/180.0),
        Eigen::Vector3d::UnitX()) * facet_rotation;
    } else if(param.secondary_facet_theta_canting() != 0) {
      facet_rotation = Eigen::AngleAxisd(
        facet_off_axis_dist*param.secondary_facet_theta_canting()*(M_PI/180.0),
        Eigen::Vector3d::UnitX()) * facet_rotation;
    }
    if(facet_rotation.vec().squaredNorm() > 0) {
      auto* facet_euler = facet->mutable_rotation();
      facet_euler->set_rotation_order(calin::ix::common_types::EulerAngles3D::YXY);
      calin::math::geometry::quaternion_to_euler(facet_euler, facet_rotation);
        // facet_frame_matrix * facet_rotation * facet_frame_matrix.transpose());
    }

    if(param.secondary_facet_spot_size_mean()>0
        and param.secondary_facet_spot_size_dispersion()>0) {
      facet->set_spot_size(rng->gamma_by_mean_and_sigma(
        param.secondary_facet_spot_size_mean(), param.secondary_facet_spot_size_dispersion()));
    } else if(param.secondary_facet_spot_size_mean()>0) {
      facet->set_spot_size(param.secondary_facet_spot_size_mean());
    }
  }
  for(auto id : param.secondary_facets_removed()) {
    if(id < secondary_facet_scheme.num_facets()) {
      telescope->mutable_secondary_facets(id)->set_removed(true);
    }
  }

  // ***************************************************************************
  // ***************************************************************************
  //
  // CAMERA SURFACE AND MODULES
  //
  // ***************************************************************************
  // ***************************************************************************

  for(unsigned i=0;i<param.camera_sag_polynomial_size();++i) {
    double p_i = param.camera_sag_polynomial(i);
    if(i == 0) {
      p_i += param.camera_distance();
    }
    telescope->add_camera_surface_polynomial(p_i);
  }

  Eigen::Vector3d camera_offset =
    calin::math::vector3d_util::from_proto(param.camera_offset());
  if(param.camera_offset_xz_dispersion() > 0) {
    camera_offset.x() += param.camera_offset_xz_dispersion()*rng->normal();
    camera_offset.z() += param.camera_offset_xz_dispersion()*rng->normal();
  }
  if(param.camera_offset_y_dispersion() > 0) {
    camera_offset.y() += param.camera_offset_y_dispersion()*rng->normal();
  }
  if(camera_offset.squaredNorm() > 0) {
    calin::math::vector3d_util::dump_as_proto(camera_offset,
      telescope->mutable_camera_offset());
  }

  Eigen::Quaterniond camera_rotation =
    calin::math::geometry::euler_to_quaternion(param.camera_rotation());
  if(param.camera_rotation_dispersion() >= 0) {
    camera_rotation =
      calin::math::geometry::euler_to_quaternion(
        calin::math::geometry::scattering_euler(param.camera_rotation_dispersion(),
          *rng, calin::ix::common_types::EulerAngles3D::YXY)) * camera_rotation;
  }
  if(camera_rotation.vec().squaredNorm() > 0) {
    telescope->mutable_camera_rotation()->set_rotation_order(
      calin::ix::common_types::EulerAngles3D::YXY);
    calin::math::geometry::quaternion_to_euler(
      telescope->mutable_camera_rotation(), camera_rotation);
  }
  telescope->set_camera_radius(param.camera_radius());

  telescope->mutable_camera_module_grid()->CopyFrom(param.camera_module_grid());
  telescope->mutable_camera_module_pixel_grid()->CopyFrom(param.camera_module_pixel_grid());

  std::set<unsigned> module_missing(
    param.camera_modules_missing().begin(),param.camera_modules_missing().end());
  unsigned pix_id = 0;
  unsigned mod_id = 0;
  unsigned grid_id = 0;
  for(unsigned iy=0; iy<param.camera_module_grid().num_side(); ++iy) {
    for(unsigned ix=0; ix<param.camera_module_grid().num_side(); ++ix) {
      double x;
      double z;
      calin::math::geometry::square_grid_site_center(x,z,grid_id,
        param.camera_module_grid().spacing(),param.camera_module_grid().num_side(),
        param.camera_module_grid().x_center(), param.camera_module_grid().y_center());
      if((param.camera_module_max_radius()<=0 or std::sqrt(x*x+z*z)<param.camera_module_max_radius())
          and (not module_missing.count(grid_id))) {
        auto* module = telescope->add_camera_modules();
        module->set_id(mod_id);
        module->set_grid_id(grid_id);
        module->set_first_pixel_id(pix_id);
        module->mutable_center_position()->set_x(x);
        module->mutable_center_position()->set_y(calin::math::least_squares::polyval(
            telescope->camera_surface_polynomial().data(),
            telescope->camera_surface_polynomial_size(), x*x+z*z));
        module->mutable_center_position()->set_z(z);

        pix_id += SQR(param.camera_module_pixel_grid().num_side());
        ++mod_id;
      }
      ++grid_id;
    }
  }

  // ***************************************************************************
  // ***************************************************************************
  //
  // OBSCURATIONS
  //
  // ***************************************************************************
  // ***************************************************************************

  for(const auto& obs_param : param.primary_obscuration()) {
    telescope->add_primary_obscuration()->CopyFrom(obs_param);
  }
  for(const auto& obs_param : param.secondary_obscuration()) {
    telescope->add_secondary_obscuration()->CopyFrom(obs_param);
  }
  for(const auto& obs_param : param.camera_obscuration()) {
    telescope->add_camera_obscuration()->CopyFrom(obs_param);
  }

  return telescope;
}

calin::ix::simulation::sct_optics::SCTArray*
calin::simulation::sct_optics::make_sct_array(
  calin::ix::simulation::sct_optics::SCTRandomArrayParameters& param,
  calin::math::rng::RNG* rng,
  calin::ix::simulation::sct_optics::SCTArray* array)
{
  calin::math::rng::RNG* my_rng = nullptr;
  if(rng == nullptr) {
    my_rng = rng = new calin::math::rng::RNG(__PRETTY_FUNCTION__, "Random SCT array generation");
  }
  if(array == nullptr) {
    array = new calin::ix::simulation::sct_optics::SCTArray();
  } else {
    array->Clear();
  }
  array->mutable_array_origin()->CopyFrom(param.array_origin());

  for(unsigned iscope=0;iscope<param.array_layout().scope_positions_size();++iscope)
  {
    const auto& pos = param.array_layout().scope_positions(iscope);
    auto* telescope = array->add_telescope();
    telescope->set_id(iscope);
    make_sct_telescope(param, pos.x(), pos.y(), pos.z(), rng, telescope);
  }
  delete my_rng;
  return array;
}

void calin::simulation::sct_optics::pixel_centers(
  const calin::ix::simulation::sct_optics::SCTTelescope& scope_params,
  Eigen::VectorXd& x_out, Eigen::VectorXd& z_out)
{
  unsigned nmodpix = SQR(scope_params.camera_module_pixel_grid().num_side());
  unsigned npix = scope_params.camera_modules_size()*nmodpix;
  x_out.resize(npix);
  z_out.resize(npix);
  for(const auto& m : scope_params.camera_modules()) {
    const double xm =
      m.center_position().x() + scope_params.camera_module_pixel_grid().x_center();
    const double ym =
      m.center_position().z() + scope_params.camera_module_pixel_grid().y_center();
    for(unsigned imodpix=0;imodpix<nmodpix;++imodpix) {
      unsigned ipix = m.first_pixel_id() + imodpix;
      if(ipix >= npix) {
        throw std::runtime_error("Pixel number out of range : "+std::to_string(ipix));
      }
      calin::math::geometry::square_grid_site_center(x_out[ipix], z_out[ipix], imodpix,
        scope_params.camera_module_pixel_grid().spacing(),
        scope_params.camera_module_pixel_grid().num_side(), xm, ym);
    }
  }
}
