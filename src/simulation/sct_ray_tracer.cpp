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
#include <simulation/sct_facet_scheme.hpp>
#include <simulation/sct_ray_tracer.hpp>

using namespace calin::simulation::sct_optics;
using namespace calin::util::log;

SCTRayTracer::SCTRayTracer(const calin::ix::simulation::sct_optics::SCTArray* array,
    calin::math::rng::RNG* rng, bool adopt_array, bool adopt_rng):
  array_(array), rng_(rng ? rng : new calin::math::rng::RNG(__PRETTY_FUNCTION__,
    "SCT ray tracer")), adopt_array_(adopt_array), adopt_rng_(rng ? adopt_rng : true)
{
  // nothing to see here
}

SCTRayTracer::~SCTRayTracer()
{
  if(adopt_array_)delete array_;
  if(adopt_rng_)delete rng_;
}

void SCTRayTracer::trace_ray_in_reflector_frame(unsigned iscope, calin::math::ray::Ray& ray)
{
  
}

calin::ix::simulation::sct_optics::SCTTelescope*
calin::simulation::sct_optics::make_sct_telescope(
  calin::ix::simulation::sct_optics::SCTRandomArrayParameters& param,
  double x, double y, double z,
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
    x += param.telescope_position_xy_dispersion()*rng->normal();
    y += param.telescope_position_xy_dispersion()*rng->normal();
  }
  if(param.telescope_position_z_dispersion() > 0) {
    z += param.telescope_position_z_dispersion()*rng->normal();
  }

  telescope->mutable_position()->set_x(x);
  telescope->mutable_position()->set_y(y);
  telescope->mutable_position()->set_z(z);

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
          *rng, calin::ix::common_types::EulerAngles3D::XYX)) * primary_rotation;
  }
  if(primary_rotation.vec().squaredNorm() > 0) {
    telescope->mutable_primary_rotation()->set_rotation_order(
      calin::ix::common_types::EulerAngles3D::XYX);
    calin::math::geometry::quaternion_to_euler(
      telescope->mutable_primary_rotation(), primary_rotation);
  }
  telescope->mutable_primary_facet_scheme()->CopyFrom(param.primary_facet_scheme());

  auto primary_facet_scheme = calin::simulation::sct_optics::SCTPrimaryFacetScheme(
    param.primary_facet_scheme());
  for(unsigned ifacet=0; ifacet<primary_facet_scheme.num_facets(); ++ifacet) {
    auto* facet = telescope->add_primary_facets();
    facet->set_id(ifacet);
    Eigen::Vector3d facet_position(0,0,0);
    primary_facet_scheme.facet_centroid(ifacet, facet_position.x(), facet_position.z());
    facet_position.y() = calin::math::least_squares::polyval(
      telescope->primary_surface_polynomial().data(),
      telescope->primary_surface_polynomial_size(), facet_position.squaredNorm());
    if(param.primary_facet_offset_xz_dispersion() > 0) {
      facet_position.x() += param.primary_facet_offset_xz_dispersion()*rng->normal();
      facet_position.z() += param.primary_facet_offset_xz_dispersion()*rng->normal();
    }
    if(param.primary_facet_offset_y_dispersion() > 0) {
      facet_position.y() += param.primary_facet_offset_y_dispersion()*rng->normal();
    }
    calin::math::vector3d_util::dump_as_proto(facet_position, facet->mutable_position());
    if(param.primary_facet_rotation_dispersion() > 0) {
      auto* facet_euler = facet->mutable_rotation();
      facet_euler->set_rotation_order(calin::ix::common_types::EulerAngles3D::XYX);
      calin::math::geometry::scattering_euler(facet_euler,
        param.primary_facet_rotation_dispersion(), *rng);
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
          *rng, calin::ix::common_types::EulerAngles3D::XYX)) * secondary_rotation;
  }
  if(secondary_rotation.vec().squaredNorm() > 0) {
    telescope->mutable_secondary_rotation()->set_rotation_order(
      calin::ix::common_types::EulerAngles3D::XYX);
    calin::math::geometry::quaternion_to_euler(
      telescope->mutable_secondary_rotation(), secondary_rotation);
  }
  telescope->mutable_secondary_facet_scheme()->CopyFrom(param.secondary_facet_scheme());

  auto secondary_facet_scheme = calin::simulation::sct_optics::SCTSecondaryFacetScheme(
    param.secondary_facet_scheme());
  for(unsigned ifacet=0; ifacet<secondary_facet_scheme.num_facets(); ++ifacet) {
    auto* facet = telescope->add_secondary_facets();
    facet->set_id(ifacet);
    Eigen::Vector3d facet_position(0,0,0);
    secondary_facet_scheme.facet_centroid(ifacet, facet_position.x(), facet_position.z());
    facet_position.y() = calin::math::least_squares::polyval(
      telescope->secondary_surface_polynomial().data(),
      telescope->secondary_surface_polynomial_size(), facet_position.squaredNorm());
    if(param.secondary_facet_offset_xz_dispersion() > 0) {
      facet_position.x() += param.secondary_facet_offset_xz_dispersion()*rng->normal();
      facet_position.z() += param.secondary_facet_offset_xz_dispersion()*rng->normal();
    }
    if(param.secondary_facet_offset_y_dispersion() > 0) {
      facet_position.y() += param.secondary_facet_offset_y_dispersion()*rng->normal();
    }
    calin::math::vector3d_util::dump_as_proto(facet_position, facet->mutable_position());
    if(param.secondary_facet_rotation_dispersion() > 0) {
      auto* facet_euler = facet->mutable_rotation();
      facet_euler->set_rotation_order(calin::ix::common_types::EulerAngles3D::XYX);
      calin::math::geometry::scattering_euler(facet_euler,
        param.secondary_facet_rotation_dispersion(), *rng);
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
          *rng, calin::ix::common_types::EulerAngles3D::XYX)) * camera_rotation;
  }
  if(camera_rotation.vec().squaredNorm() > 0) {
    telescope->mutable_camera_rotation()->set_rotation_order(
      calin::ix::common_types::EulerAngles3D::XYX);
    calin::math::geometry::quaternion_to_euler(
      telescope->mutable_camera_rotation(), camera_rotation);
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
