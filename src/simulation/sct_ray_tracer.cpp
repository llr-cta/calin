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
#include <math/geometry.hpp>
#include <math/vector3d_util.hpp>
#include <simulation/sct_facet_scheme.hpp>
#include <simulation/sct_ray_tracer.hpp>

using namespace calin::simulation::sct_optics;
using namespace calin::util::log;

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

  for(unsigned i=0;i<param.primary_sag_polynomial_size();++i) {
    double p_i = param.primary_sag_polynomial(i);
    telescope->add_primary_surface_polynomial(p_i);
  }
  Eigen::Vector3d primary_offset =
    calin::math::vector3d_util::from_proto(telescope->primary_offset());
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
    Eigen::Vector3d facet_offset(0,0,0);
    if(param.primary_facet_offset_xz_dispersion() > 0) {
      facet_offset.x() += param.primary_facet_offset_xz_dispersion()*rng->normal();
      facet_offset.z() += param.primary_facet_offset_xz_dispersion()*rng->normal();
      calin::math::vector3d_util::dump_as_proto(facet_offset, facet->mutable_offset());
    }
    if(param.primary_facet_offset_y_dispersion() > 0) {
      facet_offset.y() += param.primary_facet_offset_y_dispersion()*rng->normal();
      calin::math::vector3d_util::dump_as_proto(facet_offset, facet->mutable_offset());
    }
    if(param.primary_facet_rotation_dispersion() > 0) {
      auto* facet_euler = facet->mutable_rotation();
      facet_euler->set_rotation_order(calin::ix::common_types::EulerAngles3D::XYX);
      calin::math::geometry::scattering_euler(facet_euler,
        param.primary_facet_rotation_dispersion(), *rng);
    }
  }
  for(auto id : param.primary_facets_removed()) {
    if(id < primary_facet_scheme.num_facets()) {
      telescope->mutable_primary_facets(id)->set_removed(true);
    }
  }

/*
  message SCTFacet
  {
    uint32 id                                                = 1
      [(CFO).desc = "Sequential facet identification number, starting at zero."];
    bool removed                                             = 2
      [(CFO).desc = "True if this mirror facet has been removed."];
    double spot_size                                         = 3 [
      (CFO).desc = "Parameter describing spot size from this mirror in "
        "perfectly aligned system.", (CFO).units="cm" ];
    calin.ix.common_types.Vector3D offset                    = 4 [
      (CFO).desc = "Facet offset, if needed.", (CFO).units="cm" ];
    calin.ix.common_types.EulerAngles3D rotation             = 5 [
      (CFO).desc = "Facet rotation angles, if needed." ];
  };
*/

  // calin.ix.common_types.Matrix3D primary_rotation          = 12 [
  //   (CFO).desc = "Primary rotation matrix, if needed." ];
  // SCTFacetScheme primary_facet_scheme                      = 13 [
  //   (CFO).desc = "Primary facet scheme." ];
  // repeated SCTFacet primary_facets                         = 14 [
  //   (CFO).desc = "Primary facet parameters." ];
  //
  //
  // repeated double primary_surface_polynomial               = 10 [
  //   (CFO).desc = "Coefficients describing primary surface sag. Coefficients "
  //     "are function of squared-distance from optic axis in cm.",
  //   (CFO).units = "cm" ];
  // calin.ix.common_types.Vector3D primary_offset            = 11 [
  //   (CFO).desc = "Primary origin in reflector frame, if needed.", (CFO).units="cm" ];
  // calin.ix.common_types.Matrix3D primary_rotation          = 12 [
  //   (CFO).desc = "Primary rotation matrix, if needed." ];
  // SCTFacetScheme primary_facet_scheme                      = 13 [
  //   (CFO).desc = "Primary facet scheme." ];
  // repeated SCTFacet primary_facets                         = 14 [
  //   (CFO).desc = "Primary facet parameters." ];


  return telescope;
}

calin::ix::simulation::sct_optics::SCTArray*
calin::simulation::sct_optics::make_sct_array(
  calin::ix::simulation::sct_optics::SCTRandomArrayParameters& param,
  calin::math::rng::RNG* rng,
  calin::ix::simulation::sct_optics::SCTArray* array)
{
  if(rng == nullptr) {
    rng = new calin::math::rng::RNG(__PRETTY_FUNCTION__, "Random SCT array generation");
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

  return array;
}
