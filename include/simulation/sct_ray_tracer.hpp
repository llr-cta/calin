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
#include <simulation/sct_optics.pb.h>
#include <simulation/sct_facet_scheme.hpp>

namespace calin { namespace simulation { namespace sct_optics {

class SCTRayTracer
{
public:
  SCTRayTracer(const calin::ix::simulation::sct_optics::SCTArray* array,
    calin::math::rng::RNG* rng = nullptr, bool adopt_array = false, bool adopt_rng = false);
  ~SCTRayTracer();
  void trace_ray_in_reflector_frame(unsigned iscope, calin::math::ray::Ray& ray);
private:
  const calin::ix::simulation::sct_optics::SCTArray* array_ = nullptr;
  calin::math::rng::RNG* rng_ = nullptr;
  bool adopt_array_ = false;
  bool adopt_rng_ = false;
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
