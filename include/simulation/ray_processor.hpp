/*

   calin/simulation/ray_processor.hpp -- Stephen Fegan -- 2017-01-16

   Multi-purpose ray processor.

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#pragma once

#include <math/ray.hpp>

namespace calin { namespace simulation { namespace ray_processor {

class RayProcessor
{
public:
  virtual ~RayProcessor();
  virtual void start_processing();
  virtual void process_ray(unsigned scope_id, const calin::math::ray::Ray& ray,
    double pe_weight);
  virtual void finish_processing();
};

} } } // namespace calin::simulation::ray_processor
