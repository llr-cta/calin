/*

   calin/simulation/ray_processor.cpp -- Stephen Fegan -- 2017-01-16

   Multi-purpose ray (weight, scope, ray) processor.

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

#include <simulation/ray_processor.hpp>

using namespace calin::simulation::ray_processor;

RayProcessor::~RayProcessor()
{
  // nothing to see here
}

void RayProcessor::start_processing()
{
  // nothing to see here
}

void RayProcessor::process_ray(unsigned scope_id,
  const calin::math::ray::Ray& ray, double pe_weight)
{
  // nothing to see here
}

void RayProcessor::finish_processing()
{
  // nothing to see here
}
