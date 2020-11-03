/*

   calin/simulation/vcl_raytracer.cpp -- Stephen Fegan -- 2020-11-30

   Some test functions for VCL raytracer components

   Copyright 2020, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <simulation/vcl_ray_processor.hpp>

using namespace calin::simulation::vcl_ray_processor;

SingleRayVCLScopeTraceInfoProcessor::~SingleRayVCLScopeTraceInfoProcessor()
{
  // nothing to see here
}

void SingleRayVCLScopeTraceInfoProcessor::start_processing()
{
  // nothing to see here
}

void SingleRayVCLScopeTraceInfoProcessor::process_vcl_scope_trace_info(
  unsigned scope_id, bool hit_fplane,
  const SingleRayVCLScopeTraceInfo& trace_info, double pe_weight)
{
  // nothing to see here
}

void SingleRayVCLScopeTraceInfoProcessor::finish_processing()
{
  // nothing to see here
}
