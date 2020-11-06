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

void SingleRayVCLScopeTraceInfoProcessor::
process_vcl_scope_trace_info(const SingleRayVCLScopeTraceInfo& trace_info)
{
  // nothing to see here
}

void SingleRayVCLScopeTraceInfoProcessor::finish_processing()
{
  // nothing to see here
}

RecordingSingleRayVCLScopeTraceInfoProcessor::RecordingSingleRayVCLScopeTraceInfoProcessor():
  SingleRayVCLScopeTraceInfoProcessor()
{
  // nothing to see here
}

RecordingSingleRayVCLScopeTraceInfoProcessor::~RecordingSingleRayVCLScopeTraceInfoProcessor()
{
  // nothing to see here
}

void RecordingSingleRayVCLScopeTraceInfoProcessor::start_processing()
{
  traces_.clear();
}

void RecordingSingleRayVCLScopeTraceInfoProcessor::
process_vcl_scope_trace_info(const SingleRayVCLScopeTraceInfo& trace_info)
{
  traces_.emplace_back(trace_info);
}

RayMapSingleRayVCLScopeTraceInfoProcessor::
RayMapSingleRayVCLScopeTraceInfoProcessor(double xside, unsigned nside, MapQuantity map_quantity):
  SingleRayVCLScopeTraceInfoProcessor(),
  dx_inv_(double(nside)/xside), xside_2_(0.5*xside), nside_(nside),
  map_quantity_(map_quantity), hist_(nside, nside), weight_(nside, nside)
{
  // nothing to see here
}

RayMapSingleRayVCLScopeTraceInfoProcessor::~RayMapSingleRayVCLScopeTraceInfoProcessor()
{
  // nothing to see here
}

void RayMapSingleRayVCLScopeTraceInfoProcessor::start_processing()
{
  hist_.setZero();
  weight_.setZero();
}

void RayMapSingleRayVCLScopeTraceInfoProcessor::
process_vcl_scope_trace_info(const SingleRayVCLScopeTraceInfo& trace_info)
{
  double x;
  double y;
  double q;
  switch(map_quantity_) {
  case MQ_WEIGHT_BY_FOCAL_PLANE_POSITION:
    x = trace_info.fplane_x;
    y = trace_info.fplane_z;
    q = 1.0;
    break;
  case MQ_WEIGHT_BY_REFLECTOR_SPHERE_POSITION:
    x = trace_info.reflec_x;
    y = trace_info.reflec_z;
    q = 1.0;
    break;
  case MQ_WEIGHT_BY_MIRROR_FACET_POSITION:
    x = trace_info.mirror_x;
    y = trace_info.mirror_z;
    q = 1.0;
    break;
  case MQ_MIRROR_INCIDENCE_ANGLE_BY_FOCAL_PLANE_POSITION:
    x = trace_info.fplane_x;
    y = trace_info.fplane_z;
    q = trace_info.mirror_n_dot_u;
    break;
  case MQ_MIRROR_INCIDENCE_ANGLE_BY_MIRROR_FACET_POSITION:
    x = trace_info.mirror_x;
    y = trace_info.mirror_z;
    q = trace_info.mirror_n_dot_u;
    break;
  case MQ_FOCAL_PLANE_INCIDENCE_ANGLE_BY_FOCAL_PLANE_POSITION:
    x = trace_info.fplane_x;
    y = trace_info.fplane_z;
    q = trace_info.fplane_uy;
    break;
  case MQ_FOCAL_PLANE_INCIDENCE_ANGLE_BY_MIRROR_FACET_POSITION:
    x = trace_info.mirror_x;
    y = trace_info.mirror_z;
    q = trace_info.fplane_uy;
    break;
  }
  int ix = std::round((x + xside_2_)*dx_inv_ );
  int iy = std::round((y + xside_2_)*dx_inv_);
  if(ix>=0 and iy>=0 and ix<nside_ and iy<nside_) {
    hist_(ix,iy) += q * trace_info.ray_weight;
    weight_(ix,iy) += trace_info.ray_weight;
  }
}

void RayMapSingleRayVCLScopeTraceInfoProcessor::finish_processing()
{
  switch(map_quantity_) {
  case MQ_WEIGHT_BY_FOCAL_PLANE_POSITION:
  case MQ_WEIGHT_BY_REFLECTOR_SPHERE_POSITION:
  case MQ_WEIGHT_BY_MIRROR_FACET_POSITION:
    // no post-processing needed
    break;
  case MQ_MIRROR_INCIDENCE_ANGLE_BY_FOCAL_PLANE_POSITION:
  case MQ_MIRROR_INCIDENCE_ANGLE_BY_MIRROR_FACET_POSITION:
  case MQ_FOCAL_PLANE_INCIDENCE_ANGLE_BY_FOCAL_PLANE_POSITION:
  case MQ_FOCAL_PLANE_INCIDENCE_ANGLE_BY_MIRROR_FACET_POSITION:
    hist_.array() = (hist_.array().abs() / weight_.array()).acos()*180.0/M_PI;
    break;
  }
}
