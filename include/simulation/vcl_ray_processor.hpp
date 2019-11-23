/*

   calin/simulation/vcl_ray_processor.hpp -- Stephen Fegan -- 2019-11-15

   Class for processing rays through VCL raytacer

   Copyright 2019, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <algorithm>
#include <limits>

#include <util/vcl.hpp>
#include <math/ray_vcl.hpp>
#include <math/rng_vcl.hpp>
#include <math/geometry_vcl.hpp>
#include <simulation/vcl_raytracer.hpp>
#include <simulation/ray_processor.hpp>

namespace calin { namespace simulation { namespace vcl_ray_processor {

template<typename VCLRealType> class VCLRayTracerRayProcessor:
  public calin::simulation::ray_processor::RayProcessor
{
public:
  CALIN_TYPEALIAS(real_t, typename VCLRealType::real_t);
  CALIN_TYPEALIAS(real_vt, typename VCLRealType::real_vt);
  CALIN_TYPEALIAS(real_at, typename VCLRealType::real_at);
  CALIN_TYPEALIAS(bool_vt, typename VCLRealType::bool_vt);
  CALIN_TYPEALIAS(int_vt, typename VCLRealType::int_vt);
  CALIN_TYPEALIAS(int_at, typename VCLRealType::int_at);

  CALIN_TYPEALIAS(RNG, calin::math::rng::VCLRealRNG<VCLRealType>);
  CALIN_TYPEALIAS(RayTracer, calin::simulation::vcl_raytracer::VCLScopeRayTracer<VCLRealType>);

  VCLRayTracerRayProcessor(calin::simulation::vs_optics::VSOArray* array,
      calin::simulation::pe_processor::PEProcessor* visitor, RNG* rng = nullptr,
      bool adopt_array = false, bool adopt_visitor = false,
      bool adopt_rng = false):
    array_(array), adopt_array_(adopt_array),
    visitor_(visitor), adopt_visitor_(adopt_visitor),
    rng_(rng==nullptr ? new RNG(__PRETTY_FUNCTION__) : rng),
    adopt_rng_(rng==nullptr ? true : adopt_rng),
    ray_tracer_(array->numTelescopes()), nray_(array->numTelescopes()),
    ray_ux_(array->numTelescopes() * buffer_depth_ * VCLRealType::num_real),
    ray_uy_(array->numTelescopes() * buffer_depth_ * VCLRealType::num_real),
    ray_uz_(array->numTelescopes() * buffer_depth_ * VCLRealType::num_real),
    ray_e_(array->numTelescopes() * buffer_depth_ * VCLRealType::num_real),
    ray_x_(array->numTelescopes() * buffer_depth_ * VCLRealType::num_real),
    ray_y_(array->numTelescopes() * buffer_depth_ * VCLRealType::num_real),
    ray_z_(array->numTelescopes() * buffer_depth_ * VCLRealType::num_real),
    ray_ct_(array->numTelescopes() * buffer_depth_ * VCLRealType::num_real),
    ray_w_(array->numTelescopes() * buffer_depth_ * VCLRealType::num_real)
  {
    for(unsigned iscope=0;iscope<array->numTelescopes();++iscope) {
      ray_tracer_[iscope] = new RayTracer(
        array->telescope(iscope), ref_index_, rng_, /* adopt_rng= */ false);
    }
  }

  virtual ~VCLRayTracerRayProcessor()
  {
    if(adopt_array_)delete array_;
    if(adopt_visitor_)delete visitor_;
    if(adopt_rng_)delete rng_;
  }

  std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere> detector_spheres() override
  {
    // Copied from VSORayProcessor::detector_spheres() - Ugh!
    std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere> s;
    for(unsigned iscope=0; iscope<array_->numTelescopes(); iscope++)
    {
      auto* scope = array_->telescope(iscope);
      Eigen::Vector3d sphere_center = scope->reflectorIPCenter();
      scope->reflectorToGlobal_pos(sphere_center);
      s.emplace_back(sphere_center, 0.5*scope->reflectorIP());
    }
    return s;
  }

  void start_processing() override
  {
    visitor_->start_processing();
  }

  void process_ray(unsigned scope_id, const calin::math::ray::Ray& ray,
    double pe_weight) override
  {
    unsigned iray = scope_id * buffer_depth_ * VCLRealType::num_real + nray_[scope_id];
    ++nray_[scope_id];
    ray_ux_[iray] = ray.ux();
    ray_uy_[iray] = ray.uy();
    ray_uz_[iray] = ray.uz();
    ray_e_[iray] = ray.energy();
    ray_x_[iray] = ray.x();
    ray_y_[iray] = ray.y();
    ray_z_[iray] = ray.z();
    ray_ct_[iray] = ray.ct();
    ray_w_[iray] = pe_weight;
    if(nray_[scope_id] == buffer_depth_ * VCLRealType::num_real) {
      process_buffered_rays(scope_id);
    }
  }

  void finish_processing() override
  {
    for(unsigned iscope=0; iscope<array_->numTelescopes(); iscope++) {
      process_buffered_rays(iscope);
    }
    visitor_->finish_processing();
  }

private:
  void process_buffered_rays(unsigned scope_id) {
    using Ray = calin::math::ray::VCLRay<VCLRealType>;
    using TraceInfo = calin::simulation::vcl_raytracer::VCLScopeTraceInfo<VCLRealType>;
    unsigned nray = nray_[scope_id];
    unsigned iray = scope_id * buffer_depth_ * VCLRealType::num_real;
    nray_[scope_id] = 0;
    while(nray) {
      bool_vt mask = true;
      Ray ray;
      TraceInfo info;
      ray.mutable_ux().load(ray_ux_.data() + iray);
      ray.mutable_uy().load(ray_uy_.data() + iray);
      ray.mutable_uz().load(ray_uz_.data() + iray);
      ray.mutable_energy().load(ray_e_.data() + iray);
      ray.mutable_x().load(ray_x_.data() + iray);
      ray.mutable_y().load(ray_y_.data() + iray);
      ray.mutable_z().load(ray_z_.data() + iray);
      ray.mutable_ct() = 0;
      mask = ray_tracer_[scope_id]->
        trace_global_frame(mask, ray, info, /* do_derotation = */ false);

      real_at fp_x;
      real_at fp_z;
      real_at fp_t;
      int_at fp_pixel_id;
      unsigned imask = vcl::to_bits(mask);

      info.fplane_x.store(fp_x);
      info.fplane_z.store(fp_z);
      info.fplane_t.store(fp_t);
      info.pixel_id.store(fp_pixel_id);

      unsigned mray = std::min(VCLRealType::num_real, nray);
      for(unsigned jray=0; jray<mray; jray++) {
        if(imask & (0x1<<jray)) {
          visitor_->process_focal_plane_hit(scope_id, fp_pixel_id[jray],
            double(fp_x[jray]), double(fp_z[jray]),
            double(fp_t[jray]) + ray_ct_[iray+jray]*math::constants::cgs_1_c,
            ray_w_[iray+jray]);
        }
      }

      iray += VCLRealType::num_real;
      nray -= mray;
    }
  }

  calin::simulation::vs_optics::VSOArray* array_ = nullptr;
  bool adopt_array_ = false;
  calin::simulation::pe_processor::PEProcessor* visitor_ = nullptr;
  bool adopt_visitor_ = false;
  RNG* rng_ = nullptr;
  bool adopt_rng_ = false;
  std::vector<RayTracer*> ray_tracer_;
  double ref_index_ = 1.0;
  unsigned buffer_depth_ = 1;
  std::vector<unsigned> nray_;
  std::vector<real_t> ray_ux_;
  std::vector<real_t> ray_uy_;
  std::vector<real_t> ray_uz_;
  std::vector<real_t> ray_e_;
  std::vector<real_t> ray_x_;
  std::vector<real_t> ray_y_;
  std::vector<real_t> ray_z_;
  std::vector<double> ray_ct_;
  std::vector<double> ray_w_;
};

} } } // namespace calin::simulations::vcl_ray_processor
