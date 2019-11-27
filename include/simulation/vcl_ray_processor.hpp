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
#include <util/log.hpp>

namespace calin { namespace simulation { namespace vcl_ray_processor {

template<typename VCLArchitecture> class VCLRayTracerRayProcessorDouble:
  public calin::simulation::ray_processor::RayProcessor
{
public:
  CALIN_TYPEALIAS(double_vt, typename VCLArchitecture::double_vt);
  CALIN_TYPEALIAS(double_at, typename VCLArchitecture::double_at);
  CALIN_TYPEALIAS(int64_vt, typename VCLArchitecture::int64_vt);
  CALIN_TYPEALIAS(int64_at, typename VCLArchitecture::int64_at);
  CALIN_TYPEALIAS(bool_vt, typename VCLArchitecture::double_bvt);

  CALIN_TYPEALIAS(VCLRealType, calin::util::vcl::VCLDoubleReal<VCLArchitecture>);
  CALIN_TYPEALIAS(ArchRNG, calin::math::rng::VCLRNG<VCLArchitecture>);
  CALIN_TYPEALIAS(RNG, calin::math::rng::VCLRealRNG<VCLRealType>);
  CALIN_TYPEALIAS(RayTracer, calin::simulation::vcl_raytracer::VCLScopeRayTracer<VCLRealType>);

  VCLRayTracerRayProcessorDouble(calin::simulation::vs_optics::VSOArray* array,
      calin::simulation::pe_processor::PEProcessor* visitor, ArchRNG* rng = nullptr,
      bool adopt_array = false, bool adopt_visitor = false,
      bool adopt_rng = false):
    array_(array), adopt_array_(adopt_array),
    visitor_(visitor), adopt_visitor_(adopt_visitor),
    rng_(new RNG(rng==nullptr ? new ArchRNG(__PRETTY_FUNCTION__) : rng,
      rng==nullptr ? true : adopt_rng)), adopt_rng_(true),
    ray_tracer_(array->numTelescopes()), nray_(array->numTelescopes()),
    ray_ux_(array->numTelescopes() * buffer_depth_ * VCLArchitecture::num_double),
    ray_uy_(array->numTelescopes() * buffer_depth_ * VCLArchitecture::num_double),
    ray_uz_(array->numTelescopes() * buffer_depth_ * VCLArchitecture::num_double),
    ray_e_(array->numTelescopes() * buffer_depth_ * VCLArchitecture::num_double),
    ray_x_(array->numTelescopes() * buffer_depth_ * VCLArchitecture::num_double),
    ray_y_(array->numTelescopes() * buffer_depth_ * VCLArchitecture::num_double),
    ray_z_(array->numTelescopes() * buffer_depth_ * VCLArchitecture::num_double),
    ray_ct_(array->numTelescopes() * buffer_depth_ * VCLArchitecture::num_double),
    ray_w_(array->numTelescopes() * buffer_depth_ * VCLArchitecture::num_double)
  {
    for(unsigned iscope=0;iscope<array->numTelescopes();++iscope) {
      ray_tracer_[iscope] = new RayTracer(
        array->telescope(iscope), ref_index_, rng_, /* adopt_rng= */ false);
    }
  }

  virtual ~VCLRayTracerRayProcessorDouble()
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
    nhit_ = 0;
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

  uint64_t nhit() const { return nhit_; }

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

      double_at fp_x;
      double_at fp_z;
      double_at fp_t;
      int64_at fp_pixel_id;
      unsigned imask = vcl::to_bits(mask);

      info.fplane_x.store(fp_x);
      info.fplane_z.store(fp_z);
      info.fplane_t.store(fp_t);
      info.pixel_id.store(fp_pixel_id);

      unsigned mray = std::min(VCLRealType::num_real, nray);
      for(unsigned jray=0; jray<mray; jray++) {
        if(imask & (0x1<<jray)) {
          ++nhit_;
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
  std::vector<double> ray_ux_;
  std::vector<double> ray_uy_;
  std::vector<double> ray_uz_;
  std::vector<double> ray_e_;
  std::vector<double> ray_x_;
  std::vector<double> ray_y_;
  std::vector<double> ray_z_;
  std::vector<double> ray_ct_;
  std::vector<double> ray_w_;
  uint64_t nhit_ = 0;
};

template<typename VCLArchitecture> class VCLRayTracerRayProcessorFloat:
  public calin::simulation::ray_processor::RayProcessor
{
public:
  CALIN_TYPEALIAS(float_vt, typename VCLArchitecture::float_vt);
  CALIN_TYPEALIAS(double_vt, typename VCLArchitecture::double_vt);
  CALIN_TYPEALIAS(int32_vt, typename VCLArchitecture::int32_vt);
  CALIN_TYPEALIAS(int32_at, typename VCLArchitecture::int32_at);
  CALIN_TYPEALIAS(bool_vt, typename VCLArchitecture::float_bvt);

  CALIN_TYPEALIAS(FloatType, calin::util::vcl::VCLFloatReal<VCLArchitecture>);
  CALIN_TYPEALIAS(DoubleType, calin::util::vcl::VCLDoubleReal<VCLArchitecture>);
  CALIN_TYPEALIAS(ArchRNG, calin::math::rng::VCLRNG<VCLArchitecture>);
  CALIN_TYPEALIAS(FloatRNG, calin::math::rng::VCLRealRNG<FloatType>);
  CALIN_TYPEALIAS(FloatRayTracer, calin::simulation::vcl_raytracer::VCLScopeRayTracer<FloatType>);
  CALIN_TYPEALIAS(DoubleRayTracer, calin::simulation::vcl_raytracer::VCLScopeRayTracer<DoubleType>);

  VCLRayTracerRayProcessorFloat(calin::simulation::vs_optics::VSOArray* array,
      calin::simulation::pe_processor::PEProcessor* visitor, ArchRNG* rng = nullptr,
      bool adopt_array = false, bool adopt_visitor = false,
      bool adopt_rng = false):
    array_(array), adopt_array_(adopt_array),
    visitor_(visitor), adopt_visitor_(adopt_visitor),
    rng_(new FloatRNG(rng==nullptr ? new ArchRNG(__PRETTY_FUNCTION__) : rng,
      rng==nullptr ? true : adopt_rng)), adopt_rng_(true),
    ray_tracer_(array->numTelescopes()), nray_(array->numTelescopes()),
    ray_ux_(array->numTelescopes() * buffer_depth_ * VCLArchitecture::num_float),
    ray_uy_(array->numTelescopes() * buffer_depth_ * VCLArchitecture::num_float),
    ray_uz_(array->numTelescopes() * buffer_depth_ * VCLArchitecture::num_float),
    ray_e_(array->numTelescopes() * buffer_depth_ * VCLArchitecture::num_float),
    ray_x_(array->numTelescopes() * buffer_depth_ * VCLArchitecture::num_float),
    ray_y_(array->numTelescopes() * buffer_depth_ * VCLArchitecture::num_float),
    ray_z_(array->numTelescopes() * buffer_depth_ * VCLArchitecture::num_float),
    ray_ct_(array->numTelescopes() * buffer_depth_ * VCLArchitecture::num_float),
    ray_w_(array->numTelescopes() * buffer_depth_ * VCLArchitecture::num_float)
  {
    for(unsigned iscope=0;iscope<array->numTelescopes();++iscope) {
      ray_tracer_[iscope] = new FloatRayTracer(
        array->telescope(iscope), ref_index_, rng_, /* adopt_rng= */ false);
    }
  }

  virtual ~VCLRayTracerRayProcessorFloat()
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
    nhit_ = 0;
    visitor_->start_processing();
  }

  void process_ray(unsigned scope_id, const calin::math::ray::Ray& ray,
    double pe_weight) override
  {
    unsigned iray = scope_id * buffer_depth_ * VCLArchitecture::num_float + nray_[scope_id];
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
    if(nray_[scope_id] == buffer_depth_ * VCLArchitecture::num_float) {
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

  uint64_t nhit() const { return nhit_; }

private:
  void store(double* xa, float_vt xv) const {
    double_vt xlo = extend_low(xv);
    double_vt xhi = extend_high(xv);
    xlo.store(xa);
    xhi.store(xa + VCLArchitecture::num_double);
  }

  void process_buffered_rays(unsigned scope_id) {
    const auto* scope = array_->telescope(scope_id);

    using FloatRay = calin::math::ray::VCLRay<FloatType>;
    using DoubleRay = calin::math::ray::VCLRay<DoubleType>;

    using TraceInfo = calin::simulation::vcl_raytracer::VCLScopeTraceInfo<FloatType>;
    unsigned nray = nray_[scope_id];
    unsigned iray = scope_id * buffer_depth_ * VCLArchitecture::num_float;
    nray_[scope_id] = 0;
    while(nray) {
      DoubleRay dbl_ray_lo;
      dbl_ray_lo.mutable_ux().load(ray_ux_.data() + iray);
      dbl_ray_lo.mutable_uy().load(ray_uy_.data() + iray);
      dbl_ray_lo.mutable_uz().load(ray_uz_.data() + iray);
      dbl_ray_lo.mutable_energy().load(ray_e_.data() + iray);
      dbl_ray_lo.mutable_x().load(ray_x_.data() + iray);
      dbl_ray_lo.mutable_y().load(ray_y_.data() + iray);
      dbl_ray_lo.mutable_z().load(ray_z_.data() + iray);
      dbl_ray_lo.mutable_ct().load(ray_ct_.data() + iray);
      if(propagate_to_scope_environs_) {
        dbl_ray_lo.propagate_to_point_nearly_closest_approach_with_mask(/*mask=*/ true,
          scope->pos().cast<double_vt>(), -2*scope->focalPlanePosition().y(),
          /* time_reversal_ok = */ false, ref_index_);
        dbl_ray_lo.ct().store(ray_ct_.data() + iray);
      }
      DoubleRayTracer::transform_to_scope_reflector_frame(dbl_ray_lo, scope);
      iray += VCLArchitecture::num_double;

      DoubleRay dbl_ray_hi;
      dbl_ray_hi.mutable_ux().load(ray_ux_.data() + iray);
      dbl_ray_hi.mutable_uy().load(ray_uy_.data() + iray);
      dbl_ray_hi.mutable_uz().load(ray_uz_.data() + iray);
      dbl_ray_hi.mutable_energy().load(ray_e_.data() + iray);
      dbl_ray_hi.mutable_x().load(ray_x_.data() + iray);
      dbl_ray_hi.mutable_y().load(ray_y_.data() + iray);
      dbl_ray_hi.mutable_z().load(ray_z_.data() + iray);
      dbl_ray_hi.mutable_ct().load(ray_ct_.data() + iray);
      if(propagate_to_scope_environs_) {
        dbl_ray_hi.propagate_to_point_nearly_closest_approach_with_mask(/*mask=*/ true,
          scope->pos().cast<double_vt>(), -2*scope->focalPlanePosition().y(),
          /* time_reversal_ok = */ false, ref_index_);
        dbl_ray_hi.ct().store(ray_ct_.data() + iray);
      }
      DoubleRayTracer::transform_to_scope_reflector_frame(dbl_ray_hi, scope);
      iray -= VCLArchitecture::num_double;

      FloatRay flt_ray;
      flt_ray.mutable_ux() = compress(dbl_ray_lo.ux(), dbl_ray_hi.ux());
      flt_ray.mutable_uy() = compress(dbl_ray_lo.uy(), dbl_ray_hi.uy());
      flt_ray.mutable_uz() = compress(dbl_ray_lo.uz(), dbl_ray_hi.uz());
      flt_ray.mutable_energy() = compress(dbl_ray_lo.energy(), dbl_ray_hi.energy());
      flt_ray.mutable_x() = compress(dbl_ray_lo.x(), dbl_ray_hi.x());
      flt_ray.mutable_y() = compress(dbl_ray_lo.y(), dbl_ray_hi.y());
      flt_ray.mutable_z() = compress(dbl_ray_lo.z(), dbl_ray_hi.z());
      flt_ray.mutable_ct() = 0;

      bool_vt mask = true;
      TraceInfo info;

      mask = ray_tracer_[scope_id]->trace_reflector_frame(mask, flt_ray, info);

      double fp_x[VCLArchitecture::num_float];
      double fp_z[VCLArchitecture::num_float];
      double fp_t[VCLArchitecture::num_float];
      int32_t fp_pixel_id[VCLArchitecture::num_float];
      unsigned imask = vcl::to_bits(mask);

      store(fp_x, info.fplane_x);
      store(fp_z, info.fplane_z);
      store(fp_t, info.fplane_t);
      info.pixel_id.store(fp_pixel_id);

      unsigned mray = std::min(VCLArchitecture::num_float, nray);
      for(unsigned jray=0; jray<mray; jray++) {
        if(imask & (0x1<<jray)) {
          ++nhit_;
          visitor_->process_focal_plane_hit(scope_id, fp_pixel_id[jray],
            fp_x[jray], fp_z[jray], fp_t[jray] + ray_ct_[iray+jray]*math::constants::cgs_1_c,
            ray_w_[iray+jray]);
        }
      }

      iray += VCLArchitecture::num_float;
      nray -= mray;
    }
  }

  calin::simulation::vs_optics::VSOArray* array_ = nullptr;
  bool adopt_array_ = false;
  calin::simulation::pe_processor::PEProcessor* visitor_ = nullptr;
  bool adopt_visitor_ = false;
  FloatRNG* rng_ = nullptr;
  bool adopt_rng_ = false;
  std::vector<FloatRayTracer*> ray_tracer_;
  double ref_index_ = 1.0;
  double ref_index_inv_ = 1.0;
  unsigned buffer_depth_ = 1;
  bool propagate_to_scope_environs_ = true;
  std::vector<unsigned> nray_;
  std::vector<double> ray_ux_;
  std::vector<double> ray_uy_;
  std::vector<double> ray_uz_;
  std::vector<double> ray_e_;
  std::vector<double> ray_x_;
  std::vector<double> ray_y_;
  std::vector<double> ray_z_;
  std::vector<double> ray_ct_;
  std::vector<double> ray_w_;
  uint64_t nhit_ = 0;
};

} } } // namespace calin::simulations::vcl_ray_processor
