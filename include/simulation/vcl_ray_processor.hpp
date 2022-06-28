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
#include <simulation/vso_array.hpp>
#include <simulation/vcl_raytracer.hpp>
#include <simulation/ray_processor.hpp>
#include <util/log.hpp>

namespace calin { namespace simulation { namespace vcl_ray_processor {

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//     THESE CLASSES DEFINE AND IMPLEMENT THE VCL VECORIZED RAY PROCESSOR     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<typename VCLArchitecture> struct alignas(VCLArchitecture::vec_bytes) VCLFocalPlaneParameters
{
#ifndef SWIG
  using int64_vt    = typename VCLArchitecture::int64_vt;
  using double_bvt  = typename VCLArchitecture::double_bvt;
  using double_vt   = typename VCLArchitecture::double_vt;
  using Vector3d_vt = typename VCLArchitecture::Vector3d_vt;
  using Ray_vt      = typename calin::math::ray::VCLRay<VCLArchitecture>;
#endif // not defined SWIG
  double_vt fplane_x;
  double_vt fplane_y;
  double_vt fplane_z;
  double_vt fplane_ux;
  double_vt fplane_uy;
  double_vt fplane_uz;
  double_vt fplane_t;
  int64_vt pixel_id;
};

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCLFocalPlaneRayPropagator
{
public:
#ifndef SWIG
  using int64_vt    = typename VCLArchitecture::int64_vt;
  using double_bvt  = typename VCLArchitecture::double_bvt;
  using double_vt   = typename VCLArchitecture::double_vt;
  using Vector3d_vt = typename VCLArchitecture::Vector3d_vt;
  using Real_vt     = calin::util::vcl::VCLDoubleReal<VCLArchitecture>;
  using Ray_vt      = calin::math::ray::VCLRay<Real_vt>;
#endif // not defined SWIG

  virtual ~VCLFocalPlaneRayPropagator() {
    // nothing to see here
  }

  virtual std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere> detector_spheres() {
    return {};
  }

  virtual void start_propagating() {
    // nothing to see here
  }

  virtual double_bvt propagate_rays_to_focal_plane(
      unsigned scope_id, Ray_vt& ray, double_bvt ray_mask,
      VCLFocalPlaneParameters<VCLArchitecture>& fp_parameters) {
    return false;
  }

  virtual void finish_propagating() {
    // nothing to see here
  }
};

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) DaviesCottonVCLFocalPlaneRayPropagator:
  public VCLFocalPlaneRayPropagator<VCLArchitecture>
{
public:
#ifndef SWIG
  using int64_vt    = typename VCLArchitecture::int64_vt;
  using double_bvt  = typename VCLArchitecture::double_bvt;
  using double_vt   = typename VCLArchitecture::double_vt;
  using Vector3d_vt = typename VCLArchitecture::Vector3d_vt;
  using Real_vt     = calin::util::vcl::VCLDoubleReal<VCLArchitecture>;
  using Ray_vt      = typename calin::math::ray::VCLRay<Real_vt>;

  using ArchRNG     = calin::math::rng::VCLRNG<VCLArchitecture>;
  using RealRNG     = calin::math::rng::VCLRealRNG<Real_vt>;
  using RayTracer   = calin::simulation::vcl_raytracer::VCLScopeRayTracer<Real_vt>;
  using TraceInfo   = calin::simulation::vcl_raytracer::VCLScopeTraceInfo<Real_vt>;
#endif // not defined SWIG

  DaviesCottonVCLFocalPlaneRayPropagator(calin::simulation::vs_optics::VSOArray* array,
      ArchRNG* rng = nullptr,
      bool adopt_array = false, bool adopt_visitor = false,
      bool adopt_rng = false):
    array_(array), adopt_array_(adopt_array),
    rng_(new RealRNG(rng==nullptr ? new ArchRNG(__PRETTY_FUNCTION__) : rng,
      rng==nullptr ? true : adopt_rng)), adopt_rng_(true),
    ray_tracer_(array->numTelescopes())
  {
    for(unsigned iscope=0;iscope<array->numTelescopes();++iscope) {
      ray_tracer_[iscope] = new RayTracer(
        array->telescope(iscope), ref_index_, rng_, /* adopt_rng= */ false);
    }
  }

  virtual ~DaviesCottonVCLFocalPlaneRayPropagator() {
    for(auto* rt : ray_tracer_)delete rt;
  }

  virtual std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere> detector_spheres() {
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

  virtual void start_propagating() {
    // nothing to see here
  }

  virtual double_bvt propagate_rays_to_focal_plane(
      unsigned scope_id, Ray_vt& ray, double_bvt ray_mask,
      VCLFocalPlaneParameters<VCLArchitecture>& fp_parameters) {
    TraceInfo info;
    ray_mask = ray_tracer_[scope_id]->
      trace_global_frame(ray_mask, ray, info, /* do_derotation = */ false);
    fp_parameters.fplane_x     = info.fplane_x;
    fp_parameters.fplane_y     = info.fplane_y;
    fp_parameters.fplane_z     = info.fplane_z;
    fp_parameters.fplane_ux    = info.fplane_ux;
    fp_parameters.fplane_uy    = info.fplane_uy;
    fp_parameters.fplane_uz    = info.fplane_uz;
    fp_parameters.fplane_t     = info.fplane_t;
    fp_parameters.pixel_id     = info.pixel_id;
    return ray_mask;
  }

  virtual void finish_propagating() {
    // nothing to see here
  }

private:
  calin::simulation::vs_optics::VSOArray* array_ = nullptr;
  bool adopt_array_ = false;
  unsigned vcl_sti_visitor_status_min_ = 0;
  RealRNG* rng_ = nullptr;
  bool adopt_rng_ = false;
  std::vector<RayTracer*> ray_tracer_;
  double ref_index_ = 1.0;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//   THESE CLASSES IMPLEMENT THE SCALAR RayProcessor USING THE VCLRayTracer   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class SingleRayVCLScopeTraceInfo
{
public:
  bool             hit_fplane;
  double           ray_weight;
  unsigned         scope_id;

  int64_t          status;

  double           reflec_x;
  double           reflec_y;
  double           reflec_z;

  uint64_t         pre_reflection_obs_hitmask;
  uint64_t         post_reflection_obs_hitmask;
  uint64_t         camera_obs_hitmask;

  int64_t          mirror_hexid;
  int64_t          mirror_id;
  double           mirror_x;
  double           mirror_y;
  double           mirror_z;
  double           mirror_n_dot_u;

  double           fplane_x;
  double           fplane_z;
  double           fplane_t;

  double           fplane_ux;
  double           fplane_uy;
  double           fplane_uz;

  int64_t          pixel_hexid;
  int64_t          pixel_id;
};

class SingleRayVCLScopeTraceInfoProcessor
{
public:
  virtual ~SingleRayVCLScopeTraceInfoProcessor();
  virtual void start_processing();
  virtual void process_vcl_scope_trace_info(const SingleRayVCLScopeTraceInfo& trace_info);
  virtual void finish_processing();
};

class RecordingSingleRayVCLScopeTraceInfoProcessor: public SingleRayVCLScopeTraceInfoProcessor
{
public:
  RecordingSingleRayVCLScopeTraceInfoProcessor();
  virtual ~RecordingSingleRayVCLScopeTraceInfoProcessor();
  void start_processing() override;
  void process_vcl_scope_trace_info(const SingleRayVCLScopeTraceInfo& trace_info) override;
  unsigned ntrace() const { return traces_.size(); }
  SingleRayVCLScopeTraceInfo trace(unsigned itrace) const { return traces_.at(itrace); }
private:
  std::vector<SingleRayVCLScopeTraceInfo> traces_;
};

enum MapQuantity {
  MQ_WEIGHT_BY_FOCAL_PLANE_POSITION,                        // 0
  MQ_WEIGHT_BY_REFLECTOR_SPHERE_POSITION,                   // 1
  MQ_WEIGHT_BY_MIRROR_FACET_POSITION,                       // 2
  MQ_MIRROR_INCIDENCE_ANGLE_BY_FOCAL_PLANE_POSITION,        // 3
  MQ_MIRROR_INCIDENCE_ANGLE_BY_MIRROR_FACET_POSITION,       // 4
  MQ_FOCAL_PLANE_INCIDENCE_ANGLE_BY_FOCAL_PLANE_POSITION,   // 5
  MQ_FOCAL_PLANE_INCIDENCE_ANGLE_BY_MIRROR_FACET_POSITION   // 6
};

class RayMapSingleRayVCLScopeTraceInfoProcessor: public SingleRayVCLScopeTraceInfoProcessor
{
public:
  RayMapSingleRayVCLScopeTraceInfoProcessor(double xside, unsigned nside,
    MapQuantity map_quantity = MQ_WEIGHT_BY_FOCAL_PLANE_POSITION);
  virtual ~RayMapSingleRayVCLScopeTraceInfoProcessor();
  void start_processing() override;
  void process_vcl_scope_trace_info(const SingleRayVCLScopeTraceInfo& trace_info) override;
  void finish_processing() override;
  const Eigen::MatrixXd& hist() const { return hist_; }
  const Eigen::MatrixXd& weight() const { return weight_; }
private:
  double dx_inv_;
  double xside_2_;
  int nside_;
  MapQuantity map_quantity_;
  Eigen::MatrixXd hist_;
  Eigen::MatrixXd weight_;
};

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCLRayTracerRayProcessorDouble:
  public calin::simulation::ray_processor::RayProcessor
{
public:
  CALIN_TYPEALIAS(double_vt, typename VCLArchitecture::double_vt);
  CALIN_TYPEALIAS(double_at, typename VCLArchitecture::double_at);
  CALIN_TYPEALIAS(int64_vt, typename VCLArchitecture::int64_vt);
  CALIN_TYPEALIAS(int64_at, typename VCLArchitecture::int64_at);
  CALIN_TYPEALIAS(uint64_vt, typename VCLArchitecture::uint64_vt);
  CALIN_TYPEALIAS(uint64_at, typename VCLArchitecture::uint64_at);
  CALIN_TYPEALIAS(bool_vt, typename VCLArchitecture::double_bvt);

  CALIN_TYPEALIAS(VCLReal, calin::util::vcl::VCLDoubleReal<VCLArchitecture>);
  CALIN_TYPEALIAS(ArchRNG, calin::math::rng::VCLRNG<VCLArchitecture>);
  CALIN_TYPEALIAS(RNG, calin::math::rng::VCLRealRNG<VCLReal>);
  CALIN_TYPEALIAS(RayTracer, calin::simulation::vcl_raytracer::VCLScopeRayTracer<VCLReal>);

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

  VCLRayTracerRayProcessorDouble(calin::simulation::vs_optics::VSOArray* array,
    SingleRayVCLScopeTraceInfoProcessor* vcl_sti_visitor,
    unsigned vcl_sti_visitor_status_min = 0,
    ArchRNG* rng = nullptr,
    bool adopt_array = false, bool adopt_visitor = false,
    bool adopt_rng = false): VCLRayTracerRayProcessorDouble(array,
      (calin::simulation::pe_processor::PEProcessor*)nullptr, rng, adopt_array,
      adopt_visitor, adopt_rng)
  {
    vcl_sti_visitor_ = vcl_sti_visitor;
    vcl_sti_visitor_status_min_ = vcl_sti_visitor_status_min;
  }

  virtual ~VCLRayTracerRayProcessorDouble()
  {
    if(adopt_array_)delete array_;
    if(adopt_visitor_)delete visitor_;
    if(adopt_visitor_)delete vcl_sti_visitor_;
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
    for(unsigned iscope=0; iscope<array_->numTelescopes(); iscope++) {
      nray_[iscope] = 0;
    }
    if(vcl_sti_visitor_) {
      vcl_sti_visitor_->start_processing();
    }
    if(visitor_) {
      visitor_->start_processing();
    }
  }

  void process_ray(unsigned scope_id, const calin::math::ray::Ray& ray,
    double pe_weight) override
  {
    unsigned iray = scope_id * buffer_depth_ * VCLArchitecture::num_double + nray_[scope_id];
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
    if(nray_[scope_id] == buffer_depth_ * VCLArchitecture::num_double) {
      process_buffered_rays(scope_id);
    }
  }

  void finish_processing() override
  {
    for(unsigned iscope=0; iscope<array_->numTelescopes(); iscope++) {
      process_buffered_rays(iscope);
    }
    if(vcl_sti_visitor_) {
      vcl_sti_visitor_->finish_processing();
    }
    if(visitor_) {
      visitor_->finish_processing();
    }
  }

  uint64_t nhit() const { return nhit_; }

private:
  void process_buffered_rays(unsigned scope_id) {
    using Ray = calin::math::ray::VCLRay<VCLReal>;
    using TraceInfo = calin::simulation::vcl_raytracer::VCLScopeTraceInfo<VCLReal>;
    unsigned nray = nray_[scope_id];
    unsigned iray = scope_id * buffer_depth_ * VCLArchitecture::num_double;
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
      double_at fplane_ux;
      double_at fplane_uz;
      int64_at fp_pixel_id;
      unsigned imask = vcl::to_bits(mask);

      info.fplane_x.store(fp_x);
      info.fplane_z.store(fp_z);
      info.fplane_t.store(fp_t);
      info.fplane_ux.store(fplane_ux);
      info.fplane_uz.store(fplane_uz);
      info.pixel_id.store(fp_pixel_id);

      unsigned mray = std::min(VCLArchitecture::num_double, nray);

      if(vcl_sti_visitor_) {
        SingleRayVCLScopeTraceInfo single_info;

        int64_at status;
        double_at reflec_x;
        double_at reflec_y;
        double_at reflec_z;
        uint64_at pre_reflection_obs_hitmask;
        uint64_at post_reflection_obs_hitmask;
        uint64_at camera_obs_hitmask;
        int64_at mirror_hexid;
        int64_at mirror_id;
        double_at mirror_x;
        double_at mirror_y;
        double_at mirror_z;
        double_at mirror_n_dot_u;
        double_at fplane_uy;
        int64_at pixel_hexid;

        info.status.store(status);
        info.reflec_x.store(reflec_x);
        info.reflec_y.store(reflec_y);
        info.reflec_z.store(reflec_z);
        info.pre_reflection_obs_hitmask.store(pre_reflection_obs_hitmask);
        info.post_reflection_obs_hitmask.store(post_reflection_obs_hitmask);
        info.camera_obs_hitmask.store(camera_obs_hitmask);
        info.mirror_hexid.store(mirror_hexid);
        info.mirror_id.store(mirror_id);
        info.mirror_x.store(mirror_x);
        info.mirror_y.store(mirror_y);
        info.mirror_z.store(mirror_z);
        info.mirror_n_dot_u.store(mirror_n_dot_u);
        info.fplane_uy.store(fplane_uy);
        info.pixel_hexid.store(pixel_hexid);

        for(unsigned jray=0; jray<mray; jray++) {
          if(status[jray] >= vcl_sti_visitor_status_min_) {
            ++nhit_;

            single_info.hit_fplane = imask & (0x1<<jray);
            single_info.ray_weight = ray_w_[iray+jray];
            single_info.scope_id = scope_id;
            single_info.status = status[jray];
            single_info.reflec_x = reflec_x[jray];
            single_info.reflec_y = reflec_y[jray];
            single_info.reflec_z = reflec_z[jray];
            single_info.pre_reflection_obs_hitmask = pre_reflection_obs_hitmask[jray];
            single_info.post_reflection_obs_hitmask = post_reflection_obs_hitmask[jray];
            single_info.camera_obs_hitmask = camera_obs_hitmask[jray];
            single_info.mirror_hexid = mirror_hexid[jray];
            single_info.mirror_id = mirror_id[jray];
            single_info.mirror_x = mirror_x[jray];
            single_info.mirror_y = mirror_y[jray];
            single_info.mirror_z = mirror_z[jray];
            single_info.mirror_n_dot_u = mirror_n_dot_u[jray];
            single_info.fplane_x = fp_x[jray];
            single_info.fplane_z = fp_z[jray];
            single_info.fplane_t = fp_t[jray] + ray_ct_[iray+jray]*math::constants::cgs_1_c;
            single_info.fplane_ux = fplane_ux[jray];
            single_info.fplane_uy = fplane_uy[jray];
            single_info.fplane_uz = fplane_uz[jray];
            single_info.pixel_hexid = pixel_hexid[jray];
            single_info.pixel_id = fp_pixel_id[jray];

            vcl_sti_visitor_->process_vcl_scope_trace_info(single_info);
          }
        }
      }

      if(visitor_) {
        for(unsigned jray=0; jray<mray; jray++) {
          if(imask & (0x1<<jray)) {
            ++nhit_;
            visitor_->process_focal_plane_hit(scope_id, fp_pixel_id[jray],
              double(fp_x[jray]), double(fp_z[jray]),
              double(fplane_ux[jray]), double(fplane_uz[jray]),
              double(fp_t[jray]) + ray_ct_[iray+jray]*math::constants::cgs_1_c,
              ray_w_[iray+jray]);
          }
        }
      }

      iray += VCLArchitecture::num_double;
      nray -= mray;
    }
  }

  calin::simulation::vs_optics::VSOArray* array_ = nullptr;
  bool adopt_array_ = false;
  calin::simulation::pe_processor::PEProcessor* visitor_ = nullptr;
  SingleRayVCLScopeTraceInfoProcessor* vcl_sti_visitor_ = nullptr;
  unsigned vcl_sti_visitor_status_min_ = 0;
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

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCLRayTracerRayProcessorFloat:
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
      double fp_ux[VCLArchitecture::num_float];
      double fp_uz[VCLArchitecture::num_float];

      int32_t fp_pixel_id[VCLArchitecture::num_float];
      unsigned imask = vcl::to_bits(mask);

      store(fp_x, info.fplane_x);
      store(fp_z, info.fplane_z);
      store(fp_t, info.fplane_t);
      store(fp_ux, info.fplane_ux);
      store(fp_uz, info.fplane_uz);
      info.pixel_id.store(fp_pixel_id);

      unsigned mray = std::min(VCLArchitecture::num_float, nray);
      for(unsigned jray=0; jray<mray; jray++) {
        if(imask & (0x1<<jray)) {
          ++nhit_;
          visitor_->process_focal_plane_hit(scope_id, fp_pixel_id[jray],
            fp_x[jray], fp_z[jray], fp_ux[jray], fp_uz[jray],
            fp_t[jray] + ray_ct_[iray+jray]*math::constants::cgs_1_c,
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
