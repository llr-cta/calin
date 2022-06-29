/*

   calin/simulation/vcl_ray_processor.hpp -- Stephen Fegan -- 2022-06-29

   Class for propagating rays through VCL raytracers

   Copyright 2022, Stephen Fegan <sfegan@llr.in2p3.fr>
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

namespace calin { namespace simulation { namespace vcl_ray_propagator {

template<typename VCLArchitecture> struct alignas(VCLArchitecture::vec_bytes) VCLFocalPlaneParameters
{
#ifndef SWIG
  using int64_vt    = typename VCLArchitecture::int64_vt;
  using double_bvt  = typename VCLArchitecture::double_bvt;
  using double_vt   = typename VCLArchitecture::double_vt;
  using Vector3d_vt = typename VCLArchitecture::Vector3d_vt;
  using Ray_vt      = typename calin::math::ray::VCLRay<VCLArchitecture>;

  double_vt fplane_x;
  double_vt fplane_y;
  double_vt fplane_z;
  double_vt fplane_ux;
  double_vt fplane_uy;
  double_vt fplane_uz;
  double_vt fplane_t;
  int64_vt pixel_id;
#endif // not defined SWIG
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

#ifndef SWIG
  virtual double_bvt propagate_rays_to_focal_plane(
      unsigned scope_id, Ray_vt& ray, double_bvt ray_mask,
      VCLFocalPlaneParameters<VCLArchitecture>& fp_parameters) {
    return false;
  }
#endif

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

  using RayTracer   = calin::simulation::vcl_raytracer::VCLScopeRayTracer<Real_vt>;
  using TraceInfo   = calin::simulation::vcl_raytracer::VCLScopeTraceInfo<Real_vt>;
  using RealRNG     = calin::math::rng::VCLRealRNG<Real_vt>;
#endif // not defined SWIG

  CALIN_TYPEALIAS(ArchRNG, calin::math::rng::VCLRNG<VCLArchitecture>);

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

  void start_propagating() final {
    // nothing to see here
  }

#ifndef SWIG
  double_bvt propagate_rays_to_focal_plane(
      unsigned scope_id, Ray_vt& ray, double_bvt ray_mask,
      VCLFocalPlaneParameters<VCLArchitecture>& fp_parameters) final {
    TraceInfo info;
    ray_mask = ray_tracer_[scope_id]->
      trace_global_frame(ray_mask, ray, info, /* do_derotation = */ false);
    fp_parameters.fplane_x     = info.fplane_x;
    fp_parameters.fplane_y     = 0.0; // info.fplane_y;
    fp_parameters.fplane_z     = info.fplane_z;
    fp_parameters.fplane_ux    = info.fplane_ux;
    fp_parameters.fplane_uy    = info.fplane_uy;
    fp_parameters.fplane_uz    = info.fplane_uz;
    fp_parameters.fplane_t     = info.fplane_t;
    fp_parameters.pixel_id     = info.pixel_id;
    return ray_mask;
  }
#endif

  void finish_propagating() final {
    // nothing to see here
  }

#ifndef SWIG
private:
  calin::simulation::vs_optics::VSOArray* array_ = nullptr;
  bool adopt_array_ = false;
  unsigned vcl_sti_visitor_status_min_ = 0;
  RealRNG* rng_ = nullptr;
  bool adopt_rng_ = false;
  std::vector<RayTracer*> ray_tracer_;
  double ref_index_ = 1.0;
#endif
};

} } } // namespace calin::simulations::vcl_ray_propagator
