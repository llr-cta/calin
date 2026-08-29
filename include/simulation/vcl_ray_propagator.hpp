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
#include <math/spline_interpolation.hpp>
#include <simulation/vso_array.hpp>
#include <simulation/vcl_raytracer.hpp>
#include <simulation/vcl_panoseti_raytracer.hpp>
#include <simulation/ray_processor.hpp>
#include <simulation/vso_ray_processor.hpp>
#include <util/log.hpp>
#include <util/string.hpp>

#include <simulation/detector_efficiency.hpp>

namespace calin { namespace simulation { namespace vcl_ray_propagator {

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCLRayColorizer
{
public:
#ifndef SWIG
  using double_bvt  = typename VCLArchitecture::double_bvt;
  using double_vt   = typename VCLArchitecture::double_vt;
  using Real_vt     = calin::util::vcl::VCLDoubleReal<VCLArchitecture>;
  using Ray_vt      = typename calin::math::ray::VCLRay<Real_vt>;
  using RealRNG     = calin::math::rng::VCLRealRNG<Real_vt>;
#endif // not defined SWIG
  CALIN_TYPEALIAS(ArchRNG, calin::math::rng::VCLRNG<VCLArchitecture>);

  VCLRayColorizer(ArchRNG* rng = nullptr, bool adopt_rng = false):
    rng_(new RealRNG(rng==nullptr ? new ArchRNG(__PRETTY_FUNCTION__) : rng, rng==nullptr ? true : adopt_rng))
  {
    // nothing to see here
  }

  virtual ~VCLRayColorizer() 
  { 
    delete rng_; // adopt_rng is handled by RealRNG destructor
  }

#ifndef SWIG
  virtual void colorize(Ray_vt& ray, double_bvt ray_mask) = 0;
#endif // not defined SWIG
  virtual std::string banner(const std::string& indent0 = "", const std::string& indentN = "") const = 0;

#ifndef SWIG
protected:
  RealRNG* rng_;
#endif
};

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCLInverseCDFRayColorizer:
  public VCLRayColorizer<VCLArchitecture>
{
public:
#ifndef SWIG
  using double_bvt  = typename VCLArchitecture::double_bvt;
  using double_vt   = typename VCLArchitecture::double_vt;
  using Real_vt     = calin::util::vcl::VCLDoubleReal<VCLArchitecture>;
  using Ray_vt      = typename calin::math::ray::VCLRay<Real_vt>;
  using RealRNG     = calin::math::rng::VCLRealRNG<Real_vt>;
#endif // not defined SWIG
  CALIN_TYPEALIAS(ArchRNG, calin::math::rng::VCLRNG<VCLArchitecture>);

  VCLInverseCDFRayColorizer(const Eigen::VectorXd& color_inv_cdf, ArchRNG* rng = nullptr, bool adopt_rng = false):
    VCLRayColorizer<VCLArchitecture>(rng, adopt_rng),
    color_inv_cdf_(color_inv_cdf)
  {
    // nothing to see here
  }

#ifndef SWIG
  void colorize(Ray_vt& ray, double_bvt ray_mask) override {
    if(!vcl::horizontal_or(ray_mask)) return;
    double_vt energy = this->rng_->from_inverse_cdf(color_inv_cdf_.data(), color_inv_cdf_.size());
    ray.set_energy(vcl::select(ray_mask, energy, ray.energy()));
  }
#endif // not defined SWIG

  std::string banner(const std::string& indent0 = "", const std::string& indentN = "") const override {
    return indent0 + "Inverse CDF Ray Colorizer";
  }

private:
  Eigen::VectorXd color_inv_cdf_;
};

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCLLogitSplineRayColorizer:
  public VCLRayColorizer<VCLArchitecture>
{
public:
#ifndef SWIG
  using double_bvt  = typename VCLArchitecture::double_bvt;
  using double_vt   = typename VCLArchitecture::double_vt;
  using Real_vt     = calin::util::vcl::VCLDoubleReal<VCLArchitecture>;
  using Ray_vt      = typename calin::math::ray::VCLRay<Real_vt>;
  using RealRNG     = calin::math::rng::VCLRealRNG<Real_vt>;
#endif // not defined SWIG
  CALIN_TYPEALIAS(ArchRNG, calin::math::rng::VCLRNG<VCLArchitecture>);

  VCLLogitSplineRayColorizer(const calin::math::spline_interpolation::CubicSpline* spline, bool adopt_spline = false, ArchRNG* rng = nullptr, bool adopt_rng = false):
    VCLRayColorizer<VCLArchitecture>(rng, adopt_rng),
    spline_(adopt_spline ? new calin::math::spline_interpolation::CubicSpline(*spline) : spline),
    adopt_spline_(adopt_spline)
  {
    // nothing to see here
  }

  ~VCLLogitSplineRayColorizer() {
    if(adopt_spline_) delete spline_;
  }

#ifndef SWIG
  void colorize(Ray_vt& ray, double_bvt ray_mask) override {
    if(!vcl::horizontal_or(ray_mask)) return;
    double_vt energy = this->rng_->from_inverse_cdf_logit([this](const double_vt& x) {
      return spline_->template vcl_value<VCLArchitecture>(x);
    });
    ray.set_energy(vcl::select(ray_mask, energy, ray.energy()));
  }
#endif // not defined SWIG

  std::string banner(const std::string& indent0 = "", const std::string& indentN = "") const override {
    return indent0 + "Logit Spline Ray Colorizer";
  }

private:
  const calin::math::spline_interpolation::CubicSpline* spline_;
  bool adopt_spline_;
};

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCLAtmosphericAbsorptionRayColorizer:
  public VCLRayColorizer<VCLArchitecture>
{
public:
#ifndef SWIG
  using double_bvt  = typename VCLArchitecture::double_bvt;
  using double_vt   = typename VCLArchitecture::double_vt;
  using Real_vt     = calin::util::vcl::VCLDoubleReal<VCLArchitecture>;
  using Ray_vt      = typename calin::math::ray::VCLRay<Real_vt>;
  using RealRNG     = calin::math::rng::VCLRealRNG<Real_vt>;
  using double_at   = typename VCLArchitecture::double_at;
#endif // not defined SWIG
  CALIN_TYPEALIAS(ArchRNG, calin::math::rng::VCLRNG<VCLArchitecture>);

  VCLAtmosphericAbsorptionRayColorizer(VCLRayColorizer<VCLArchitecture>* base_colorizer, 
    const calin::simulation::detector_efficiency::AtmosphericAbsorption* atm_abs, 
    double zobs, bool adopt_colorizer = false, ArchRNG* rng = nullptr, bool adopt_rng = false):
    VCLRayColorizer<VCLArchitecture>(rng, adopt_rng),
    base_colorizer_(base_colorizer), atm_abs_(atm_abs), zobs_(zobs), adopt_colorizer_(adopt_colorizer)
  {
    // nothing to see here
  }

  ~VCLAtmosphericAbsorptionRayColorizer() {
    if(adopt_colorizer_) delete base_colorizer_;
  }

#ifndef SWIG
  void colorize(Ray_vt& ray, double_bvt ray_mask) override {
    unsigned active_bits = vcl::to_bits(ray_mask);

    double_at ray_energy_at;
    ray.energy().store(ray_energy_at);

    double_at ray_z_at;
    ray.z().store(ray_z_at);

    double_vt ray_uz_inv = 1.0/abs(ray.uz());
    double_at ray_uz_inv_at;
    ray_uz_inv.store(ray_uz_inv_at);
    
    unsigned iray = 0;
    while(iray<VCLArchitecture::num_double) {
      // Propose energies for all lanes
      base_colorizer_->colorize(ray, true);

      double_at proposed_energy_at;
      ray.energy().store(proposed_energy_at);

      double_vt rand_vt = this->rng_->uniform_double();
      double_at rand_at;
      rand_vt.store(rand_at);

      unsigned ilane = 0;
      while(iray<VCLArchitecture::num_double and ilane<VCLArchitecture::num_double) {
        if(not((active_bits >> iray) & 1)) {
          ++iray;
          continue;
        }

        double tau_vertical = atm_abs_->optical_depth_for_altitude_and_energy(ray_z_at[iray], proposed_energy_at[ilane]) 
                            - atm_abs_->optical_depth_for_altitude_and_energy(zobs_, proposed_energy_at[ilane]);
        double prob = std::exp(-std::max(0.0, tau_vertical * ray_uz_inv_at[iray]));
        double rand = rand_at[ilane];

        if(rand < prob) {
          // Color is accepted, set the ray energy as proposed
          ray_energy_at[iray] = proposed_energy_at[ilane];
          ++iray;
        } else {
          // Color is rejected, go again with next lane
        }

        ++ilane;
      }
    }

    ray.mutable_energy().load(ray_energy_at);
  }
#endif // not defined SWIG

  std::string banner(const std::string& indent0 = "", const std::string& indentN = "") const override {
    return indent0 + "Atmospheric Absorption Ray Colorizer\n" + base_colorizer_->banner(indent0 + "  ", indentN + "  ");
  }

private:
  VCLRayColorizer<VCLArchitecture>* base_colorizer_;
  const calin::simulation::detector_efficiency::AtmosphericAbsorption* atm_abs_;
  double zobs_;
  bool adopt_colorizer_;
};

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
  double_vt detection_prob;
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

  virtual void point_telescope_az_el_phi_deg(unsigned iscope,
      double az_deg, double el_deg, double phi_deg=0.0) {
    // nothing to see here
  }

  virtual std::string banner(const std::string& indent0 = "", const std::string& indentN = "") const {
    return {};
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
      ArchRNG* rng = nullptr, double ref_index = 1.0,
      bool adopt_array = false, bool adopt_rng = false):
    array_(adopt_array ? new calin::simulation::vs_optics::VSOArray(*array) : array),
    rng_(new RealRNG(rng==nullptr ? new ArchRNG(__PRETTY_FUNCTION__) : rng,
      rng==nullptr ? true : adopt_rng)),
    ray_tracer_(array->numTelescopes()), ref_index_(ref_index)
  {
    for(unsigned iscope=0;iscope<array->numTelescopes();++iscope) {
      ray_tracer_[iscope] = new RayTracer(
        array->telescope(iscope), ref_index_, rng_, /* adopt_rng= */ false);
    }
  }

  virtual ~DaviesCottonVCLFocalPlaneRayPropagator() {
    for(auto* rt : ray_tracer_)delete rt;
    delete array_;
  }

  virtual std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere> detector_spheres() {
    return calin::simulation::vso_ray_processor::VSORayProcessor::detector_spheres_for_array(array_);
  }

  void start_propagating() final {
    // nothing to see here
  }

  void point_telescope_az_el_phi_deg(unsigned iscope,
      double az_deg, double el_deg, double phi_deg) final {
    if(iscope >= array_->numTelescopes()) {
      throw std::out_of_range("Telescope number out of range");
    }
    array_->telescope(iscope)->pointTelescopeAzElPhi(
      az_deg/180.0*M_PI, el_deg/180.0*M_PI, phi_deg/180.0*M_PI);
    ray_tracer_[iscope]->point_telescope(array_->telescope(iscope));
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
    fp_parameters.detection_prob = 1.0;
    return ray_mask;
  }
#endif

  void finish_propagating() final {
    // nothing to see here
  }

  std::string banner(const std::string& indent0 = "", const std::string& indentN = "") const final {
    return array_->banner(indent0, indentN);
  }

  const calin::simulation::vs_optics::VSOArray* array() const { return array_; }

#ifndef SWIG
private:
  calin::simulation::vs_optics::VSOArray* array_ = nullptr;
  unsigned vcl_sti_visitor_status_min_ = 0;
  RealRNG* rng_ = nullptr;
  std::vector<RayTracer*> ray_tracer_;
  double ref_index_ = 1.0;
#endif
};

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) PerfectOpticsVCLFocalPlaneRayPropagator:
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

  using TraceInfo   = calin::simulation::vcl_raytracer::VCLScopeTraceInfo<Real_vt>;
#endif // not defined SWIG

  CALIN_TYPEALIAS(ArchRNG, calin::math::rng::VCLRNG<VCLArchitecture>);

  PerfectOpticsVCLFocalPlaneRayPropagator(double ref_index = 1.0):
    ref_index_(ref_index)
  {
    // nothing to see here
  }

  virtual ~PerfectOpticsVCLFocalPlaneRayPropagator() {
    for(auto* scope : scopes_) {
      delete scope;
    }
  }

  void add_telescope(const Eigen::Vector3d& r0, double radius, unsigned iobs, 
      double focal_length = 0, double field_of_view_radius = M_PI/2) {
    TelescopeDetails* scope = new TelescopeDetails;
    scope->r0                   = r0;
    scope->radius               = radius;
    scope->radius_squared       = radius*radius;
    scope->iobs                 = iobs;
    scope->obs_dir              = Eigen::Vector3d::UnitY();
    scope->field_of_view_radius = field_of_view_radius;
    scope->field_of_view_uycut  = std::cos(field_of_view_radius);
    scope->global_to_fp         = Eigen::Matrix3d::Identity();
    scope->focal_length         = focal_length;
    scopes_.push_back(scope);
  }

  virtual std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere> detector_spheres() {
    std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere> spheres;
    spheres.reserve(scopes_.size());
    for(const auto* scope : scopes_) {
      spheres.emplace_back(scope->r0, scope->radius, scope->obs_dir, scope->field_of_view_radius,
        scope->iobs);
    }
    return spheres;
  }

  void start_propagating() final {
    // nothing to see here
  }

  void point_telescope_az_el_phi_deg(unsigned iscope,
      double az_deg, double el_deg, double phi_deg) final {
    if(iscope >= scopes_.size()) {
      throw std::out_of_range("Telescope number out of range");
    }
    auto* scope = scopes_[iscope];
    scope->global_to_fp = 
      Eigen::AngleAxisd(-phi_deg*M_PI/180.0, Eigen::Vector3d::UnitY()) *
      Eigen::AngleAxisd(-el_deg*M_PI/180.0, Eigen::Vector3d::UnitX()) *
      Eigen::AngleAxisd(az_deg*M_PI/180.0, Eigen::Vector3d::UnitZ());
    scope->obs_dir = scope->global_to_fp.transpose() * Eigen::Vector3d::UnitY();
  }

#ifndef SWIG
  double_bvt propagate_rays_to_focal_plane(
      unsigned scope_id, Ray_vt& ray, double_bvt ray_mask,
      VCLFocalPlaneParameters<VCLArchitecture>& fp_parameters) final {
    if(scope_id >= scopes_.size()) {
      throw std::out_of_range("Telescope number out of range");
    }
    const auto* scope = scopes_[scope_id];
    ray.translate_origin(scope->r0.template cast<double_vt>());
    ray.rotate(scope->global_to_fp.template cast<double_vt>());
    double_vt uy = -ray.uy();
    ray_mask = ray.propagate_to_y_plane_with_mask(ray_mask, /* d= */ 0, /* time_reversal_ok = */ false, ref_index_);
    ray_mask &= ray.x()*ray.x() + ray.z()*ray.z() < scope->radius_squared;
    ray_mask &= uy >= scope->field_of_view_uycut;
    double_vt F_over_uy = scope->focal_length/uy;
    TraceInfo info;
    fp_parameters.fplane_x       = select(ray_mask, F_over_uy * ray.ux(), 0);
    fp_parameters.fplane_y       = 0.0;
    fp_parameters.fplane_z       = select(ray_mask, F_over_uy * ray.uz(), 0);
    fp_parameters.fplane_ux      = select(ray_mask, ray.ux(), 0);
    fp_parameters.fplane_uy      = select(ray_mask, uy, 0);
    fp_parameters.fplane_uz      = select(ray_mask, ray.uz(), 0);
    fp_parameters.fplane_t       = select(ray_mask, ray.time(), 0);
    fp_parameters.pixel_id       = 0;
    fp_parameters.detection_prob = 1.0;
    return ray_mask;
  }
#endif

  void finish_propagating() final {
    // nothing to see here
  }

  std::string banner(const std::string& indent0 = "", const std::string& indentN = "") const final {
    using calin::util::string::double_to_string_with_commas;
    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymin = std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();
    for(const auto* scope : scopes_) {
      xmin = std::min(xmin, scope->r0[0]);
      xmax = std::max(xmax, scope->r0[0]);
      ymin = std::min(ymin, scope->r0[1]);
      ymax = std::max(ymax, scope->r0[1]);
    }
    double A = (xmax-xmin)*(ymax-ymin)*1e-10;
    std::ostringstream stream;
    stream << indent0 << "Array of " << this->scopes_.size()
      << " pseudo-telescopes covering "
      << double_to_string_with_commas(A,3) << " km^2.";
    return stream.str();
  }

#ifndef SWIG
private:

  struct TelescopeDetails
  {
    Eigen::Vector3d r0;            // Center of detector sphere [cm]
    double radius;                 // Radius of sphere [cm]
    double radius_squared;         // Squared radius of sphere [cm^2]
    unsigned iobs;                 // Observation layer associated with this detector
    Eigen::Vector3d obs_dir;       // Pointing direction of detector
    double field_of_view_radius;   // Field of view of detector [radians]
    double field_of_view_uycut;    // Value of cos(uy) that coresponds to FoV radius
    Eigen::Matrix3d global_to_fp;  // Rotation matrix from global to focal plane
    double focal_length;           // Focal length [cm]
  };

  double ref_index_ = 1.0;
  std::vector<TelescopeDetails*> scopes_;
#endif
};

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) AllSkyVCLFocalPlaneRayPropagator:
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

  using TraceInfo   = calin::simulation::vcl_raytracer::VCLScopeTraceInfo<Real_vt>;
#endif // not defined SWIG

  CALIN_TYPEALIAS(ArchRNG, calin::math::rng::VCLRNG<VCLArchitecture>);

  AllSkyVCLFocalPlaneRayPropagator(unsigned observation_level, 
      const Eigen::VectorXd& sphere_center, double sphere_radius, double field_of_view_radius = M_PI/2,
      double ref_index = 1.0):
    observation_level_(observation_level), sphere_center_(sphere_center), 
    sphere_radius_(sphere_radius), sphere_radius_squared_(sphere_radius * sphere_radius),
    sphere_field_of_view_radius_(field_of_view_radius), 
    sphere_field_of_view_uycut_(-std::cos(field_of_view_radius)),
    global_to_fp_(Eigen::Matrix3d::Identity()), ref_index_(ref_index)
  {
    // nothing to see here
  }

  virtual ~AllSkyVCLFocalPlaneRayPropagator() {
    // nothing to see here
  }

  virtual std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere> detector_spheres() {
    std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere> spheres;
    spheres.emplace_back(sphere_center_, sphere_radius_, 
      global_to_fp_.transpose() * Eigen::Vector3d::UnitY(), M_PI/2, observation_level_);
    return spheres;
  }

  void start_propagating() final {
    // nothing to see here
  }

  void point_telescope_az_el_phi_deg(unsigned iscope,
      double az_deg, double el_deg, double phi_deg) final {
    if(iscope >= 1) {
      throw std::out_of_range("Telescope number out of range");
    }
    global_to_fp_ = 
      Eigen::AngleAxisd(phi_deg*M_PI/180.0, Eigen::Vector3d::UnitZ()) *
      Eigen::AngleAxisd(-el_deg*M_PI/180.0, Eigen::Vector3d::UnitX()) *
      Eigen::AngleAxisd(az_deg*M_PI/180.0, Eigen::Vector3d::UnitZ());
  }

#ifndef SWIG
  double_bvt propagate_rays_to_focal_plane(
      unsigned scope_id, Ray_vt& ray, double_bvt ray_mask,
      VCLFocalPlaneParameters<VCLArchitecture>& fp_parameters) final {
    if(scope_id >= 1) {
      throw std::out_of_range("Telescope number out of range");
    }
    ray.translate_origin(sphere_center_.template cast<double_vt>());
    ray_mask = ray.propagate_to_z_plane_with_mask(ray_mask, /* d= */ 0, 
      /* time_reversal_ok = */ false, ref_index_);
    double_vt x = ray.x();
    double_vt y = ray.y();
    ray_mask &= x*x+y*y < sphere_radius_squared_;
    ray.set_direction(global_to_fp_.template cast<double_vt>() * ray.direction());
    ray_mask &= ray.uy() <= sphere_field_of_view_uycut_;
    TraceInfo info;
    fp_parameters.fplane_x       = select(ray_mask, x, 0);
    fp_parameters.fplane_y       = select(ray_mask, 0, 0); // ** UNUSED **
    fp_parameters.fplane_z       = select(ray_mask, y, 0);
    fp_parameters.fplane_ux      = select(ray_mask, ray.ux(), 0);
    fp_parameters.fplane_uy      = select(ray_mask, ray.uy(), 0); // ** UNUSED **
    fp_parameters.fplane_uz      = select(ray_mask, ray.uz(), 0);
    fp_parameters.fplane_t       = select(ray_mask, ray.time(), 0);
    fp_parameters.pixel_id       = 0;
    fp_parameters.detection_prob = 1.0;
    return ray_mask;
  }
#endif

  void finish_propagating() final {
    // nothing to see here
  }

  std::string banner(const std::string& indent0 = "", const std::string& indentN = "") const final {
    std::ostringstream stream;
    stream << indent0 << "All-sky detector.";
    return stream.str();
  }

#ifndef SWIG
private:
  unsigned observation_level_;          // Observation layer associated with this detector
  Eigen::Vector3d sphere_center_;       // Center of detector sphere [cm]
  double sphere_radius_;                // Radius of sphere [cm]
  double sphere_radius_squared_;        // Squared radius of sphere [cm^2]
  double sphere_field_of_view_radius_;  // Field of view of detector [radians]
  double sphere_field_of_view_uycut_;   // Value of cos(uy) that coresponds to FoV radius
  Eigen::Matrix3d global_to_fp_;        // Rotation matrix from global to focal plane

  double ref_index_ = 1.0;
#endif
};

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) PanosetiVCLFocalPlaneRayPropagator:
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

  using RayTracer   = calin::simulation::vcl_raytracer::VCLPanosetiThinLensScopeRayTracer<Real_vt>;
  using TraceInfo   = calin::simulation::vcl_raytracer::VCLPanosetiScopeTraceInfo<Real_vt>;
  using RealRNG     = calin::math::rng::VCLRealRNG<Real_vt>;
#endif // not defined SWIG

  CALIN_TYPEALIAS(ArchRNG, calin::math::rng::VCLRNG<VCLArchitecture>);

  PanosetiVCLFocalPlaneRayPropagator(const calin::ix::simulation::panoseti_optics::ArrayParameters& array_params,
      const calin::math::spline_interpolation::CubicSpline* lens_refractive_index_spline,
      VCLRayColorizer<VCLArchitecture>* colorizer = nullptr,
      ArchRNG* rng = nullptr, double air_refractive_index = 1.0, unsigned observation_layer = 0,
      bool adopt_lens_refractive_index_spline = false, bool adopt_colorizer = false, bool adopt_rng = false):
    array_params_(array_params),
    lens_refractive_index_spline_(adopt_lens_refractive_index_spline ? new calin::math::spline_interpolation::CubicSpline(*lens_refractive_index_spline) : lens_refractive_index_spline),
    adopt_lens_refractive_index_spline_(adopt_lens_refractive_index_spline),
    colorizer_(colorizer), adopt_colorizer_(adopt_colorizer),
    rng_(new RealRNG(rng==nullptr ? new ArchRNG(__PRETTY_FUNCTION__) : rng, rng==nullptr ? true : adopt_rng)),
    ray_tracer_(array_params.scope_positions_size()), air_ref_index_(air_refractive_index)
  {
    for(unsigned iscope=0;iscope<unsigned(array_params.scope_positions_size());++iscope) {
      ray_tracer_[iscope] = new RayTracer(
        array_params_, iscope, lens_refractive_index_spline_,
        air_ref_index_, observation_layer, rng_, /* adopt_lens_refractive_index_spline= */ false, /* adopt_rng= */ false);
    }
  }

  virtual ~PanosetiVCLFocalPlaneRayPropagator() {
    for(auto* rt : ray_tracer_)delete rt;
    if(adopt_lens_refractive_index_spline_) delete lens_refractive_index_spline_;
    if(adopt_colorizer_) delete colorizer_;
    delete rng_; // adopt_rng_ is handled by RealRNG
  }

  virtual std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere> detector_spheres() {
    std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere> spheres;
    spheres.reserve(ray_tracer_.size());
    for(const auto* rt : ray_tracer_) {
      spheres.push_back(rt->detector_sphere());
    }
    return spheres;
  }

  void start_propagating() final {
    // nothing to see here
  }

  void point_telescope_az_el_phi_deg(unsigned iscope,
      double az_deg, double el_deg, double phi_deg) final {
    if(iscope >= ray_tracer_.size()) {
      throw std::out_of_range("Telescope number out of range");
    }
    ray_tracer_[iscope]->point_telescope_az_el_phi(
      az_deg/180.0*M_PI, el_deg/180.0*M_PI, phi_deg/180.0*M_PI);
  }

#ifndef SWIG
  double_bvt propagate_rays_to_focal_plane(
      unsigned scope_id, Ray_vt& ray, double_bvt ray_mask,
      VCLFocalPlaneParameters<VCLArchitecture>& fp_parameters) final {
    colorizer_->colorize(ray, ray_mask);

    TraceInfo info;
    ray_mask = ray_tracer_[scope_id]->
      trace_global_frame(ray_mask, ray, info, /* do_derotation = */ false);
    fp_parameters.fplane_x     = info.detector_x;
    fp_parameters.fplane_y     = 0.0;
    fp_parameters.fplane_z     = info.detector_z;
    fp_parameters.fplane_ux    = info.detector_ux;
    fp_parameters.fplane_uy    = info.detector_uy;
    fp_parameters.fplane_uz    = info.detector_uz;
    fp_parameters.fplane_t     = info.detector_t;
    fp_parameters.pixel_id     = info.pixel_id.template cast<int64_vt>();
    fp_parameters.detection_prob = 1.0;
    return ray_mask;
  }
#endif

  void finish_propagating() final {
    // nothing to see here
  }

  std::string banner(const std::string& indent0 = "", const std::string& indentN = "") const final {
    std::ostringstream stream;
    stream << indent0 << "Array of " << this->ray_tracer_.size() << " PANOSETI telescopes.";
    return stream.str();
  }

#ifndef SWIG
private:
  calin::ix::simulation::panoseti_optics::ArrayParameters array_params_;
  const calin::math::spline_interpolation::CubicSpline* lens_refractive_index_spline_ = nullptr;
  bool adopt_lens_refractive_index_spline_ = false;
  VCLRayColorizer<VCLArchitecture>* colorizer_ = nullptr;
  bool adopt_colorizer_ = false;
  RealRNG* rng_ = nullptr;
  std::vector<RayTracer*> ray_tracer_;
  double air_ref_index_ = 1.0;
#endif
};

} } } // namespace calin::simulations::vcl_ray_propagator
