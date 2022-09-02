/*

   calin/simulation/vcl_iact.hpp -- Stephen Fegan -- 2019-02-19

   Class for imaging atmospheric cherenkov technique - produce rays from charged
   tracks in the atmosphere, propagate them to the ground and trace them through
   optics.

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

#include <util/vcl.hpp>
#include <math/rng.hpp>
#include <math/rng_vcl.hpp>
#include <math/special.hpp>
#include <math/geometry_vcl.hpp>
#include <math/ray_vcl.hpp>
#include <simulation/tracker.hpp>
#include <simulation/atmosphere.hpp>
#include <util/log.hpp>
#include <simulation/air_cherenkov_tracker.hpp>
#include <simulation/vcl_iact.pb.h>

namespace calin { namespace simulation { namespace vcl_iact {

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCLIACTTrackVisitor:
  public calin::simulation::tracker::TrackVisitor
{
public:
#ifndef SWIG
  using uint32_vt  = typename VCLArchitecture::uint32_vt;
  using int32_vt   = typename VCLArchitecture::int32_vt;
  using uint64_vt  = typename VCLArchitecture::uint64_vt;
  using int64_vt   = typename VCLArchitecture::int64_vt;
  using float_vt   = typename VCLArchitecture::float_vt;
  using double_vt  = typename VCLArchitecture::double_vt;
  using uint64_bvt  = typename VCLArchitecture::uint64_bvt;
  using float_bvt   = typename VCLArchitecture::float_bvt;
  using double_bvt  = typename VCLArchitecture::double_bvt;
  using Vector3f_vt = typename VCLArchitecture::Vector3f_vt;
  using Vector3d_vt = typename VCLArchitecture::Vector3d_vt;
  using float_real  = typename VCLArchitecture::float_real;
  using double_real = typename VCLArchitecture::double_real;
  using float_at    = typename VCLArchitecture::float_at;
  using double_at   = typename VCLArchitecture::double_at;
#endif // not defined SWIG

  VCLIACTTrackVisitor(calin::simulation::atmosphere::LayeredRefractiveAtmosphere* atm,
    const calin::ix::simulation::vcl_iact::VCLIACTConfiguration& config = default_config(),
    calin::math::rng::VCLRNG<VCLArchitecture>* rng = nullptr,
    bool adopt_atm = false, bool adopt_rng = false);
  virtual ~VCLIACTTrackVisitor();
  void visit_event(const calin::simulation::tracker::Event& event, bool& kill_event) override;
  void visit_track(const calin::simulation::tracker::Track& track, bool& kill_track) override;
  void leave_event() override;

  uint64_t num_tracks() const { return num_tracks_; }
  uint64_t num_steps() const { return num_steps_; }
  uint64_t num_rays() const { return num_rays_; }

  double total_yield() const { return vcl::horizontal_add(sum_yield_); }
  double total_yield_by_height() const { return vcl::horizontal_add(sum_yield_h_); }
  double total_yield_by_height_squared() const { return vcl::horizontal_add(sum_yield_h2_); }

  static calin::ix::simulation::vcl_iact::VCLIACTConfiguration default_config() {
    calin::ix::simulation::vcl_iact::VCLIACTConfiguration config;
    config.set_bandwidth(3.0);
    config.set_enable_forced_cherenkov_angle_mode(false);
    config.set_forced_cherenkov_angle(-1.0);
    return config;
  }

#ifndef SWIG
  virtual void propagate_rays(calin::math::ray::VCLRay<double_real> ray, double_bvt ray_mask,
    double_vt bandwidth, double_vt weight);

  void set_fixed_pe_bandwidth_mode(double bandwidth);
  void set_fixed_photon_bandwidth_mode(double bandwidth, double min_cherenkov_energy);
  void set_height_dependent_pe_bandwidth_mode(calin::math::spline_interpolation::CubicSpline* spline,
    bool adopt_spline = false);

protected:
  inline int insert_track(const Eigen::Vector3d& x, const double t, const double g,
    const double yield, const Eigen::Vector3d& u, const double dx, const double dt_dx,
    const double dg_dx, const double weight, bool valid = true);
  void generate_mc_rays(bool drain_tracks = false);

  calin::simulation::atmosphere::LayeredRefractiveAtmosphere* atm_ = nullptr;
  bool adopt_atm_ = false;
  calin::math::rng::VCLRNG<VCLArchitecture>* rng_ = nullptr;
  bool adopt_rng_ = false;

  double_bvt track_valid_; // track : is this a valid track or not

  Vector3d_vt track_x_;   // track : position at track point
  double_vt track_t_;     // track : time at track point
  double_vt track_g_;     // track : gamma at track point
  double_vt track_yield_const_; // track : cherenkov yield constant for this track

  Vector3d_vt track_u_;   // track : normalized direction vector
  double_vt track_dx_;    // track : distance reamaining to end
  double_vt track_dt_dx_; // track : rate of change of time per unit track length
  double_vt track_dg_dx_; // track : rate of change of gamma of particle along track

  double_vt track_weight_; // track : weight for thinning

  double_vt sum_yield_;
  double_vt sum_yield_h_;
  double_vt sum_yield_h2_;

  uint64_t num_tracks_ = 0;
  uint64_t num_steps_ = 0;
  uint64_t num_rays_ = 0;

  double min_cherenkov_energy_ = 1.5;
  double fixed_bandwidth_ = 3.0;
  bool do_color_photons_ = true;
  calin::math::spline_interpolation::CubicSpline* variable_bandwidth_spline_ = nullptr;
  double adopt_variable_bandwidth_spline_ = false;
  double forced_sin2theta_ = -1.0;
  double cherenkov_weight_ = 1.0;
#endif // not defined SWIG
};

#ifndef SWIG

template<typename VCLArchitecture> VCLIACTTrackVisitor<VCLArchitecture>::
VCLIACTTrackVisitor(
    calin::simulation::atmosphere::LayeredRefractiveAtmosphere* atm,
    const calin::ix::simulation::vcl_iact::VCLIACTConfiguration& config,
    calin::math::rng::VCLRNG<VCLArchitecture>* rng,
    bool adopt_atm, bool adopt_rng):
  calin::simulation::tracker::TrackVisitor(),
  atm_(atm), adopt_atm_(adopt_atm),
  rng_(rng ? rng : new calin::math::rng::VCLRNG<VCLArchitecture>()),
  adopt_rng_(rng ? adopt_rng : true),
  fixed_bandwidth_(config.bandwidth()),
  forced_sin2theta_(config.enable_forced_cherenkov_angle_mode()?
    calin::math::special::SQR(std::sin(config.forced_cherenkov_angle()/180.0*M_PI)) : -1.0)
{
  // nothing to see here
}

template<typename VCLArchitecture> VCLIACTTrackVisitor<VCLArchitecture>::
~VCLIACTTrackVisitor()
{
  // do some tests that there are not unprocessed tracks or rays hanging around
  if(adopt_rng_)delete rng_;
}

template<typename VCLArchitecture> void VCLIACTTrackVisitor<VCLArchitecture>::
set_fixed_pe_bandwidth_mode(double bandwidth)
{
  delete(variable_bandwidth_spline_);
  variable_bandwidth_spline_ = nullptr;
  fixed_bandwidth_ = bandwidth;
  do_color_photons_ = false;
}

template<typename VCLArchitecture> void VCLIACTTrackVisitor<VCLArchitecture>::
set_fixed_photon_bandwidth_mode(double bandwidth, double min_cherenkov_energy)
{
  delete(variable_bandwidth_spline_);
  variable_bandwidth_spline_ = nullptr;
  fixed_bandwidth_ = bandwidth;
  min_cherenkov_energy_ = min_cherenkov_energy;
  do_color_photons_ = true;
}

template<typename VCLArchitecture> void VCLIACTTrackVisitor<VCLArchitecture>::
set_height_dependent_pe_bandwidth_mode(calin::math::spline_interpolation::CubicSpline* spline,
  bool adopt_spline)
{
  delete(variable_bandwidth_spline_);
  if(adopt_spline) {
    variable_bandwidth_spline_ = spline;
  } else {
    variable_bandwidth_spline_ = new calin::math::spline_interpolation::CubicSpline(*spline);
  }
  fixed_bandwidth_ = 0;
  min_cherenkov_energy_ = 0;
  do_color_photons_ = false;
}

template<typename VCLArchitecture> void VCLIACTTrackVisitor<VCLArchitecture>::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
  num_tracks_ = 0;
  num_steps_ = 0;
  num_rays_ = 0;
  track_dx_ = 0.0;
  track_valid_ = false;
  sum_yield_ = 0;
  sum_yield_h_ = 0;
  sum_yield_h2_ = 0;
}

template<typename VCLArchitecture> void VCLIACTTrackVisitor<VCLArchitecture>::
visit_track(const calin::simulation::tracker::Track& track, bool& kill_track)
{
  using namespace calin::util::log;

  double zobs = atm_->zobs(0);
  double dx = track.dx;

  bool z0_at_or_below_ground = track.x0.z() <= zobs;
  bool z1_below_ground = track.x1.z() < zobs;

  kill_track = z0_at_or_below_ground;

  // it's store policy: no charge, no track, no photons
  if(track.q == 0 or dx <= 0 or z0_at_or_below_ground)return;

  if(z1_below_ground) {
    dx = (zobs - track.x0.z())/track.dx_hat.z();
    if(dx <= 0)return;
  }

  using calin::math::special::SQR;
  using calin::simulation::air_cherenkov_tracker::YIELD_CONST;

  // Watch out for gamma^2 slightly less than 1.0
  vcl::Vec2d g = vcl::max(vcl::Vec2d(track.e0, track.e1)/track.mass, 1.0);

  double dx_inv = 1.0/dx;

  ++num_tracks_;

  int unused_track_count =
    insert_track(track.x0, track.t0, g[0], YIELD_CONST*SQR(track.q),
      track.dx_hat, dx, track.dt*dx_inv, (g[1]-g[0])*dx_inv, track.weight);

  if(unused_track_count == 1)
  {
    // All vector slots are full now
    generate_mc_rays();
  }
}

template<typename VCLArchitecture> inline int VCLIACTTrackVisitor<VCLArchitecture>::
insert_track(const Eigen::Vector3d& x, const double t, const double g,
  const double yield, const Eigen::Vector3d& u, const double dx, const double dt_dx,
  const double dg_dx, const double weight, bool valid)
{
  using calin::util::vcl::insert_into_vec3_with_mask;
  using calin::util::vcl::insert_into_with_mask;
  using namespace calin::util::log;

  double_bvt unused_tracks = track_dx_ <= 0;
  int insert_index = vcl::horizontal_find_first(unused_tracks);

  if(insert_index == -1) {
    throw std::logic_error(
      calin::util::vcl::templated_class_name<VCLArchitecture>("VCLIACTTrackVisitor")
      + "::insert_track: No free SIMD vector slot");
  }

  double_bvt insert_mask = false;
  insert_mask.insert(insert_index, true);

  if(valid)track_valid_ = track_valid_ || insert_mask;
  else track_valid_ = vcl::andnot(track_valid_, insert_mask);

  insert_into_vec3_with_mask<double_real>(track_x_, x, insert_mask);
  insert_into_with_mask<double_real>(track_t_, t, insert_mask);
  insert_into_with_mask<double_real>(track_g_, g, insert_mask);
  insert_into_with_mask<double_real>(track_yield_const_, yield, insert_mask);

  insert_into_vec3_with_mask<double_real>(track_u_, u, insert_mask);
  insert_into_with_mask<double_real>(track_dx_, dx, insert_mask);
  insert_into_with_mask<double_real>(track_dt_dx_, dt_dx, insert_mask);
  insert_into_with_mask<double_real>(track_dg_dx_, dg_dx, insert_mask);

  insert_into_with_mask<double_real>(track_weight_, weight, insert_mask);

  return vcl::horizontal_count(unused_tracks);
}


template<typename VCLArchitecture> void VCLIACTTrackVisitor<VCLArchitecture>::
leave_event()
{
  while(horizontal_or(track_valid_)) {
    int index = horizontal_find_first(track_dx_ > 0);

    if(index == -1)
      throw std::logic_error(
        calin::util::vcl::templated_class_name<VCLArchitecture>("VCLIACTTrackVisitor")
        + "::leave_event: No active track found");

    Eigen::Vector3d x(track_x_.x()[index], track_x_.y()[index], track_x_.z()[index]);
    Eigen::Vector3d u(track_u_.x()[index], track_u_.y()[index], track_u_.z()[index]);

    while(horizontal_or(track_dx_ <= 0)) {
      insert_track(x, track_t_[index], track_g_[index], track_yield_const_[index],
        u, track_dx_[index], track_dg_dx_[index], track_dg_dx_[index], false);
    }

    generate_mc_rays(/* drain_tracks = */ true);
  }
}

template<typename VCLArchitecture> void VCLIACTTrackVisitor<VCLArchitecture>::
generate_mc_rays(bool drain_tracks)
{
  using calin::math::special::SQR;
  using namespace calin::util::log;

  do {
    ++num_steps_;

    double_vt sin2thetac; // sin squared of the Cherenkov angle

    if(forced_sin2theta_ > 0) {
      sin2thetac = forced_sin2theta_;
    } else {
      double_vt nmo = atm_->vcl_n_minus_one<VCLArchitecture>(track_x_.z());
      double_vt n2inv = 1.0/SQR(nmo + 1.0);
      double_vt g2 = SQR(track_g_);
      double_vt b2inv = g2/(g2-1.0);
      sin2thetac = vcl::max(1.0 - b2inv * n2inv, 0.0);
    }

    double_vt bandwidth;
    if(variable_bandwidth_spline_) {
      bandwidth = variable_bandwidth_spline_->vcl_value<VCLArchitecture>(track_x_.z());
    } else {
      bandwidth = fixed_bandwidth_;
    }

    double_vt yield_per_ev = track_yield_const_ * sin2thetac;
    double_vt mfp = cherenkov_weight_/(bandwidth * yield_per_ev);
    double_vt dx_emission = vcl::min(mfp * rng_->exponential_double(), track_dx_);

    double_vt sum_yield_z_to_the_n = dx_emission*yield_per_ev;
    sum_yield_ += sum_yield_z_to_the_n;
    sum_yield_z_to_the_n *= track_x_.z();
    sum_yield_h_ += sum_yield_z_to_the_n;
    sum_yield_z_to_the_n *= track_x_.z();
    sum_yield_h2_ += sum_yield_z_to_the_n;

    track_dx_ -= dx_emission;
    track_valid_ &= track_dx_>0;

    track_x_ += track_u_ * dx_emission;
    track_t_ += track_dt_dx_ * dx_emission;
    track_g_ += track_dg_dx_ * dx_emission;

    unsigned nvalid = horizontal_count(track_valid_);
    if(nvalid) {
      num_rays_ += nvalid;

      double_vt cos_thetac = vcl::sqrt(1.0 - sin2thetac);
      double_vt sin_thetac = vcl::sqrt(sin2thetac);

      double_vt cos_phi;
      double_vt sin_phi;
      rng_->sincos_double(sin_phi, cos_phi);

      Vector3d_vt v(cos_phi*sin_thetac, sin_phi*sin_thetac,cos_thetac);
      calin::math::geometry::VCL<double_real>::rotate_in_place_z_to_u_Rzy(v, track_u_);

      calin::math::ray::VCLRay<double_real> rays(track_x_, v, track_t_);

      propagate_rays(rays, track_valid_, bandwidth, cherenkov_weight_*track_weight_);
    }
  } while(vcl::horizontal_and(track_valid_) or
      (drain_tracks and vcl::horizontal_or(track_valid_)));
}

template<typename VCLArchitecture> void VCLIACTTrackVisitor<VCLArchitecture>::
propagate_rays(calin::math::ray::VCLRay<double_real> ray, double_bvt ray_mask,
  double_vt bandwidth, double_vt weight)
{
  // default does nothing
}

#endif // not defined SWIG

} } } // namespace calin::simulation::vcl_iact
