/*

   calin/simulation/vcl_iact.hpp -- Stephen Fegan -- 2019-02-19

   Class for imaging atmospheric cherenkov technique - produce rays from charged
   tracks in the atmosphere, propagate them to the ground and trace them through
   optics.

   Copyright 2019, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

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

namespace calin { namespace simulation { namespace vcl_iact {

template<typename VCLArchitecture> class VCLIACTTrackVisitor:
  public VCLArchitecture, public calin::simulation::tracker::TrackVisitor
{
public:
#ifndef SWIG
  using typename VCLArchitecture::uint32_vt;
  using typename VCLArchitecture::int32_vt;
  using typename VCLArchitecture::uint64_vt;
  using typename VCLArchitecture::int64_vt;
  using typename VCLArchitecture::float_vt;
  using typename VCLArchitecture::double_vt;
  using typename VCLArchitecture::uint64_bvt;
  using typename VCLArchitecture::float_bvt;
  using typename VCLArchitecture::double_bvt;
  using typename VCLArchitecture::Vector3f_vt;
  using typename VCLArchitecture::Vector3d_vt;
  using typename VCLArchitecture::float_real;
  using typename VCLArchitecture::double_real;
  using typename VCLArchitecture::float_at;
  using typename VCLArchitecture::double_at;
#endif // not defined SWIG

  VCLIACTTrackVisitor(calin::simulation::atmosphere::LayeredRefractiveAtmosphere* atm,
    calin::math::rng::VCLRNG<VCLArchitecture>* rng = nullptr,
    bool adopt_atm = false, bool adopt_rng = false);
  virtual ~VCLIACTTrackVisitor();
  void visit_event(const calin::simulation::tracker::Event& event, bool& kill_event) override;
  void visit_track(const calin::simulation::tracker::Track& track, bool& kill_track) override;
  void leave_event() override;

  uint64_t num_tracks() const { return num_tracks_; }
  uint64_t num_steps() const { return num_steps_; }
  uint64_t num_rays() const { return num_rays_; }

  const std::vector<double>& xgnd() const { return xgnd_; }
  const std::vector<double>& ygnd() const { return ygnd_; }
  const std::vector<double>& uxgnd() const { return uxgnd_; }
  const std::vector<double>& uygnd() const { return uygnd_; }

#ifndef SWIG
public:
  inline int insert_track(const Eigen::Vector3d& x, const double t, const double g,
    const double yield, const Eigen::Vector3d& u, const double dx, const double dt_dx,
    const double dg_dx, bool valid = true);
  void generate_mc_rays();
  void propagate_rays(double_vt sin2thetac);

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

  uint64_t num_tracks_ = 0;
  uint64_t num_steps_ = 0;
  uint64_t num_rays_ = 0;

  double bandwidth_ = 3.0;
  double forced_sin2theta_ = -1.0;

  std::vector<double> xgnd_;
  std::vector<double> ygnd_;
  std::vector<double> uxgnd_;
  std::vector<double> uygnd_;
#endif // not defined SWIG
};

#ifndef SWIG

template<typename VCLArchitecture> VCLIACTTrackVisitor<VCLArchitecture>::
VCLIACTTrackVisitor(
    calin::simulation::atmosphere::LayeredRefractiveAtmosphere* atm,
    calin::math::rng::VCLRNG<VCLArchitecture>* rng,
    bool adopt_atm, bool adopt_rng):
  VCLArchitecture(), calin::simulation::tracker::TrackVisitor(),
  atm_(atm), adopt_atm_(adopt_atm),
  rng_(rng ? rng : new calin::math::rng::VCLRNG<VCLArchitecture>()),
  adopt_rng_(rng ? adopt_rng : true)
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
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
  num_tracks_ = 0;
  num_steps_ = 0;
  num_rays_ = 0;
  track_dx_ = 0.0;
  track_valid_ = false;

  xgnd_.clear();
  ygnd_.clear();
  uxgnd_.clear();
  uygnd_.clear();
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
    insert_track(track.x0, track.t0, g[0], bandwidth_*YIELD_CONST*SQR(track.q),
      track.dx_hat, dx, track.dt*dx_inv, (g[1]-g[0])*dx_inv);

  if(unused_track_count == 1)
  {
    // All vector slots are full now
    generate_mc_rays();
  }
}

template<typename VCLArchitecture> inline int VCLIACTTrackVisitor<VCLArchitecture>::
insert_track(const Eigen::Vector3d& x, const double t, const double g,
  const double yield, const Eigen::Vector3d& u, const double dx, const double dt_dx,
  const double dg_dx, bool valid)
{
  using calin::util::vcl::insert_into_vec3_with_mask;
  using calin::util::vcl::insert_into_with_mask;
  using namespace calin::util::log;

  double_bvt unused_tracks = track_dx_ <= 0;
  int insert_index = vcl::horizontal_find_first(unused_tracks);

  if(insert_index == -1)
    throw std::logic_error(
      calin::util::vcl::templated_class_name<VCLArchitecture>("VCLIACTTrackVisitor")
      + "::insert_track: No free SIMD vector slot");

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

    generate_mc_rays();
  }
}

template<typename VCLArchitecture> void VCLIACTTrackVisitor<VCLArchitecture>::
generate_mc_rays()
{
  using calin::math::special::SQR;
  using namespace calin::util::log;

  double_vt nmo; // n minus one
  double_vt dlognmo_dx; // rate of change of log(n minus one) with track length
  double_vt sin2thetac; // sin squared of the Cherenkov angle

  if(forced_sin2theta_ > 0) {
    nmo = 0.0;
    dlognmo_dx = 0.0;
    sin2thetac = forced_sin2theta_;
  } else {
    dlognmo_dx =
      atm_->vcl_dlognmo_dz<VCLArchitecture>(track_x_.z(), nmo) * track_u_.z();
    double_vt n2inv = 1.0/SQR(nmo + 1.0);
    double_vt g2 = SQR(track_g_);
    double_vt b2inv = g2/(g2-1.0);
    sin2thetac = vcl::max(1.0 - b2inv * n2inv, 0.0);
  }

  double max_loop = 10000;
  do {
    ++num_steps_;

    double_vt yield = track_yield_const_ * sin2thetac;
    double_vt mfp = 1.0/yield;
    double_vt dx_emission = vcl::min(mfp * rng_->exponential_double(), track_dx_);

    track_dx_ -= dx_emission;
    track_valid_ &= track_dx_>0;

    track_x_ += track_u_ * dx_emission;
    track_t_ += track_dt_dx_ * dx_emission;
    track_g_ += track_dg_dx_ * dx_emission;

    nmo *= (1.0 + dlognmo_dx * dx_emission);

    unsigned nvalid = horizontal_count(track_valid_);
    if(nvalid) {
      if(forced_sin2theta_ > 0) {
        // nothing to see here
      } else {
        double_vt n2inv = 1.0/SQR(nmo + 1.0);
        double_vt g2 = SQR(track_g_);
        double_vt b2inv = g2/(g2-1.0);
        sin2thetac = vcl::max(1.0 - b2inv * n2inv, 0.0);
      }

      num_rays_ += nvalid;

      propagate_rays(sin2thetac);
    }

    if(--max_loop == 0)
    {
      LOG(WARNING) << " valid: " << track_valid_ << '\n'
                   << "     x: " << track_x_.x() << '\n'
                   << "     y: " << track_x_.y() << '\n'
                   << "     z: " << track_x_.z() << '\n'
                   << "    ux: " << track_u_.x() << '\n'
                   << "    uy: " << track_u_.y() << '\n'
                   << "    uz: " << track_u_.z() << '\n'
                   << "   n-1: " << nmo << '\n'
                   << "    dx: " << track_dx_ << '\n'
                   << "    Yc: " << track_yield_const_ << '\n'
                   << " gamma: " << track_g_ << '\n'
                   << " sin2c: " << sin2thetac << '\n'
                   << "     Y: " << yield << '\n'
                   << "   MFP: " << mfp;
      throw std::runtime_error("Maximum loop exceeded");
    }
  } while(vcl::horizontal_and(track_dx_ > 0));
}

template<typename VCLArchitecture> void VCLIACTTrackVisitor<VCLArchitecture>::
propagate_rays(double_vt sin2thetac)
{
  double_vt cos_thetac = vcl::sqrt(1.0 - sin2thetac);
  double_vt sin_thetac = vcl::sqrt(sin2thetac);

  double_vt cos_phi;
  double_vt sin_phi;
  rng_->sincos_double(sin_phi, cos_phi);

  Vector3d_vt v(cos_phi*sin_thetac, sin_phi*sin_thetac,cos_thetac);
  calin::math::geometry::VCL<double_real>::rotate_in_place_z_to_u_Rzy(v, track_u_);

  calin::math::ray::VCLRay<double_real> ray(track_x_, v, track_t_);
  double_bvt mask = ray.propagate_to_z_plane_with_mask(track_valid_, atm_->zobs(0), false);

  double_at x;
  double_at y;
  double_at ux;
  double_at uy;
  ray.x().store(x);
  ray.y().store(y);
  ray.ux().store(ux);
  ray.uy().store(uy);

  for(unsigned iv=0; iv<VCLArchitecture::num_double; iv++) {
    if(mask[iv]) {
      xgnd_.push_back(x[iv]);
      ygnd_.push_back(y[iv]);
      uxgnd_.push_back(ux[iv]);
      uygnd_.push_back(uy[iv]);
    }
  }
}

#endif // not defined SWIG

} } } // namespace calin::simulation::vcl_iact
