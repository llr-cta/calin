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
#endif // not defined SWIG

  VCLIACTTrackVisitor(calin::simulation::atmosphere::LayeredRefractiveAtmosphere* atm,
    calin::math::rng::VCLRNG<VCLArchitecture>* rng = nullptr,
    bool adopt_atm = false, bool adopt_rng = false);
  virtual ~VCLIACTTrackVisitor();
  void visit_event(const calin::simulation::tracker::Event& event, bool& kill_event) override;
  void visit_track(const calin::simulation::tracker::Track& track, bool& kill_track) override;
  void leave_event() override;

#ifndef SWIG
public:
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

  double bandwidth_ = 3.0;
  double forced_sin2theta_ = -1.0;
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
  track_dx_ = 0.0;
}

template<typename VCLArchitecture> void VCLIACTTrackVisitor<VCLArchitecture>::
visit_track(const calin::simulation::tracker::Track& track, bool& kill_track)
{
  if(track.q == 0 or track.dx <= 0)return; // it's store policy: no charge, no track = no radiation

  using calin::util::vcl::insert_into_vec3_with_mask;
  using calin::util::vcl::insert_into_with_mask;
  using calin::math::special::SQR;
  using calin::simulation::air_cherenkov_tracker::YIELD_CONST;

  double_bvt unused_tracks = track_dx_ <= 0;
  int unused_track_count = vcl::horizontal_count(unused_tracks);

  int insert_index = vcl::horizontal_find_first(unused_tracks);

  if(insert_index == -1)
    throw std::logic_error(
      calin::util::vcl::templated_class_name<VCLArchitecture>("VCLIACTTrackVisitor")
      + "::visit_track: No free SIMD vector slot");

  double_bvt insert_mask = false;
  insert_mask.insert(insert_index, true);

  vcl::Vec2d g = vcl::max(vcl::Vec2d(track.e0, track.e1)/track.mass, 1.0);

  double dx_inv = 1.0/track.dx;

  // Watch out for gamma^2 slightly less than 1.0
  // const double g2 = SQR(std::max(track.e0/track.mass,1.0)); // gamma^2
  // const double b2 = 1.0 - 1.0/g2; // beta^2

  track_valid_ |= insert_mask;

  insert_into_vec3_with_mask<double_real>(track_x_, track.x0, insert_mask);
  insert_into_with_mask<double_real>(track_t_, track.t0, insert_mask);
  insert_into_with_mask<double_real>(track_g_, g[0], insert_mask);
  insert_into_with_mask<double_real>(track_yield_const_, bandwidth_*YIELD_CONST*SQR(track.q), insert_mask);

  insert_into_vec3_with_mask<double_real>(track_u_, track.dx_hat, insert_mask);
  insert_into_with_mask<double_real>(track_dx_, track.dx, insert_mask);
  insert_into_with_mask<double_real>(track_dt_dx_, track.dt*dx_inv, insert_mask);
  insert_into_with_mask<double_real>(track_dg_dx_, (g[1]-g[0])*dx_inv, insert_mask);

  if(unused_track_count == 1)
  {
    // All vector slots are full now
    generate_mc_rays();
  }

/*

  if(enable_forced_cherenkov_angle_mode_)
    cherenkov.sin2_thetac  = forced_sin2_thetac_;
  else
    cherenkov.sin2_thetac  = 1.0 - 1.0/(b2*SQR(cherenkov.n));
  if(cherenkov.sin2_thetac <= 0.0)return;
  cherenkov.yield_density  =
    YIELD_CONST*SQR(track.q)*cherenkov.sin2_thetac*cherenkov.dx;
  cherenkov.cos_thetac     = std::sqrt(1.0 - cherenkov.sin2_thetac);
  cherenkov.sin_thetac     = std::sqrt(cherenkov.sin2_thetac);

  cherenkov.dx             = track.dx;
  cherenkov.de             = track.de;
  cherenkov.dt             = track.dt;


  Vector3d_vt track_x0; // track : position at start
  Vector3d_vt track_u;  // track : normalized direction vector
  double_vt track_dx_;  // track : distance from x0 to end
  double_vt track_b2_;  // track : beta^2 of particle
*/
}

template<typename VCLArchitecture> void VCLIACTTrackVisitor<VCLArchitecture>::
leave_event()
{
  // nothing to see here
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
    double_vt yield = track_yield_const_ * sin2thetac;
    double_vt mfp = 1.0/yield;
    double_vt dx_emission = vcl::min(mfp * rng_->exponential_double(), track_dx_);

    track_dx_ -= dx_emission;
    track_valid_ &= track_dx_>0;

    track_x_ += track_u_ * dx_emission;
    track_t_ += track_dt_dx_ * dx_emission;
    track_g_ += track_dg_dx_ * dx_emission;

    nmo *= (1.0 + dlognmo_dx * dx_emission);

    if(horizontal_or(track_valid_)) {
      if(forced_sin2theta_ > 0) {
        // nothing to see here
      } else {
        double_vt n2inv = 1.0/SQR(nmo + 1.0);
        double_vt g2 = SQR(track_g_);
        double_vt b2inv = g2/(g2-1.0);
        sin2thetac = vcl::max(1.0 - b2inv * n2inv, 0.0);
      }
      propagate_rays(sin2thetac);
    }

    if(--max_loop == 0)
    {
      LOG(WARNING) << "IN   z: " << track_x_.z() << '\n'
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

}

#endif // not defined SWIG

} } } // namespace calin::simulation::vcl_iact
