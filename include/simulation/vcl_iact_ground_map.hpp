/*

   calin/simulation/vcl_iact_ground_map.hpp -- Stephen Fegan -- 2022-06-25

   Class for imaging atmospheric cherenkov technique - produce rays from charged
   tracks in the atmosphere, propagate them to the ground and make a map of
   where they arrive and from which direction.

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

#include <vector>

#include <math/special.hpp>
#include <simulation/vcl_iact.hpp>
#include <simulation/vcl_iact.pb.h>

namespace calin { namespace simulation { namespace vcl_iact {

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCLIACTGroundMap:
  public VCLIACTTrackVisitor<VCLArchitecture>
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

  VCLIACTGroundMap(calin::simulation::atmosphere::LayeredRefractiveAtmosphere* atm,
    const calin::ix::simulation::vcl_iact::VCLIACTGroundMapConfiguration& config = default_config(),
    calin::math::rng::VCLRNG<VCLArchitecture>* rng = nullptr,
    bool adopt_atm = false, bool adopt_rng = false);
  virtual ~VCLIACTGroundMap();

  static calin::ix::simulation::vcl_iact::VCLIACTGroundMapConfiguration default_config();

  const calin::ix::simulation::vcl_iact::VCLIACTGroundMapConfiguration& config() const { return config_; }
  double dzatm_profile() const { return 1.0/dzatm_profile_inv_; }

  const std::vector<double>& xatm(unsigned idetector) const { 
    return detector_.at(idetector)->xatm; }
  const std::vector<double>& yatm(unsigned idetector) const { 
    return detector_.at(idetector)->yatm; }
  const std::vector<double>& zatm(unsigned idetector) const { 
    return detector_.at(idetector)->zatm; }
  const std::vector<double>& tgnd(unsigned idetector) const { 
    return detector_.at(idetector)->tgnd; }
  const std::vector<double>& xgnd(unsigned idetector) const { 
    return detector_.at(idetector)->xgnd; }
  const std::vector<double>& ygnd(unsigned idetector) const { 
    return detector_.at(idetector)->ygnd; }
  const std::vector<double>& uxgnd(unsigned idetector) const { 
    return detector_.at(idetector)->uxgnd; }
  const std::vector<double>& uygnd(unsigned idetector) const { 
    return detector_.at(idetector)->uygnd; }

  double num_cherenkov() const { return ncherenkov_; }
  const std::vector<double>& zatm_profile() const { return zatm_profile_; }

#ifndef SWIG
  void visit_event(const calin::simulation::tracker::Event& event, bool& kill_event) final;
  void propagate_rays(calin::math::ray::VCLRay<double_real> ray, double_bvt ray_mask,
    double_vt bandwidth, double_vt ray_weight) final;

protected:
  struct Detector {
    double r2 = 0;
    const calin::ix::simulation::vcl_iact::VCLIACTGroundMapDetectorConfiguration* config = nullptr;
    std::vector<double> xatm;
    std::vector<double> yatm;
    std::vector<double> zatm;
    std::vector<double> tgnd;
    std::vector<double> xgnd;
    std::vector<double> ygnd;
    std::vector<double> uxgnd;
    std::vector<double> uygnd;
    void clear() {
      xatm.clear();
      yatm.clear();
      zatm.clear();
      tgnd.clear();
      xgnd.clear();
      ygnd.clear();
      uxgnd.clear();
      uygnd.clear();
    }
  };

  calin::ix::simulation::vcl_iact::VCLIACTGroundMapConfiguration config_;

  double ncherenkov_;
  std::vector<double> zatm_profile_;
  double dzatm_profile_inv_;

  std::vector<Detector*> detector_;

  unsigned iobs_ = 0;
  double zobs_ = 0.0;
#endif
};

#ifndef SWIG

template<typename VCLArchitecture> VCLIACTGroundMap<VCLArchitecture>::
VCLIACTGroundMap(
    calin::simulation::atmosphere::LayeredRefractiveAtmosphere* atm,
    const calin::ix::simulation::vcl_iact::VCLIACTGroundMapConfiguration& config,
    calin::math::rng::VCLRNG<VCLArchitecture>* rng, 
    bool adopt_atm, bool adopt_rng):
  VCLIACTTrackVisitor<VCLArchitecture>(atm, config.base_config(), rng, adopt_atm, adopt_rng),
  config_(config), dzatm_profile_inv_(1.0/config.dzatm_profile()),
  iobs_(config.observation_level()), zobs_(atm->zobs(iobs_))
{
  using calin::math::special::SQR;

  zatm_profile_.resize(std::ceil(atm->top_of_atmosphere() * dzatm_profile_inv_));
  detector_.resize(config_.detector_size());
  for(int idetector=0;idetector<config_.detector_size(); idetector++) {
    detector_[idetector] = new Detector;
    detector_[idetector]->r2 = SQR(config_.detector(idetector).r_gnd());
    detector_[idetector]->config = &config_.detector(idetector);
  }
}

template<typename VCLArchitecture> VCLIACTGroundMap<VCLArchitecture>::
~VCLIACTGroundMap()
{
  for(auto* idetector : detector_) {
    delete idetector;
  }
}

template<typename VCLArchitecture> 
calin::ix::simulation::vcl_iact::VCLIACTGroundMapConfiguration 
VCLIACTGroundMap<VCLArchitecture>::default_config()
{
  calin::ix::simulation::vcl_iact::VCLIACTGroundMapConfiguration config;
  config.mutable_base_config()->CopyFrom(VCLIACTTrackVisitor<VCLArchitecture>::default_config());
  auto* detector = config.add_detector();
  detector->set_x_gnd(0);
  detector->set_y_gnd(0);
  detector->set_r_gnd(1e5 /* 1km */);
  detector->set_store_position(true);
  detector->set_store_direction(true);
  detector->set_store_time(true);
  detector->set_store_emission_point(false);
  detector->set_store_fraction(1.0);
  config.set_dzatm_profile(5000.0);
  config.set_observation_level(0);
  return config;
}

template<typename VCLArchitecture> void VCLIACTGroundMap<VCLArchitecture>::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
  ncherenkov_ = 0.0;
  std::fill(zatm_profile_.begin(), zatm_profile_.end(), 0.0);
  for(auto& idetector : detector_) {
    idetector->clear();
  }
  return VCLIACTTrackVisitor<VCLArchitecture>::visit_event(event, kill_event);
}

template<typename VCLArchitecture> void VCLIACTGroundMap<VCLArchitecture>::
propagate_rays(calin::math::ray::VCLRay<double_real> ray, double_bvt ray_mask,
  double_vt bandwidth, double_vt ray_weight)
{
  double_at atm_x;
  double_at atm_y;
  double_at atm_z;
  double_at ray_w;
  ray.x().store(atm_x);
  ray.y().store(atm_y);
  ray.z().store(atm_z);
  ray_weight.store(ray_w);

  for(unsigned iv=0; iv<VCLArchitecture::num_double; iv++) {
    if(ray_mask[iv]) {
      ncherenkov_ += ray_w[iv];
      zatm_profile_.at(atm_z[iv] * dzatm_profile_inv_) += ray_w[iv];
    }
  }

  ray_mask &= (ray.z()>zobs_) & (ray.uz()<0);

  // Note that uz is negative so we subtract delay
  ray.mutable_ct() -= vcl::select(ray_mask,
    this->atm_->template vcl_propagation_ct_correction_to_iobs<VCLArchitecture>(ray.z(), 0)/ray.uz(),
    0);

  ray_mask = ray.propagate_to_z_plane_with_mask(ray_mask, zobs_, false);

  double_at t;
  double_at x;
  double_at y;
  double_at ux;
  double_at uy;
  ray.time().store(t);
  ray.x().store(x);
  ray.y().store(y);
  ray.ux().store(ux);
  ray.uy().store(uy);

  double_vt u = VCLIACTTrackVisitor<VCLArchitecture>::rng_->uniform_double();

  for(unsigned idetector=0;idetector<config_.detector_size();idetector++) {
    auto* detector = detector_[idetector];
    double_vt xrel = ray.x() - detector->config->x_gnd();
    double_vt yrel = ray.y() - detector->config->y_gnd();
    double_bvt store_mask = ray_mask && (xrel*xrel + yrel*yrel)<detector->r2 && u<=detector->config->store_fraction();

      for(unsigned iv=0; iv<VCLArchitecture::num_double; iv++) {
        if(store_mask[iv]) {
          if(detector->config->store_position()) {
            detector->xgnd.push_back(x[iv]);
            detector->ygnd.push_back(y[iv]);
          }
          if(detector->config->store_direction()) {
            detector->uxgnd.push_back(ux[iv]);
            detector->uygnd.push_back(uy[iv]);
          }
          if(detector->config->store_time()) {
            detector->tgnd.push_back(t[iv]);
          }
          if(detector->config->store_time()) {
            detector->xatm.push_back(atm_x[iv]);
            detector->yatm.push_back(atm_y[iv]);
            detector->zatm.push_back(atm_z[iv]);
          }
        }
      }
  }
}

#endif // not defined SWIG

} } } // namespace calin::simulation::vcl_iact
