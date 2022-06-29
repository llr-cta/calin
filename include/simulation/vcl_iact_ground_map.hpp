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

#include <simulation/vcl_iact.hpp>

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
    const calin::ix::simulation::vcl_iact::VCLIACTConfiguration& config = VCLIACTTrackVisitor<VCLArchitecture>::default_config(),
    calin::math::rng::VCLRNG<VCLArchitecture>* rng = nullptr,
    bool adopt_atm = false, bool adopt_rng = false);
  virtual ~VCLIACTGroundMap();

  const std::vector<double>& xgnd() const { return xgnd_; }
  const std::vector<double>& ygnd() const { return ygnd_; }
  const std::vector<double>& uxgnd() const { return uxgnd_; }
  const std::vector<double>& uygnd() const { return uygnd_; }

#ifndef SWIG
  void visit_event(const calin::simulation::tracker::Event& event, bool& kill_event) final;
  void propagate_rays(calin::math::ray::VCLRay<double_real> ray, double_bvt ray_mask, double_vt ray_weight) final;

protected:
  std::vector<double> xgnd_;
  std::vector<double> ygnd_;
  std::vector<double> uxgnd_;
  std::vector<double> uygnd_;
#endif
};

#ifndef SWIG

template<typename VCLArchitecture> VCLIACTGroundMap<VCLArchitecture>::
VCLIACTGroundMap(
    calin::simulation::atmosphere::LayeredRefractiveAtmosphere* atm,
    const calin::ix::simulation::vcl_iact::VCLIACTConfiguration& config,
    calin::math::rng::VCLRNG<VCLArchitecture>* rng,
    bool adopt_atm, bool adopt_rng):
  VCLIACTTrackVisitor<VCLArchitecture>(atm, config, rng, adopt_atm, adopt_rng)
{
  // nothing to see here
}

template<typename VCLArchitecture> VCLIACTGroundMap<VCLArchitecture>::
~VCLIACTGroundMap()
{
  // nothing to see here
}

template<typename VCLArchitecture> void VCLIACTGroundMap<VCLArchitecture>::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
  xgnd_.clear();
  ygnd_.clear();
  uxgnd_.clear();
  uygnd_.clear();

  return VCLIACTTrackVisitor<VCLArchitecture>::visit_event(event, kill_event);
}

template<typename VCLArchitecture> void VCLIACTGroundMap<VCLArchitecture>::
propagate_rays(calin::math::ray::VCLRay<double_real> ray, double_bvt ray_mask, double_vt ray_weight)
{
  ray_mask = ray.propagate_to_z_plane_with_mask(ray_mask,
    VCLIACTTrackVisitor<VCLArchitecture>::atm_->zobs(0), false);

  double_at x;
  double_at y;
  double_at ux;
  double_at uy;
  ray.x().store(x);
  ray.y().store(y);
  ray.ux().store(ux);
  ray.uy().store(uy);

  for(unsigned iv=0; iv<VCLArchitecture::num_double; iv++) {
    if(ray_mask[iv]) {
      xgnd_.push_back(x[iv]);
      ygnd_.push_back(y[iv]);
      uxgnd_.push_back(ux[iv]);
      uygnd_.push_back(uy[iv]);
    }
  }
}

#endif // not defined SWIG

} } } // namespace calin::simulation::vcl_iact
