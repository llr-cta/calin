/*

   calin/simulation/vcl_iact_array.hpp -- Stephen Fegan -- 2022-06-29

   Class for imaging atmospheric cherenkov technique - produce rays from charged
   tracks in the atmosphere, propagate them to the ground and trace them through
   telescope ptics.

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
#include <simulation/vcl_raytracer.hpp>

namespace calin { namespace simulation { namespace vcl_iact {

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCLIACTArray:
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

  using Ray         = calin::math::ray::VCLRay<double_real>;

  using VCLFocalPlaneParameters = calin::simulation::vcl_ray_propagator::VCLFocalPlaneParameters<VCLArchitecture>;

  using calin::simulation::ray_processor::RayProcessorDetectorSphere;
#endif // not defined SWIG

  CALIN_TYPEALIAS(VCLFocalPlaneRayPropagator, calin::simulation::vcl_ray_propagator::VCLFocalPlaneRayPropagator<VCLArchitecture>);

  VCLIACTArray(calin::simulation::atmosphere::LayeredRefractiveAtmosphere* atm,
    const calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration& config = default_config(),
    calin::math::rng::VCLRNG<VCLArchitecture>* rng = nullptr,
    bool adopt_atm = false, bool adopt_rng = false);
  virtual ~VCLIACTArray();

  static calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration default_config();

#ifndef SWIG
  void visit_event(const calin::simulation::tracker::Event& event, bool& kill_event) final;
  void propagate_rays(calin::math::ray::VCLRay<double_real> ray, double_bvt ray_mask, double_vt ray_weight) final;
  void leave_event() final;

protected:
  static calin::ix::simulation::vcl_iact::VCLIACTConfiguration base_config(
    const calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration& config);

  struct PropagatorInfo {
    VCLFocalPlaneRayPropagator* propagator;
    unsigned detector0;
    unsigned ndetector;
  };

  struct DetectorInfo {
    RayProcessorDetectorSphere sphere;
    double squared_radius;
    double squared_safety_radius;
    VCLFocalPlaneRayPropagator* propagator;
    unsigned propagator_iscope;
    unsigned global_iscope;
    double ref_index;

    VCLRayArray rays_to_refract;
    unsigned nrays_to_refract;
    double_at ray_weights_to_refract;

    VCLRayArray rays_to_propagate;
    unsigned nrays_to_propagate;
    double_at ray_weights_to_propagate;
  };

  void do_propagate_rays_for_detector(DetectorInfo& idetector);
  void do_refract_rays_for_detector(DetectorInfo& idetector);

  std::vector<PropagatorInfo> propoagator_;
  std::vector<DetectorInfo> detector_;
  calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration config_;
  double zobs_;
  double ref_index_;
  double ref_index_inv_;
  double ref_index_correction_;
#endif
};

#ifndef SWIG

template<typename VCLArchitecture> VCLIACTArray<VCLArchitecture>::
VCLIACTArray(
    calin::simulation::atmosphere::LayeredRefractiveAtmosphere* atm,
    const calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration& config,
    calin::math::rng::VCLRNG<VCLArchitecture>* rng,
    bool adopt_atm, bool adopt_rng):
  VCLIACTTrackVisitor<VCLArchitecture>(atm, base_config(config), rng, adopt_atm, adopt_rng),
  config_(config), ref_index_correction_(config_.ref)
{
  if(config_.observation_level >= atm.num_obs_levels) {
    throw std::out_of_range("Request observation level out of range.");
  }
  zobs_ = atm.zobs(config_.observation_level);
  double n_minus_one = atm.n_minus_one(zobs_);
  ref_index_ = 1.0 + n_minus_one;
  ref_index_correction_ = -n_minus_one;
}

template<typename VCLArchitecture> VCLIACTGroundMap<VCLArchitecture>::
~VCLIACTArray()
{
  // nothing to see here
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
  for(const auto& ipropoagator in propoagator) {
    auto spheres = ipropoagator.propagator.detector_spheres();
    if(spheres.size() != ipropoagator.ndetector) {
      // this should never happen
      throw std::runtime_error("Number of detectors proposed by propagator must remain constant over events.");
    }
    for(unsigned isphere=0; isphere<spheres.size(); ++isphere) {
      if(sphere.iobs != config_.observation_level()) {
        throw std::runtime_error("Detector observation level does not match configured value.");
      }
      detector[ipropoagator.detector0+isphere].sphere = spheres[isphere];
      detector[ipropoagator.detector0+isphere].nrays_to_refract = 0;
      detector[ipropoagator.detector0+isphere].nrays_to_propagate = 0;
    }
  }

  return VCLIACTTrackVisitor<VCLArchitecture>::visit_event(event, kill_event);
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
propagate_rays(calin::math::ray::VCLRay<double_real> ray, double_bvt ray_mask, double_vt ray_weight)
{
  ray_mask &= (ray.z()>zobs_) & (ray.uz()<0);

  switch(double_vt dz = ray.z-zobs_; config_.refraction_mode()) {
  case calin::ix::simulation::vcl_iact::REFRACT_NO_RAYS:
    // Note use of "-=" below because uz() is negative
    ray.mutable_ct() -= vcl::select(ray_mask,
      (dz*ref_index_correction_ + atm_.propagation_ct_correction_to_iobs(ray.z(), config_.observation_level()))/ray.uz(),
      0);
    break;
  case calin::ix::simulation::vcl_iact::REFRACT_ALL_RAYS:
    atm.vcl_propagate_ray_with_refraction_and_mask(ray, ray_mask, config_.observation_level());
    // Note propagation backwards by distance since uz() is negative
    ray.propagate_dist_with_mask(dz/ray.uz(), ray_mask, ref_index_);
    break;
  case calin::ix::simulation::vcl_iact::REFRACT_ONLY_CLOSE_RAYS:
    break;
  }

  calin::math::ray::VCLRayArray<double_real> ray_array { ray };
  double_at ray_weight_array
  ray_weight.store(ray_weight_array);

  for(auto& idetector : detector_) {
    auto intersecting_rays = ray_mask & (ray.squared_distance_at_closest_approach() < detector.squared_safety_radius);
    unsigned intersecting_rays_bitmask = vcl::to_bits(intersecting_rays);
    if(intersecting_rays_bitmask) {
      for(unsigned iray=0; iray<VCLArchitecture::num_double; ++iray) {
        if(intersecting_rays_bitmask & 1) {
          switch(double_vt dz = ray.z-zobs_; config_.refraction_mode()) {
          case calin::ix::simulation::vcl_iact::REFRACT_NO_RAYS:
          case calin::ix::simulation::vcl_iact::REFRACT_ALL_RAYS:
            idetector.rays_to_propagate.insert(detector.nrays_to_propagate, ray_array.extract(iray));
            idetector.ray_weights_to_propagate[detector.nrays_to_propagate] = ray_weight_array[iray];
            ++idetector.nrays_to_propagate;
            if(idetector.nrays_to_propagate == VCLArchitecture::num_double) {
              do_propagate_rays_for_detector(idetector)
            }
            break;
          case calin::ix::simulation::vcl_iact::REFRACT_ONLY_CLOSE_RAYS:
            idetector.rays_to_refract.insert(detector.nrays_to_refract, ray_array.extract(iray));
            idetector.ray_weights_to_refract[detector.nrays_to_refract] = ray_weight_array[iray];
            ++idetector.nrays_to_refract;
            if(idetector.nrays_to_refract == VCLArchitecture::num_double) {
              do_refract_rays_for_detector(idetector)
            }
            break;

        }
        intersecting_rays_bitmask >>= 1;
      }
    }
  }

  ray_mask = ray.propagate_to_z_plane_with_mask(ray_mask,
    VCLIACTTrackVisitor<VCLArchitecture>::atm_->zobs(0), false);
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
do_refract_rays_for_detector(DetectorInfo& idetector)
{
  Ray ray;
  idetector.rays_to_refract.get_rays(ray);
  double_bvt ray_mask = VCLArchitecture::double_iota()<idetector.nrays_to_refract;
  idetector.nrays_to_refract = 0;

  double_vt dz = ray.z-zobs_;
  atm.vcl_propagate_ray_with_refraction_and_mask(ray, ray_mask, config_.observation_level());
  // Note propagation backwards by distance since uz() is negative
  ray.propagate_dist_with_mask(dz/ray.uz(), ray_mask, ref_index_);

  idetector.rays_to_refract.set_rays(ray);

  auto intersecting_rays = ray_mask & (ray.squared_distance_at_closest_approach() < detector.squared_radius);
  unsigned intersecting_rays_bitmask = vcl::to_bits(intersecting_rays);
  if(intersecting_rays_bitmask) {
    for(unsigned iray=0; iray<VCLArchitecture::num_double; ++iray) {
      if(intersecting_rays_bitmask & 1) {
        idetector.rays_to_propagate.insert(detector.nrays_to_propagate, idetector.rays_to_refract.extract(iray));
        idetector.ray_weights_to_propagate[detector.nrays_to_propagate] = idetector.ray_weights_to_refract[iray];
        ++idetector.nrays_to_propagate;
        if(idetector.nrays_to_propagate == VCLArchitecture::num_double) {
          do_propagate_rays_for_detector(idetector)
        }
      }
    }
  }
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
do_propagate_rays_for_detector(DetectorInfo& idetector)
{
  Ray ray;
  idetector.rays_to_propagate.get_rays(ray);
  double_bvt ray_mask = VCLArchitecture::double_iota()<idetector.nrays_to_propagate;
  idetector.nrays_to_propagate = 0;

  VCLFocalPlaneParameters<VCLArchitecture> fp_parameters;
  ray_mask = idetector.propagator->propagate_rays_to_focal_plane(
    idetector.propagator_iscope, ray, ray_mask, fp_parameters);

  
}

template<typename VCLArchitecture> calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration
VCLIACTArray<VCLArchitecture>::default_config()
{
  calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration config;
  config.set_bandwidth(3.0);
  return config;
}

template<typename VCLArchitecture> calin::ix::simulation::vcl_iact::VCLIACTConfiguration
VCLIACTArray<VCLArchitecture>::base_config(const calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration& config)
{
  calin::ix::simulation::vcl_iact::VCLIACTConfiguration bconfig;
  bconfig.set_bandwidth(config.bandwidth());
  bconfig.set_enable_forced_cherenkov_angle_mode(config.enable_forced_cherenkov_angle_mode());
  bconfig.set_forced_cherenkov_angle(config.forced_cherenkov_angle());
  return bconfig;
}

#endif // not defined SWIG

} } } // namespace calin::simulation::vcl_iact
