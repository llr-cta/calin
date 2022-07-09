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

#include <util/log.hpp>
#include <math/special.hpp>
#include <simulation/vcl_iact.hpp>
#include <simulation/vcl_ray_propagator.hpp>
#include <simulation/vcl_raytracer.hpp>
#include <simulation/detector_efficiency.hpp>
#include <simulation/atmosphere.hpp>

namespace calin { namespace simulation { namespace vcl_iact {

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCLIACTArray:
  public VCLIACTTrackVisitor<VCLArchitecture>
{
public:
#ifndef SWIG
  using uint32_vt   = typename VCLArchitecture::uint32_vt;
  using int32_vt    = typename VCLArchitecture::int32_vt;
  using uint64_vt   = typename VCLArchitecture::uint64_vt;
  using int64_vt    = typename VCLArchitecture::int64_vt;
  using float_vt    = typename VCLArchitecture::float_vt;
  using double_vt   = typename VCLArchitecture::double_vt;
  using uint64_bvt  = typename VCLArchitecture::uint64_bvt;
  using float_bvt   = typename VCLArchitecture::float_bvt;
  using double_bvt  = typename VCLArchitecture::double_bvt;
  using Vector3f_vt = typename VCLArchitecture::Vector3f_vt;
  using Vector3d_vt = typename VCLArchitecture::Vector3d_vt;
  using float_real  = typename VCLArchitecture::float_real;
  using double_real = typename VCLArchitecture::double_real;
  using int64_at    = typename VCLArchitecture::int64_at;
  using uint64_at   = typename VCLArchitecture::uint64_vt;
  using float_at    = typename VCLArchitecture::float_at;
  using double_at   = typename VCLArchitecture::double_at;

  using Ray         = calin::math::ray::VCLRay<double_real>;
  using RayArray    = calin::math::ray::VCLRayArray<double_real>;
  using FocalPlaneParameters = calin::simulation::vcl_ray_propagator::VCLFocalPlaneParameters<VCLArchitecture>;

  using RayProcessorDetectorSphere = calin::simulation::ray_processor::RayProcessorDetectorSphere;
#endif // not defined SWIG

  CALIN_TYPEALIAS(PEProcessor, calin::simulation::pe_processor::PEProcessor);
  CALIN_TYPEALIAS(DetectionEfficiency, calin::simulation::detector_efficiency::DetectionEfficiency);

  CALIN_TYPEALIAS(FocalPlaneRayPropagator, calin::simulation::vcl_ray_propagator::VCLFocalPlaneRayPropagator<VCLArchitecture>);
  CALIN_TYPEALIAS(DaviesCottonVCLFocalPlaneRayPropagator, calin::simulation::vcl_ray_propagator::DaviesCottonVCLFocalPlaneRayPropagator<VCLArchitecture>);

  VCLIACTArray(calin::simulation::atmosphere::LayeredRefractiveAtmosphere* atm,
    const calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration& config = default_config(),
    calin::math::rng::VCLRNG<VCLArchitecture>* rng = nullptr,
    bool adopt_atm = false, bool adopt_rng = false);
  virtual ~VCLIACTArray();

  void add_propagator(FocalPlaneRayPropagator* propagator, PEProcessor* pe_processor,
    const DetectionEfficiency& detector_efficiency, bool adopt_propagator = false,
    bool adopt_pe_processor = false);

  DaviesCottonVCLFocalPlaneRayPropagator* add_davies_cotton_propagator(calin::simulation::vs_optics::VSOArray* array,
    PEProcessor* pe_processor, const DetectionEfficiency& detector_efficiency,
    bool adopt_array = false, bool adopt_pe_processor = false);

  DaviesCottonVCLFocalPlaneRayPropagator* add_davies_cotton_propagator(const ix::simulation::vs_optics::IsotropicDCArrayParameters& param,
    PEProcessor* pe_processor, const DetectionEfficiency& detector_efficiency,
    bool adopt_pe_processor = false);

  static calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration default_config();

  const calin::math::spline_interpolation::CubicMultiSpline& detector_efficiency_spline() const {
    return detector_efficiency_spline_;
  }

#ifndef SWIG
  void visit_event(const calin::simulation::tracker::Event& event, bool& kill_event) final;
  void propagate_rays(calin::math::ray::VCLRay<double_real> ray, double_bvt ray_mask, double_vt ray_weight) final;
  // void leave_event() final;

protected:
  static calin::ix::simulation::vcl_iact::VCLIACTConfiguration base_config(
    const calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration& config);
  std::vector<double> detector_efficiency_energy_knots();
  void add_detector_efficiency(const DetectionEfficiency& detector_efficiency, const std::string& name);

  struct PropagatorInfo {
    FocalPlaneRayPropagator* propagator;
    unsigned ipropagator;
    PEProcessor* pe_processor;
    unsigned detector0;
    unsigned ndetector;
    bool adopt_propagator;
    bool adopt_pe_processor;
  };

  struct DetectorInfo {
    RayProcessorDetectorSphere sphere;
    double squared_radius;
    double squared_safety_radius;
    FocalPlaneRayPropagator* propagator;
    unsigned propagator_iscope;
    unsigned global_iscope;
    calin::simulation::pe_processor::PEProcessor* pe_processor;
    unsigned idetector_efficiency;

    RayArray rays_to_refract;
    unsigned nrays_to_refract;
    double_at ray_weights_to_refract;

    RayArray rays_to_propagate;
    unsigned nrays_to_propagate;
    double_at ray_weights_to_propagate;
  };

  void do_propagate_rays_for_detector(DetectorInfo& idetector);
  void do_refract_rays_for_detector(DetectorInfo& idetector);

  calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration config_;
  std::vector<PropagatorInfo> propagator_;
  std::vector<DetectorInfo> detector_;
  std::vector<DetectionEfficiency> detector_efficiency_;
  calin::math::spline_interpolation::CubicMultiSpline detector_efficiency_spline_;
  double zobs_;
  double ref_index_;
  double ref_index_correction_;
  double safety_radius_;
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
  config_(config), detector_efficiency_spline_(detector_efficiency_energy_knots())
{
  if(config_.observation_level() >= this->atm_->num_obs_levels()) {
    throw std::out_of_range("Request observation level out of range.");
  }
  zobs_ = this->atm_->zobs(config_.observation_level());
  double n_minus_one = this->atm_->n_minus_one(zobs_);
  ref_index_ = 1.0 + n_minus_one;
  ref_index_correction_ = -n_minus_one;
}

template<typename VCLArchitecture> VCLIACTArray<VCLArchitecture>::
~VCLIACTArray()
{
  for(auto& ipropagator : propagator_) {
    if(ipropagator.adopt_propagator)delete ipropagator.propagator;
    if(ipropagator.adopt_pe_processor)delete ipropagator.pe_processor;
  }
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
add_propagator(FocalPlaneRayPropagator* propagator, PEProcessor* pe_processor,
  const DetectionEfficiency& detector_efficiency, bool adopt_propagator,
  bool adopt_pe_processor)
{
  using calin::math::special::SQR;

  auto sphere = propagator->detector_spheres();

  PropagatorInfo propagator_info;
  propagator_info.propagator         = propagator;
  propagator_info.ipropagator        = propagator_.size();
  propagator_info.pe_processor       = pe_processor;
  propagator_info.detector0          = detector_.size();
  propagator_info.ndetector          = sphere.size();
  propagator_info.adopt_propagator   = adopt_propagator;
  propagator_info.adopt_pe_processor = adopt_pe_processor;
  propagator_.emplace_back(propagator_info);

  for(unsigned isphere=0; isphere<sphere.size(); ++isphere) {
    DetectorInfo detector_info;
    detector_info.sphere                 = sphere[isphere];
    detector_info.squared_radius         = SQR(sphere[isphere].radius);
    detector_info.squared_safety_radius  = SQR(sphere[isphere].radius + safety_radius_);
    detector_info.propagator             = propagator;
    detector_info.propagator_iscope      = isphere;
    detector_info.global_iscope          = detector_.size();
    detector_info.pe_processor           = pe_processor;
    detector_info.idetector_efficiency   = detector_efficiency_.size();
    detector_info.nrays_to_refract       = 0;
    detector_info.nrays_to_propagate     = 0;
    detector_.emplace_back(detector_info);
  }

  add_detector_efficiency(detector_efficiency, "propagator "+std::to_string(propagator_info.ipropagator));
}

template<typename VCLArchitecture>
calin::simulation::vcl_ray_propagator::DaviesCottonVCLFocalPlaneRayPropagator<VCLArchitecture>*
VCLIACTArray<VCLArchitecture>::add_davies_cotton_propagator(
  calin::simulation::vs_optics::VSOArray* array,
  PEProcessor* pe_processor, const DetectionEfficiency& detector_efficiency,
  bool adopt_array, bool adopt_pe_processor)
{
  auto* propagator = new calin::simulation::vcl_ray_propagator::DaviesCottonVCLFocalPlaneRayPropagator<VCLArchitecture>(
    array, this->rng_, ref_index_, adopt_array, /* adopt_rng= */ false);
  add_propagator(propagator, pe_processor, detector_efficiency,
    /* adopt_propagator= */ true, adopt_pe_processor);
  return propagator;
}

template<typename VCLArchitecture>
calin::simulation::vcl_ray_propagator::DaviesCottonVCLFocalPlaneRayPropagator<VCLArchitecture>*
VCLIACTArray<VCLArchitecture>::add_davies_cotton_propagator(
  const ix::simulation::vs_optics::IsotropicDCArrayParameters& param,
  PEProcessor* pe_processor, const DetectionEfficiency& detector_efficiency,
  bool adopt_pe_processor)
{
  auto* array = new calin::simulation::vs_optics::VSOArray;
  calin::math::rng::VCLToScalarRNGCore scalar_core(this->rng_->core());
  calin::math::rng::RNG scalar_rng(&scalar_core);
  array->generateFromArrayParameters(param, scalar_rng);
  return add_davies_cotton_propagator(array, pe_processor, detector_efficiency,
    /* adopt_array= */ true, adopt_pe_processor);
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
add_detector_efficiency(const DetectionEfficiency& detector_efficiency, const std::string& name)
{
  detector_efficiency_.push_back(detector_efficiency);
  std::vector<double> detector_efficiency_xknot =
    detector_efficiency_spline_.xknot_as_stdvec();

  double full_integral = detector_efficiency.integrate();
  double partial_integral =
    detector_efficiency.integrate(config_.detector_energy_lo(), config_.detector_energy_hi());

  if(partial_integral < 0.99*full_integral) {
    calin::util::log::LOG(calin::util::log::WARNING)
      << "Configured detector energy range does contains less than 99% of detector efficiency curve.";
  }

  std::vector<double> detector_efficiency_yknot(detector_efficiency_xknot.size());
  std::transform(detector_efficiency_xknot.begin(), detector_efficiency_xknot.end(),
    detector_efficiency_yknot.begin(),
    [detector_efficiency](double e) { return detector_efficiency(e); });

  detector_efficiency_spline_.add_spline(detector_efficiency_yknot, name);
}


template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
  for(const auto& ipropoagator : propagator_) {
    auto spheres = ipropoagator.propagator->detector_spheres();
    if(spheres.size() != ipropoagator.ndetector) {
      // this should never happen
      throw std::runtime_error("Number of detectors proposed by propagator must remain constant over events.");
    }
    for(unsigned isphere=0; isphere<spheres.size(); ++isphere) {
      if(spheres[isphere].iobs != config_.observation_level()) {
        throw std::runtime_error("Detector observation level does not match configured value.");
      }
      detector_[ipropoagator.detector0+isphere].sphere = spheres[isphere];
      detector_[ipropoagator.detector0+isphere].nrays_to_refract = 0;
      detector_[ipropoagator.detector0+isphere].nrays_to_propagate = 0;
    }
  }

  return VCLIACTTrackVisitor<VCLArchitecture>::visit_event(event, kill_event);
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
propagate_rays(calin::math::ray::VCLRay<double_real> ray, double_bvt ray_mask, double_vt ray_weight)
{
  ray_mask &= (ray.z()>zobs_) & (ray.uz()<0);

  switch(double_vt dz = ray.z()-zobs_; config_.refraction_mode()) {
  case calin::ix::simulation::vcl_iact::REFRACT_NO_RAYS:
    // Note use of "-=" below because uz() is negative
    ray.mutable_ct() -= vcl::select(ray_mask,
      (dz*ref_index_correction_ + this->atm_->template vcl_propagation_ct_correction_to_iobs<VCLArchitecture>(ray.z(), config_.observation_level()))/ray.uz(),
      0);
    break;
  case calin::ix::simulation::vcl_iact::REFRACT_ALL_RAYS:
    this->atm_->template vcl_propagate_ray_with_refraction_and_mask<VCLArchitecture>(ray, ray_mask, config_.observation_level());
    // Note propagation backwards by distance since uz() is negative
    ray.propagate_dist_with_mask(ray_mask, dz/ray.uz(), ref_index_);
    break;
  case calin::ix::simulation::vcl_iact::REFRACT_ONLY_CLOSE_RAYS:
  default:
    break;
  }

  RayArray ray_array { ray };
  double_at ray_weight_array;
  ray_weight.store(ray_weight_array);

  for(auto& idetector : detector_) {
    auto intersecting_rays = ray_mask & (ray.squared_distance_at_closest_approach(idetector.sphere.r0.template cast<double_vt>()) < idetector.squared_safety_radius);
    unsigned intersecting_rays_bitmask = vcl::to_bits(intersecting_rays);
    if(intersecting_rays_bitmask) {
      for(unsigned iray=0; iray<VCLArchitecture::num_double; ++iray) {
        if(intersecting_rays_bitmask & 1) {
          switch(config_.refraction_mode()) {
          case calin::ix::simulation::vcl_iact::REFRACT_NO_RAYS:
          case calin::ix::simulation::vcl_iact::REFRACT_ALL_RAYS:
            idetector.rays_to_propagate.insert_one_ray(idetector.nrays_to_propagate, ray_array.extract_one_ray(iray));
            idetector.ray_weights_to_propagate[idetector.nrays_to_propagate] = ray_weight_array[iray];
            ++idetector.nrays_to_propagate;
            if(idetector.nrays_to_propagate == VCLArchitecture::num_double) {
              do_propagate_rays_for_detector(idetector);
            }
            break;
          case calin::ix::simulation::vcl_iact::REFRACT_ONLY_CLOSE_RAYS:
          default:
            idetector.rays_to_refract.insert_one_ray(idetector.nrays_to_refract, ray_array.extract_one_ray(iray));
            idetector.ray_weights_to_refract[idetector.nrays_to_refract] = ray_weight_array[iray];
            ++idetector.nrays_to_refract;
            if(idetector.nrays_to_refract == VCLArchitecture::num_double) {
              do_refract_rays_for_detector(idetector);
            }
            break;
          }
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

  double_vt dz = ray.z()-zobs_;
  this->atm_->template vcl_propagate_ray_with_refraction_and_mask<VCLArchitecture>(ray, ray_mask, config_.observation_level());
  // Note propagation backwards by distance since uz() is negative
  ray.propagate_dist_with_mask(ray_mask, dz/ray.uz(), ref_index_);

  idetector.rays_to_refract.set_rays(ray);

  auto intersecting_rays = ray_mask & (ray.squared_distance_at_closest_approach(idetector.sphere.r0.template cast<double_vt>()) < idetector.squared_radius);
  unsigned intersecting_rays_bitmask = vcl::to_bits(intersecting_rays);
  if(intersecting_rays_bitmask) {
    for(unsigned iray=0; iray<VCLArchitecture::num_double; ++iray) {
      if(intersecting_rays_bitmask & 1) {
        idetector.rays_to_propagate.insert_one_ray(idetector.nrays_to_propagate, idetector.rays_to_refract.extract_one_ray(iray));
        idetector.ray_weights_to_propagate[idetector.nrays_to_propagate] = idetector.ray_weights_to_refract[iray];
        ++idetector.nrays_to_propagate;
        if(idetector.nrays_to_propagate == VCLArchitecture::num_double) {
          do_propagate_rays_for_detector(idetector);
        }
      }
      intersecting_rays_bitmask >>= 1;
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

  FocalPlaneParameters fp_parameters;
  ray_mask = idetector.propagator->propagate_rays_to_focal_plane(
    idetector.propagator_iscope, ray, ray_mask, fp_parameters);

  ray_mask &= this->rng_->uniform_double() < fp_parameters.detection_prob;

  unsigned fp_rays_bitmask = vcl::to_bits(ray_mask);
  if(fp_rays_bitmask) {
    double_at fplane_x;
    double_at fplane_y;
    double_at fplane_ux;
    double_at fplane_uy;
    double_at fplane_t;
    int64_at pixel_id;

    fp_parameters.fplane_x.store(fplane_x);
    fp_parameters.fplane_z.store(fplane_y);
    fp_parameters.fplane_ux.store(fplane_ux);
    fp_parameters.fplane_uz.store(fplane_uy);
    fp_parameters.fplane_t.store(fplane_t);
    fp_parameters.pixel_id.store(pixel_id);

    for(unsigned iray=0; iray<VCLArchitecture::num_double; ++iray) {
      if(fp_rays_bitmask & 1) {
        idetector.pe_processor->process_focal_plane_hit(idetector.propagator_iscope,
          pixel_id[iray], fplane_x[iray], fplane_y[iray], fplane_ux[iray], fplane_uy[iray],
          fplane_t[iray], idetector.ray_weights_to_propagate[iray]);
      }
    }
    fp_rays_bitmask >>= 1;
  }
}

template<typename VCLArchitecture> calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration
VCLIACTArray<VCLArchitecture>::default_config()
{
  calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration config;
  config.set_detector_energy_lo(1.5);
  config.set_detector_energy_hi(4.8);
  config.set_detector_energy_bin_width(0.05);
  return config;
}

template<typename VCLArchitecture> calin::ix::simulation::vcl_iact::VCLIACTConfiguration
VCLIACTArray<VCLArchitecture>::base_config(const calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration& config)
{
  calin::ix::simulation::vcl_iact::VCLIACTConfiguration bconfig;
  bconfig.set_bandwidth(config.detector_energy_hi()-config.detector_energy_lo());
  bconfig.set_enable_forced_cherenkov_angle_mode(config.enable_forced_cherenkov_angle_mode());
  bconfig.set_forced_cherenkov_angle(config.forced_cherenkov_angle());
  return bconfig;
}

template<typename VCLArchitecture> std::vector<double>
VCLIACTArray<VCLArchitecture>::detector_efficiency_energy_knots()
{
  std::vector<double> knots;
  for(double e=config_.detector_energy_lo(); e<=config_.detector_energy_hi();
      e+=config_.detector_energy_bin_width()) {
    knots.push_back(e);
  }
  return knots;
}

#endif // not defined SWIG

} } } // namespace calin::simulation::vcl_iact