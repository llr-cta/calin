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
#include <sstream>

#include <util/log.hpp>
#include <util/string.hpp>
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
    const calin::simulation::detector_efficiency::AtmosphericAbsorption& atm_abs,
    const calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration& config = default_config(),
    calin::math::rng::VCLRNG<VCLArchitecture>* rng = nullptr,
    bool adopt_atm = false, bool adopt_rng = false);
  virtual ~VCLIACTArray();

  void add_propagator(FocalPlaneRayPropagator* propagator, PEProcessor* pe_processor,
    const DetectionEfficiency& detector_efficiency, const std::string& propagator_name = "",
    bool adopt_propagator = false, bool adopt_pe_processor = false);

  DaviesCottonVCLFocalPlaneRayPropagator* add_davies_cotton_propagator(calin::simulation::vs_optics::VSOArray* array,
    PEProcessor* pe_processor, const DetectionEfficiency& detector_efficiency, const std::string& propagator_name = "",
    bool adopt_array = false, bool adopt_pe_processor = false);

  DaviesCottonVCLFocalPlaneRayPropagator* add_davies_cotton_propagator(const ix::simulation::vs_optics::IsotropicDCArrayParameters& param,
    PEProcessor* pe_processor, const DetectionEfficiency& detector_efficiency, const std::string& propagator_name = "",
    bool adopt_pe_processor = false);

  void point_telescope_az_el_phi_deg(unsigned iscope, double az_deg, double el_deg, double phi_deg);
  void point_telescope_az_el_deg(unsigned iscope, double az_deg, double el_deg);

  void point_all_telescopes_az_el_phi_deg(const Eigen::VectorXd& az_deg,
    const Eigen::VectorXd& el_deg, const Eigen::VectorXd& phi_deg);
  void point_all_telescopes_az_el_deg(const Eigen::VectorXd& az_deg,
    const Eigen::VectorXd& el_deg);

  void point_all_telescopes_az_el_phi_deg(double az_deg, double el_deg, double phi_deg);
  void point_all_telescopes_az_el_deg(double az_deg, double el_deg);

  static calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration default_config();

  const calin::math::spline_interpolation::CubicMultiSpline& detector_efficiency_spline() const {
    return detector_efficiency_spline_;
  }

  unsigned num_scopes() const { return detector_.size(); }

  std::string banner() const;

  calin::math::spline_interpolation::CubicSpline* new_height_dependent_pe_bandwidth_spline() const;
  double fixed_pe_bandwidth() const;

#ifndef SWIG
  void visit_event(const calin::simulation::tracker::Event& event, bool& kill_event) final;
  void propagate_rays(calin::math::ray::VCLRay<double_real> ray, double_bvt ray_mask,
    double_vt bandwidth, double_vt ray_weight) final;
  // void leave_event() final;

protected:
  using VCLIACTTrackVisitor<VCLArchitecture>::set_fixed_pe_bandwidth_mode;
  using VCLIACTTrackVisitor<VCLArchitecture>::set_fixed_photon_bandwidth_mode;
  using VCLIACTTrackVisitor<VCLArchitecture>::set_height_dependent_pe_bandwidth_mode;

  void update_detector_efficiencies();
  static calin::ix::simulation::vcl_iact::VCLIACTConfiguration base_config(
    const calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration& config);
  std::vector<double> detector_efficiency_energy_knots();
  unsigned add_detector_efficiency(const DetectionEfficiency& detector_efficiency, const std::string& name);

  struct PropagatorInfo {
    FocalPlaneRayPropagator* propagator;
    unsigned ipropagator;
    PEProcessor* pe_processor;
    unsigned detector0;
    unsigned ndetector;
    bool adopt_propagator;
    bool adopt_pe_processor;
    std::string name;
  };

  struct DetectorInfo {
    RayProcessorDetectorSphere sphere;
    double squared_radius;
    double squared_safety_radius;
    FocalPlaneRayPropagator* propagator;
    unsigned ipropagator;
    unsigned propagator_iscope;
    unsigned global_iscope;
    calin::simulation::pe_processor::PEProcessor* pe_processor;
    unsigned idetector_efficiency;

    RayArray rays_to_refract;
    unsigned nrays_to_refract;
    double_at ray_weights_to_refract;
    double_at bandwidths_to_refract;

    RayArray rays_to_propagate;
    unsigned nrays_to_propagate;
    double_at ray_weights_to_propagate;
    double_at bandwidths_to_propagate;
};

  void do_propagate_rays_for_detector(DetectorInfo& idetector);
  void do_refract_rays_for_detector(DetectorInfo& idetector);

  calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration config_;
  calin::simulation::detector_efficiency::AtmosphericAbsorption atm_abs_;
  std::vector<PropagatorInfo> propagator_;
  std::vector<DetectorInfo> detector_;
  std::vector<DetectionEfficiency> detector_efficiency_;
  calin::math::spline_interpolation::CubicMultiSpline detector_efficiency_spline_;
  std::vector<calin::math::spline_interpolation::TwoDimensionalCubicSpline*> detector_bandwidth_spline_;
  double zobs_;
  double wmax_ = 1.0;
  double wmin_ = 0.0;
  double ref_index_;
  double ref_index_correction_;
  double safety_radius_;
#endif
};

#ifndef SWIG

template<typename VCLArchitecture> VCLIACTArray<VCLArchitecture>::
VCLIACTArray(
    calin::simulation::atmosphere::LayeredRefractiveAtmosphere* atm,
    const calin::simulation::detector_efficiency::AtmosphericAbsorption& atm_abs,
    const calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration& config,
    calin::math::rng::VCLRNG<VCLArchitecture>* rng,
    bool adopt_atm, bool adopt_rng):
  VCLIACTTrackVisitor<VCLArchitecture>(atm, base_config(config), rng, adopt_atm, adopt_rng),
  config_(config), atm_abs_(atm_abs), detector_efficiency_spline_(detector_efficiency_energy_knots())
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
  for(auto* ispline : detector_bandwidth_spline_) {
    delete ispline;
  }
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
add_propagator(FocalPlaneRayPropagator* propagator, PEProcessor* pe_processor,
  const DetectionEfficiency& detector_efficiency, const std::string& propagator_name,
  bool adopt_propagator, bool adopt_pe_processor)
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
  propagator_info.name               = propagator_name;
  propagator_.emplace_back(propagator_info);

  std::string name = propagator_name;
  if(name.empty()) {
    name = "propagator "+std::to_string(propagator_info.ipropagator);
  }

  unsigned iefficiency = add_detector_efficiency(detector_efficiency, name);

  for(unsigned isphere=0; isphere<sphere.size(); ++isphere) {
    DetectorInfo detector_info;
    detector_info.sphere                 = sphere[isphere];
    detector_info.squared_radius         = SQR(sphere[isphere].radius);
    detector_info.squared_safety_radius  = SQR(sphere[isphere].radius + safety_radius_);
    detector_info.propagator             = propagator;
    detector_info.ipropagator            = propagator_.size()-1;
    detector_info.propagator_iscope      = isphere;
    detector_info.global_iscope          = detector_.size();
    detector_info.pe_processor           = pe_processor;
    detector_info.idetector_efficiency   = iefficiency;
    detector_info.nrays_to_refract       = 0;
    detector_info.nrays_to_propagate     = 0;
    detector_.emplace_back(detector_info);
  }

  update_detector_efficiencies();
}

template<typename VCLArchitecture>
calin::simulation::vcl_ray_propagator::DaviesCottonVCLFocalPlaneRayPropagator<VCLArchitecture>*
VCLIACTArray<VCLArchitecture>::add_davies_cotton_propagator(
  calin::simulation::vs_optics::VSOArray* array,
  PEProcessor* pe_processor, const DetectionEfficiency& detector_efficiency,
  const std::string& propagator_name,
  bool adopt_array, bool adopt_pe_processor)
{
  auto* propagator = new calin::simulation::vcl_ray_propagator::DaviesCottonVCLFocalPlaneRayPropagator<VCLArchitecture>(
    array, this->rng_, ref_index_, adopt_array, /* adopt_rng= */ false);
  add_propagator(propagator, pe_processor, detector_efficiency, propagator_name,
    /* adopt_propagator= */ true, adopt_pe_processor);
  return propagator;
}

template<typename VCLArchitecture>
calin::simulation::vcl_ray_propagator::DaviesCottonVCLFocalPlaneRayPropagator<VCLArchitecture>*
VCLIACTArray<VCLArchitecture>::add_davies_cotton_propagator(
  const ix::simulation::vs_optics::IsotropicDCArrayParameters& param,
  PEProcessor* pe_processor, const DetectionEfficiency& detector_efficiency,
  const std::string& propagator_name,
  bool adopt_pe_processor)
{
  auto* array = new calin::simulation::vs_optics::VSOArray;
  calin::math::rng::VCLToScalarRNGCore scalar_core(this->rng_->core());
  calin::math::rng::RNG scalar_rng(&scalar_core);
  array->generateFromArrayParameters(param, scalar_rng);
  return add_davies_cotton_propagator(array, pe_processor, detector_efficiency, propagator_name,
    /* adopt_array= */ true, adopt_pe_processor);
}

template<typename VCLArchitecture> unsigned VCLIACTArray<VCLArchitecture>::
add_detector_efficiency(const DetectionEfficiency& detector_efficiency, const std::string& name)
{
  unsigned iefficiency = detector_efficiency_.size();

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

  detector_bandwidth_spline_.push_back(
    atm_abs_.integrate_bandwidth_to_spline(zobs_, detector_efficiency));

  return iefficiency;
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
point_telescope_az_el_phi_deg(unsigned iscope,
  double az_deg, double el_deg, double phi_deg)
{
  if(iscope >= detector_.size()) {
    throw std::out_of_range("Telescope ID out of range");
  }

  DetectorInfo& idetector(detector_[iscope]);
  PropagatorInfo& ipropagator(propagator_[idetector.ipropagator]);
  unsigned propagator_isphere = iscope-ipropagator.detector0;

  ipropagator.propagator->point_telescope_az_el_phi_deg(
    propagator_isphere, az_deg, el_deg, phi_deg);

  update_detector_efficiencies();
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
point_telescope_az_el_deg(unsigned iscope,  double az_deg, double el_deg)
{
  this->point_telescope_az_el_phi_deg(iscope, az_deg, el_deg, 0.0);
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
point_all_telescopes_az_el_phi_deg(const Eigen::VectorXd& az_deg,
  const Eigen::VectorXd&  el_deg, const Eigen::VectorXd&  phi_deg)
{
  for(auto& ipropagator : propagator_) {
    for(unsigned propagator_isphere=0; propagator_isphere<ipropagator.ndetector;
        ++propagator_isphere) {
      unsigned isphere = ipropagator.detector0 + propagator_isphere;
      ipropagator.propagator->point_telescope_az_el_phi_deg(
        propagator_isphere, az_deg[isphere], el_deg[isphere],
        phi_deg[isphere]);
    }
  }
  update_detector_efficiencies();
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
point_all_telescopes_az_el_deg(const Eigen::VectorXd& az_deg, const Eigen::VectorXd& el_deg)
{
  point_all_telescopes_az_el_phi_deg(az_deg, el_deg,
    Eigen::VectorXd::Constant(detector_.size(), 0.0));
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
point_all_telescopes_az_el_phi_deg(double az_deg, double el_deg, double phi_deg)
{
  point_all_telescopes_az_el_phi_deg(
    Eigen::VectorXd::Constant(detector_.size(), az_deg),
    Eigen::VectorXd::Constant(detector_.size(), el_deg),
    Eigen::VectorXd::Constant(detector_.size(), phi_deg));
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
point_all_telescopes_az_el_deg(double az_deg, double el_deg)
{
  point_all_telescopes_az_el_phi_deg(
    Eigen::VectorXd::Constant(detector_.size(), az_deg),
    Eigen::VectorXd::Constant(detector_.size(), el_deg),
    Eigen::VectorXd::Constant(detector_.size(), 0.0));
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
update_detector_efficiencies()
{
  using calin::math::special::SQR;

  double znmin = M_PI_2;
  double znmax = 0;
  for(auto& ipropagator : propagator_) {
    auto spheres = ipropagator.propagator->detector_spheres();
    if(spheres.size() != ipropagator.ndetector) {
      // this should never happen
      throw std::runtime_error("Number of detectors proposed by propagator must remain constant over events.");
    }
    for(unsigned propagator_isphere=0; propagator_isphere<ipropagator.ndetector;
      ++propagator_isphere) {
      const auto& isphere(spheres[propagator_isphere]);
      auto& idetector(detector_[ipropagator.detector0 + propagator_isphere]);
      if(isphere.iobs != config_.observation_level()) {
        throw std::runtime_error("Detector observation level does not match configured value.");
      }
      idetector.sphere = isphere;
      double zn = std::atan2(std::sqrt(SQR(isphere.obs_dir.x())+SQR(isphere.obs_dir.y())), isphere.obs_dir.z());
      znmin = std::min(znmin, std::max(zn - isphere.field_of_view_radius, 0.0));
      znmax = std::max(znmax, std::min(zn + isphere.field_of_view_radius, M_PI_2));
    }
    wmax_ = std::cos(znmin);
    wmin_ = std::cos(znmax);
    safety_radius_ = this->atm_->refraction_safety_radius(znmax, config_.observation_level());
    for(unsigned propagator_isphere=0; propagator_isphere<ipropagator.ndetector;
      ++propagator_isphere) {
      auto& idetector(detector_[ipropagator.detector0 + propagator_isphere]);
      const auto& isphere(idetector.sphere);
      idetector.squared_radius         = SQR(isphere.radius);
      idetector.squared_safety_radius  = SQR(isphere.radius + safety_radius_);
    }
  }

  switch(config_.cherenkov_mode()) {
  case calin::ix::simulation::vcl_iact::PHOTON_MODE:
  default:
    this->set_fixed_photon_bandwidth_mode(
      config_.detector_energy_hi()-config_.detector_energy_lo(),
      config_.detector_energy_lo());
    break;
  case calin::ix::simulation::vcl_iact::FIXED_BANDWIDTH_PE_MODE:
    this->set_fixed_pe_bandwidth_mode(fixed_pe_bandwidth());
    break;
  case calin::ix::simulation::vcl_iact::VARIABLE_BANDWIDTH_PE_MODE:
    this->set_height_dependent_pe_bandwidth_mode(new_height_dependent_pe_bandwidth_spline(), true);
    break;
  }
}

template<typename VCLArchitecture> double VCLIACTArray<VCLArchitecture>::
fixed_pe_bandwidth() const
{
  double bandwidth = 0;
  for(unsigned ispline=0; ispline<detector_efficiency_spline_.num_spline(); ++ispline) {
    bandwidth = std::max(bandwidth, detector_efficiency_spline_.integral(
      detector_efficiency_spline_.xmax(), ispline));
  }
  return bandwidth;
}

template<typename VCLArchitecture> calin::math::spline_interpolation::CubicSpline*
VCLIACTArray<VCLArchitecture>::new_height_dependent_pe_bandwidth_spline() const
{
  if(detector_bandwidth_spline_.size() == 0) {
    return nullptr;
  }
  std::vector<double> heights = detector_bandwidth_spline_[0]->xknot_as_stdvec();
  std::vector<double> bandwidths(heights.size(), 0.0);
  for(unsigned ispline=0; ispline<detector_bandwidth_spline_.size(); ++ispline) {
    for(unsigned iheight=0; iheight<heights.size(); ++iheight) {
      bandwidths[iheight] = std::max(bandwidths[iheight],
        detector_bandwidth_spline_[ispline]->value(heights[iheight], wmax_));
    }
  }
  return new calin::math::spline_interpolation::CubicSpline(heights, bandwidths);
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
propagate_rays(calin::math::ray::VCLRay<double_real> ray, double_bvt ray_mask,
  double_vt bandwidth, double_vt ray_weight)
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
  double_at bandwidth_array;
  double_at ray_weight_array;
  bandwidth.store(bandwidth_array);
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
            idetector.bandwidths_to_propagate[idetector.nrays_to_propagate] = bandwidth_array[iray];
            idetector.ray_weights_to_propagate[idetector.nrays_to_propagate] = ray_weight_array[iray];
            ++idetector.nrays_to_propagate;
            if(idetector.nrays_to_propagate == VCLArchitecture::num_double) {
              do_propagate_rays_for_detector(idetector);
            }
            break;
          case calin::ix::simulation::vcl_iact::REFRACT_ONLY_CLOSE_RAYS:
          default:
            idetector.rays_to_refract.insert_one_ray(idetector.nrays_to_refract, ray_array.extract_one_ray(iray));
            idetector.bandwidths_to_refract[idetector.nrays_to_refract] = bandwidth_array[iray];
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
        idetector.bandwidths_to_propagate[idetector.nrays_to_propagate] = idetector.bandwidths_to_refract[iray];
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
  config.set_cherenkov_mode(calin::ix::simulation::vcl_iact::VARIABLE_BANDWIDTH_PE_MODE);
  config.set_detector_energy_lo(1.25);
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

template<typename VCLArchitecture> std::string VCLIACTArray<VCLArchitecture>::banner() const
{
  constexpr double EV_NM = 1239.84193009239; // gunits: c/(ev/h) -> nm
  using calin::util::string::double_to_string_with_commas;
  std::ostringstream stream;
  stream
    << "Class : " << calin::util::vcl::templated_class_name<VCLArchitecture>("VCLIACTArray") << '\n'
    << "Number of focal-plane propagators : " << propagator_.size() << ", with "
    << detector_.size() << " detectors.\n"
    << "Detector zenith range : " << double_to_string_with_commas(std::acos(wmax_)/M_PI*180.0,1)
    << " to " << double_to_string_with_commas(std::acos(wmin_)/M_PI*180.0,1) << " degrees.\n"
    << "Observation level : " << double_to_string_with_commas(zobs_/1e5,3) << " km, refraction safety radius : "
    << double_to_string_with_commas(safety_radius_/100,2) << " m.\n";
  if(this->variable_bandwidth_spline_) {
    stream << "Cherenkov ray mode : PEs, with height-dependent bandwidth\n"
      << "- " << double_to_string_with_commas(zobs_/1e5,3) << ", 10, 20, 30, "
      << double_to_string_with_commas(this->atm_->top_of_atmosphere()/1e5, 0) << " km : "
      << double_to_string_with_commas(this->variable_bandwidth_spline_->value(zobs_),3) << ", "
      << double_to_string_with_commas(this->variable_bandwidth_spline_->value(10e5),3) << ", "
      << double_to_string_with_commas(this->variable_bandwidth_spline_->value(20e5),3) << ", "
      << double_to_string_with_commas(this->variable_bandwidth_spline_->value(30e5),3) << ", "
      << double_to_string_with_commas(this->variable_bandwidth_spline_->value(this->atm_->top_of_atmosphere()),3) << " eV\n";
  } else if (this->do_color_photons_) {
    stream << "Cherenkov ray mode : photons, with bandwidth "
      << double_to_string_with_commas(this->fixed_bandwidth_,3) << " eV\n"
      << "- Energy range "
      << double_to_string_with_commas(this->min_cherenkov_energy_,3) << " - "
      << double_to_string_with_commas(this->min_cherenkov_energy_+this->fixed_bandwidth_,3) << " eV ("
      << double_to_string_with_commas(EV_NM/(this->min_cherenkov_energy_+this->fixed_bandwidth_),0) << " - "
      << double_to_string_with_commas(EV_NM/this->min_cherenkov_energy_,0) << " nm)\n";
  } else {
    stream << "Cherenkov ray mode : PEs, with fixed bandwidth "
      << double_to_string_with_commas(this->fixed_bandwidth_,3) << " eV\n";
  }
  if(detector_efficiency_spline_.num_spline() > 0) {
    stream << "Detector efficiency bandwidths :\n";
    for(unsigned ispline=0; ispline<detector_efficiency_spline_.num_spline(); ++ispline) {
      stream
        << "- " << detector_efficiency_spline_.dataset_name(ispline) << " : "
        << double_to_string_with_commas(detector_efficiency_spline_.integral(detector_efficiency_spline_.xmax(), ispline),3) << " eV\n"
        << "  Absorbed from 10 km : " << double_to_string_with_commas(detector_bandwidth_spline_[ispline]->value(10e5,wmin_),3)
        << " to " << double_to_string_with_commas(detector_bandwidth_spline_[ispline]->value(10e5,wmax_),3) << " eV\n"
        << "  Absorbed from 20 km : " << double_to_string_with_commas(detector_bandwidth_spline_[ispline]->value(20e5,wmin_),3)
        << " to " << double_to_string_with_commas(detector_bandwidth_spline_[ispline]->value(20e5,wmax_),3) << " eV\n";
    }
  }
  return stream.str();
}

#endif // not defined SWIG

} } } // namespace calin::simulation::vcl_iact
