/*

   calin/simulation/vcl_iact_array.hpp -- Stephen Fegan -- 2022-06-29

   Class for imaging atmospheric cherenkov technique - produce rays from charged
   tracks in the atmosphere, propagate them to the ground and trace them through
   telescope optics.

   This has become a very complex class that supports multiple propagators with 
   multiple telescopes, each with different detector efficiency and angular
   response. It also now supports multiple propagator sets, each with a different
   scattering radius. IT IS NOT A GOOD EXAMPLE OF HOW TO WRITE CODE.

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
#include <math/hex_array.hpp>
#include <simulation/vcl_iact.hpp>
#include <simulation/vcl_ray_propagator.hpp>
#include <simulation/vcl_raytracer.hpp>
#include <simulation/detector_efficiency.hpp>
#include <simulation/vcl_detector_efficiency.hpp>
#include <simulation/atmosphere.hpp>

namespace calin { namespace simulation { namespace vcl_iact {

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCLBandwidthManager
{
public:
#ifndef SWIG
  using double_vt   = typename VCLArchitecture::double_vt;
#endif

  VCLBandwidthManager(const std::string& name): name_(name) {
    // nothing to see here
  }

  virtual ~VCLBandwidthManager() {
    // nothing to see here
  }

  virtual double bandwidth() const {
    return 0;
  }

  virtual std::vector<double> bandwidth_vs_height(const std::vector<double>& h, double w) const {
    return std::vector<double>(h.size(), 0.0);
  }

  virtual const calin::math::spline_interpolation::CubicSpline* detector_efficiency_spline() const
  {
    return nullptr;
  }

  virtual const calin::math::spline_interpolation::TwoDimensionalCubicSpline* detector_bandwidth_spline() const
  {
    return nullptr;
  }

  virtual const calin::simulation::detector_efficiency::VCLDirectionResponse<VCLArchitecture>* fp_angular_response() const {
    return nullptr;
  }

  virtual std::string banner(double wmin, double wmax, const std::string& indent_1 = "", const std::string& indent_n = "") const {
    return indent_1 + name_;
  }

#ifndef SWIG
  virtual double_vt bandwidth_for_pe(const double_vt z_emission, const double_vt uz_emission,
    const double_vt ux_fp, const double_vt uy_fp, const double_vt uz_fp) const
  {
    return 0;
  }
#endif

protected:
  static calin::math::spline_interpolation::CubicSpline* new_detector_efficiency_spline(
    const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency,
    double e_lo, double e_hi, double delta_e)
  {
    std::vector<double> xknot;
    for(double e=e_lo; e<e_hi; e+=delta_e) { xknot.push_back(e); }

    double full_integral = detector_efficiency.integrate();
    double partial_integral = detector_efficiency.integrate(e_lo, e_hi);

    if(partial_integral < 0.99*full_integral) {
      calin::util::log::LOG(calin::util::log::WARNING)
        << "Configured detector energy range does contains less than 99% of detector efficiency curve.";
    }

    std::vector<double> yknot(xknot.size());
    std::transform(xknot.begin(), xknot.end(), yknot.begin(),
      [detector_efficiency](double e) { return detector_efficiency(e); });

    return new calin::math::spline_interpolation::CubicSpline(xknot, yknot);
  }

  static calin::math::spline_interpolation::TwoDimensionalCubicSpline* new_detector_bandwidth_spline(
    const calin::simulation::detector_efficiency::DetectionEfficiency& detector_efficiency,
    const calin::simulation::detector_efficiency::AtmosphericAbsorption& atm_abs,
    double zobs)
  {
    return atm_abs.integrate_bandwidth_to_spline(zobs, detector_efficiency);
  }

  std::string name_;
};

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCLSimpleBandwidthManager:
  public VCLBandwidthManager<VCLArchitecture>
{
public:
#ifndef SWIG
  using double_vt   = typename VCLArchitecture::double_vt;
#endif

  VCLSimpleBandwidthManager(
      const calin::simulation::detector_efficiency::AtmosphericAbsorption* atm_abs,
      const calin::simulation::detector_efficiency::DetectionEfficiency& detection_efficiency,
      double zobs, double e_lo, double e_hi, double delta_e, const std::string& name):
    VCLBandwidthManager<VCLArchitecture>(name), atm_abs_(atm_abs),
    detection_efficiency_(detection_efficiency)
  {
    detector_efficiency_spline_ = VCLBandwidthManager<VCLArchitecture>::
      new_detector_efficiency_spline(detection_efficiency_, e_lo, e_hi, delta_e);
    detector_bandwidth_spline_ = VCLBandwidthManager<VCLArchitecture>::
      new_detector_bandwidth_spline(detection_efficiency_, *atm_abs_, zobs);
  }

  virtual ~VCLSimpleBandwidthManager() {
    delete detector_efficiency_spline_;
    delete detector_bandwidth_spline_;
  }

  double bandwidth() const final {
    return detector_efficiency_spline_->integral(detector_efficiency_spline_->xmax());
  }

  std::vector<double> bandwidth_vs_height(const std::vector<double>& h, double w) const final {
    std::vector<double> bw(h.size());
    std::transform(h.begin(), h.end(), bw.begin(),
      [this,w](double hh) { return detector_bandwidth_spline_->value(hh,w); });
    return bw;
  }

  const calin::math::spline_interpolation::CubicSpline* detector_efficiency_spline() const final {
    return detector_efficiency_spline_;
  }

  const calin::math::spline_interpolation::TwoDimensionalCubicSpline* detector_bandwidth_spline() const final {
    return detector_bandwidth_spline_;
  }

  virtual std::string banner(double wmin, double wmax, const std::string& indent_1 = "", const std::string& indent_n = "") const {
    // constexpr double EV_NM = 1239.84193009239; // gunits: c/(ev/h) -> nm
    using calin::util::string::double_to_string_with_commas;
    std::ostringstream stream;
    stream
      << indent_1 << this->name_ << " : "
      << double_to_string_with_commas(detector_efficiency_spline_->integral(detector_efficiency_spline_->xmax()),3) << " eV\n"
      << indent_n << "Absorbed from 5 km : " << double_to_string_with_commas(detector_bandwidth_spline_->value(5e5,wmin),3)
      << " to " << double_to_string_with_commas(detector_bandwidth_spline_->value(5e5,wmax),3) << " eV\n"
      << indent_n << "Absorbed from 10 km : " << double_to_string_with_commas(detector_bandwidth_spline_->value(10e5,wmin),3)
      << " to " << double_to_string_with_commas(detector_bandwidth_spline_->value(10e5,wmax),3) << " eV\n"
      << indent_n << "Absorbed from 15 km : " << double_to_string_with_commas(detector_bandwidth_spline_->value(15e5,wmin),3)
      << " to " << double_to_string_with_commas(detector_bandwidth_spline_->value(15e5,wmax),3) << " eV";
    return stream.str();
  }

#ifndef SWIG
  double_vt bandwidth_for_pe(const double_vt z_emission, const double_vt uz_emission,
    const double_vt ux_fp, const double_vt uy_fp, const double_vt uz_fp) const final
  {
    return detector_bandwidth_spline_->
      template vcl_value<VCLArchitecture>(z_emission, vcl::abs(uz_emission));
  }
#endif

private:
  const calin::simulation::detector_efficiency::AtmosphericAbsorption* atm_abs_ = nullptr;
  calin::simulation::detector_efficiency::DetectionEfficiency detection_efficiency_;
  calin::math::spline_interpolation::CubicSpline* detector_efficiency_spline_ = nullptr;
  calin::math::spline_interpolation::TwoDimensionalCubicSpline* detector_bandwidth_spline_ = nullptr;
};

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCLDCBandwidthManager:
  public VCLBandwidthManager<VCLArchitecture>
{
public:
#ifndef SWIG
  using double_vt   = typename VCLArchitecture::double_vt;
#endif

  VCLDCBandwidthManager(
      const calin::simulation::detector_efficiency::AtmosphericAbsorption* atm_abs,
      const calin::simulation::detector_efficiency::DetectionEfficiency& detection_efficiency,
      const calin::simulation::detector_efficiency::AngularEfficiency& fp_angular_efficiency,
      double zobs, double e_lo, double e_hi, double delta_e, const std::string& name):
    VCLBandwidthManager<VCLArchitecture>(name), atm_abs_(atm_abs),
    detection_efficiency_(detection_efficiency)
  {
    fp_angular_response_ = new calin::simulation::detector_efficiency::
      VCLUY1DSplineDirectionResponse<VCLArchitecture>(fp_angular_efficiency, true, 1.0);
    detection_efficiency_ *= fp_angular_response_->scale();
    detector_efficiency_spline_ = VCLBandwidthManager<VCLArchitecture>::
      new_detector_efficiency_spline(detection_efficiency_, e_lo, e_hi, delta_e);
    detector_bandwidth_spline_ = VCLBandwidthManager<VCLArchitecture>::
      new_detector_bandwidth_spline(detection_efficiency_, *atm_abs_, zobs);
  }

  virtual ~VCLDCBandwidthManager() {
    delete detector_efficiency_spline_;
    delete detector_bandwidth_spline_;
    delete fp_angular_response_;
  }

  double bandwidth() const final {
    return detector_efficiency_spline_->integral(detector_efficiency_spline_->xmax());
  }

  std::vector<double> bandwidth_vs_height(const std::vector<double>& h, double w) const final {
    std::vector<double> bw(h.size());
    std::transform(h.begin(), h.end(), bw.begin(),
      [this,w](double hh) { return detector_bandwidth_spline_->value(hh,w); });
    return bw;
  }

  const calin::math::spline_interpolation::CubicSpline* detector_efficiency_spline() const final {
    return detector_efficiency_spline_;
  }

  const calin::math::spline_interpolation::TwoDimensionalCubicSpline* detector_bandwidth_spline() const final {
    return detector_bandwidth_spline_;
  }

  const calin::simulation::detector_efficiency::VCLDirectionResponse<VCLArchitecture>* fp_angular_response() const final {
    return fp_angular_response_;
  }

  virtual std::string banner(double wmin, double wmax, const std::string& indent_1 = "", const std::string& indent_n = "") const {
    // constexpr double EV_NM = 1239.84193009239; // gunits: c/(ev/h) -> nm
    using calin::util::string::double_to_string_with_commas;
    std::ostringstream stream;
    stream
      << indent_1 << this->name_ << " : "
      << double_to_string_with_commas(detector_efficiency_spline_->integral(detector_efficiency_spline_->xmax()),3) << " eV\n"
      << indent_n << "Absorbed from 5 km : " << double_to_string_with_commas(detector_bandwidth_spline_->value(5e5,wmin),3)
      << " to " << double_to_string_with_commas(detector_bandwidth_spline_->value(5e5,wmax),3) << " eV\n"
      << indent_n << "Absorbed from 10 km : " << double_to_string_with_commas(detector_bandwidth_spline_->value(10e5,wmin),3)
      << " to " << double_to_string_with_commas(detector_bandwidth_spline_->value(10e5,wmax),3) << " eV\n"
      << indent_n << "Absorbed from 15 km : " << double_to_string_with_commas(detector_bandwidth_spline_->value(15e5,wmin),3)
      << " to " << double_to_string_with_commas(detector_bandwidth_spline_->value(15e5,wmax),3) << " eV";
    std::string fp_banner = fp_angular_response_->banner();
    if(not fp_banner.empty()) {
      stream << indent_n << "\nCone : " << fp_banner;
    }
    return stream.str();
  }

#ifndef SWIG
  double_vt bandwidth_for_pe(const double_vt z_emission, const double_vt uz_emission,
    const double_vt ux_fp, const double_vt uy_fp, const double_vt uz_fp) const final
  {
    return detector_bandwidth_spline_->
      template vcl_value<VCLArchitecture>(z_emission, vcl::abs(uz_emission));
  }
#endif

private:
  const calin::simulation::detector_efficiency::AtmosphericAbsorption* atm_abs_ = nullptr;
  calin::simulation::detector_efficiency::DetectionEfficiency detection_efficiency_;
  calin::math::spline_interpolation::CubicSpline* detector_efficiency_spline_ = nullptr;
  calin::math::spline_interpolation::TwoDimensionalCubicSpline* detector_bandwidth_spline_ = nullptr;
  calin::simulation::detector_efficiency::VCLDirectionResponse<VCLArchitecture>* fp_angular_response_ = nullptr;
};

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
  CALIN_TYPEALIAS(AngularEfficiency, calin::simulation::detector_efficiency::AngularEfficiency);
  CALIN_TYPEALIAS(SplinePEAmplitudeGenerator, calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator);

  CALIN_TYPEALIAS(FocalPlaneRayPropagator, calin::simulation::vcl_ray_propagator::VCLFocalPlaneRayPropagator<VCLArchitecture>);
  CALIN_TYPEALIAS(DaviesCottonVCLFocalPlaneRayPropagator, calin::simulation::vcl_ray_propagator::DaviesCottonVCLFocalPlaneRayPropagator<VCLArchitecture>);
  CALIN_TYPEALIAS(PerfectOpticsVCLFocalPlaneRayPropagator, calin::simulation::vcl_ray_propagator::PerfectOpticsVCLFocalPlaneRayPropagator<VCLArchitecture>);
  CALIN_TYPEALIAS(AllSkyVCLFocalPlaneRayPropagator, calin::simulation::vcl_ray_propagator::AllSkyVCLFocalPlaneRayPropagator<VCLArchitecture>);

  VCLIACTArray(calin::simulation::atmosphere::LayeredRefractiveAtmosphere* atm,
    const calin::simulation::detector_efficiency::AtmosphericAbsorption& atm_abs,
    const calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration& config = default_config(),
    calin::math::rng::VCLRNG<VCLArchitecture>* rng = nullptr,
    bool adopt_atm = false, bool adopt_rng = false);

  VCLIACTArray(calin::simulation::atmosphere::LayeredRefractiveAtmosphere* atm,
    const calin::simulation::detector_efficiency::AtmosphericAbsorption& atm_abs,
    const calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration& config,
    unsigned rng_seed, bool adopt_atm = false);

  virtual ~VCLIACTArray();

  unsigned add_propagator_set(double scattering_radius = 0.0, const std::string& name = "");

  DaviesCottonVCLFocalPlaneRayPropagator* add_davies_cotton_propagator(
    calin::simulation::vs_optics::VSOArray* array, PEProcessor* pe_processor,
    const DetectionEfficiency& detector_efficiency, const AngularEfficiency& fp_angular_efficiency,
    SplinePEAmplitudeGenerator* pe_generator = nullptr, const std::string& propagator_name = "",
    bool adopt_array = false, bool adopt_pe_processor = false, bool adopt_pe_generator = false);

  DaviesCottonVCLFocalPlaneRayPropagator* add_davies_cotton_propagator(
    const ix::simulation::vs_optics::IsotropicDCArrayParameters& param, PEProcessor* pe_processor,
    const DetectionEfficiency& detector_efficiency, const AngularEfficiency& fp_angular_efficiency,
    SplinePEAmplitudeGenerator* pe_generator = nullptr, const std::string& propagator_name = "",
    bool adopt_pe_processor = false, bool adopt_pe_generator = false);

  DaviesCottonVCLFocalPlaneRayPropagator* add_davies_cotton_propagator(
    calin::simulation::vs_optics::VSOArray* array, PEProcessor* pe_processor,
    const DetectionEfficiency& detector_efficiency, const AngularEfficiency& fp_angular_efficiency,
    const std::string& propagator_name, SplinePEAmplitudeGenerator* pe_generator = nullptr,
    bool adopt_array = false, bool adopt_pe_processor = false, bool adopt_pe_generator = false);

  DaviesCottonVCLFocalPlaneRayPropagator* add_davies_cotton_propagator(
    const ix::simulation::vs_optics::IsotropicDCArrayParameters& param, PEProcessor* pe_processor,
    const DetectionEfficiency& detector_efficiency, const AngularEfficiency& fp_angular_efficiency,
    const std::string& propagator_name, SplinePEAmplitudeGenerator* pe_generator = nullptr,
    bool adopt_pe_processor = false, bool adopt_pe_generator = false);

  PerfectOpticsVCLFocalPlaneRayPropagator* add_perfect_optics_propagator(
    const Eigen::VectorXd& x, const Eigen::VectorXd& y, const Eigen::VectorXd& z, 
    double radius, double focal_length, double field_of_view_radius, PEProcessor* pe_processor,
    const DetectionEfficiency& detector_efficiency,
    const std::string& propagator_name, SplinePEAmplitudeGenerator* pe_generator = nullptr,
    bool adopt_pe_processor = false, bool adopt_pe_generator = false);

  AllSkyVCLFocalPlaneRayPropagator* add_all_sky_propagator(
    Eigen::VectorXd& r0, double radius, double field_of_view_radius, PEProcessor* pe_processor,
    const DetectionEfficiency& detector_efficiency,
    const std::string& propagator_name, SplinePEAmplitudeGenerator* pe_generator = nullptr,
    bool adopt_pe_processor = false, bool adopt_pe_generator = false);

  void point_telescope_az_el_phi_deg(unsigned iscope, double az_deg, double el_deg, double phi_deg);
  void point_telescope_az_el_deg(unsigned iscope, double az_deg, double el_deg);

  void point_all_telescopes_az_el_phi_deg(const Eigen::VectorXd& az_deg,
    const Eigen::VectorXd& el_deg, const Eigen::VectorXd& phi_deg);
  void point_all_telescopes_az_el_deg(const Eigen::VectorXd& az_deg,
    const Eigen::VectorXd& el_deg);

  void point_all_telescopes_az_el_phi_deg(double az_deg, double el_deg, double phi_deg);
  void point_all_telescopes_az_el_deg(double az_deg, double el_deg);

  static calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration default_config();

  unsigned num_propagator_sets() const { return propagator_set_.size(); }
  double scattered_distance(unsigned ipropagator_set) const {
    return propagator_set_.at(ipropagator_set)->scattered_distance; }
  Eigen::Vector3d scattered_offset(unsigned ipropagator_set) const {
    return propagator_set_.at(ipropagator_set)->scattered_offset; }
  double scattering_radius(unsigned ipropagator_set) const {
    return propagator_set_.at(ipropagator_set)->scattering_radius; }

  unsigned num_propagators() const { return propagator_.size(); }
  unsigned num_scopes() const { return detector_.size(); }

  uint64_t nrays_refracted_at_detector(unsigned idetector) const {
    return detector_.at(idetector)->nrays_refracted; }
  uint64_t nrays_propagated_at_detector(unsigned idetector) const {
    return detector_.at(idetector)->nrays_propagated; }

  std::string banner() const;

  calin::math::spline_interpolation::CubicSpline* new_height_dependent_pe_bandwidth_spline() const;
  double fixed_pe_bandwidth() const;

#ifndef SWIG
  void visit_event(const calin::simulation::tracker::Event& event, bool& kill_event) final;
  void propagate_rays(calin::math::ray::VCLRay<double_real> ray, double_bvt ray_mask,
    double_vt bandwidth, double_vt ray_weight) final;
  void leave_event() final;

protected:
  using VCLIACTTrackVisitor<VCLArchitecture>::set_fixed_pe_bandwidth_mode;
  using VCLIACTTrackVisitor<VCLArchitecture>::set_fixed_photon_bandwidth_mode;
  using VCLIACTTrackVisitor<VCLArchitecture>::set_height_dependent_pe_bandwidth_mode;

  void add_propagator(FocalPlaneRayPropagator* propagator, PEProcessor* pe_processor,
    VCLBandwidthManager<VCLArchitecture>* bandwidth_manager,
    SplinePEAmplitudeGenerator* pe_generator, const std::string& propagator_name,
    bool adopt_propagator, bool adopt_pe_processor, bool adopt_pe_generator);

  void update_detector_efficiencies();
  static calin::ix::simulation::vcl_iact::VCLIACTConfiguration base_config(
    const calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration& config);

  void make_detector_grid();

  struct PropagatorInfo;

  struct PropagatorSet {
    unsigned ipropagator_set;
    double scattering_radius;
    double scattered_distance;
    Eigen::Vector3d scattered_offset;
    std::string name;
    std::vector<PropagatorInfo*> propagators;
  };

  struct DetectorInfo;

  struct PropagatorInfo {
    PropagatorSet* propagator_set;

    FocalPlaneRayPropagator* propagator;
    unsigned ipropagator;
    PEProcessor* pe_processor;
    unsigned detector0;
    unsigned ndetector;
    VCLBandwidthManager<VCLArchitecture>* bandwidth_manager;
    SplinePEAmplitudeGenerator* pe_generator;

    bool adopt_propagator;
    bool adopt_pe_processor;
    bool adopt_pe_generator;
    std::string name;
    std::vector<DetectorInfo*> detector_infos;
  };

  struct DetectorInfo {
    PropagatorInfo* propagator_info;

    RayProcessorDetectorSphere sphere;
    RayProcessorDetectorSphere unscattered_sphere;
    double squared_radius;
    double squared_safety_radius;
    FocalPlaneRayPropagator* propagator;
    unsigned ipropagator;
    unsigned propagator_iscope;
    unsigned global_iscope;
    calin::simulation::pe_processor::PEProcessor* pe_processor;
    VCLBandwidthManager<VCLArchitecture>* bandwidth_manager;
    SplinePEAmplitudeGenerator* pe_generator;

    RayArray rays_to_refract;
    unsigned nrays_to_refract;
    double_at ray_weights_to_refract;
    double_at bandwidths_to_refract;

    RayArray rays_to_propagate;
    unsigned nrays_to_propagate;
    double_at ray_weights_to_propagate;
    double_at bandwidths_to_propagate;

    uint64_t nrays_refracted;
    uint64_t nrays_propagated;
  };

  void do_propagate_rays_for_detector(DetectorInfo* detector);
  void do_refract_rays_for_detector(DetectorInfo* detector);

  calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration config_;
  calin::simulation::detector_efficiency::AtmosphericAbsorption atm_abs_;
  std::vector<PropagatorSet*> propagator_set_;
  std::vector<PropagatorInfo*> propagator_;
  std::vector<DetectorInfo*> detector_;
  std::vector<VCLBandwidthManager<VCLArchitecture>*> bandwidth_manager_;

  double zobs_;
  double wmax_ = 1.0;
  double wmin_ = 0.0;
  double ref_index_;
  double ref_index_correction_;
  double safety_radius_;

  double grid_sep_ = 0.0;
  double grid_sep_inv_ = 0.0;
  unsigned grid_ncells_ = 0;
  unsigned grid_ndetector_per_cell_ = 0;
  unsigned grid_array_size_ = 0;
  double* grid_x = nullptr;
  double* grid_y = nullptr;
  double* grid_z = nullptr;
  double* grid_ssr = nullptr;
  int64_t* grid_idetector = nullptr;

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
  config_(config), atm_abs_(atm_abs)
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
VCLIACTArray(calin::simulation::atmosphere::LayeredRefractiveAtmosphere* atm,
  const calin::simulation::detector_efficiency::AtmosphericAbsorption& atm_abs,
  const calin::ix::simulation::vcl_iact::VCLIACTArrayConfiguration& config,
  unsigned rng_seed, bool adopt_atm):
    VCLIACTArray(atm, atm_abs, config, new calin::math::rng::VCLRNG<VCLArchitecture>(rng_seed,
        __PRETTY_FUNCTION__, "RNG for Cherenkov photon generation and propagation"),
      adopt_atm, /* adopt_rng = */ true)
{
  // nothing to see here
}

template<typename VCLArchitecture> VCLIACTArray<VCLArchitecture>::
~VCLIACTArray()
{
  for(auto* propagator : propagator_) {
    if(propagator->adopt_propagator)delete propagator->propagator;
    if(propagator->adopt_pe_processor)delete propagator->pe_processor;
    delete propagator;
  }
  for(auto* detector : detector_) {
    delete detector;
  }
  for(auto* bandwidth_manager : bandwidth_manager_) {
    delete bandwidth_manager;
  }
}

template<typename VCLArchitecture> unsigned VCLIACTArray<VCLArchitecture>::
add_propagator_set(double scattering_radius, const std::string& name)
{
  auto* propagator_set = new PropagatorSet;
  propagator_set->ipropagator_set = propagator_set_.size();
  propagator_set->scattering_radius = scattering_radius;
  propagator_set->scattered_distance = 0.0;
  propagator_set->scattered_offset = Eigen::Vector3d::Zero();
  if(name == "") {
    propagator_set->name = "Propagator set "+std::to_string(propagator_set->ipropagator_set);
  } else {
    propagator_set->name = name;
  }
  propagator_set_.emplace_back(propagator_set);
  return propagator_set->ipropagator_set;
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
add_propagator(FocalPlaneRayPropagator* propagator, PEProcessor* pe_processor,
  VCLBandwidthManager<VCLArchitecture>* bandwidth_manager,
  SplinePEAmplitudeGenerator* pe_generator, const std::string& propagator_name,
  bool adopt_propagator, bool adopt_pe_processor, bool adopt_pe_generator)
{
  using calin::math::special::SQR;

  auto sphere = propagator->detector_spheres();

  auto* propagator_info = new PropagatorInfo;
  propagator_info->propagator          = propagator;
  propagator_info->ipropagator         = propagator_.size();
  propagator_info->pe_processor        = pe_processor;
  propagator_info->detector0           = detector_.size();
  propagator_info->ndetector           = sphere.size();
  propagator_info->bandwidth_manager   = bandwidth_manager;
  propagator_info->pe_generator        = pe_generator;
  propagator_info->adopt_propagator    = adopt_propagator;
  propagator_info->adopt_pe_processor  = adopt_pe_processor;
  propagator_info->adopt_pe_generator  = adopt_pe_generator;
  std::string name = propagator_name;
  if(name.empty()) {
    name = "propagator "+std::to_string(propagator_info->ipropagator);
  }
  propagator_info->name                = name;
  propagator_.emplace_back(propagator_info);

  for(unsigned isphere=0; isphere<sphere.size(); ++isphere) {
    auto detector_info = new DetectorInfo;
    detector_info->propagator_info        = propagator_info;
    detector_info->sphere                 = sphere[isphere];
    detector_info->unscattered_sphere     = sphere[isphere];
    detector_info->squared_radius         = SQR(sphere[isphere].radius);
    detector_info->squared_safety_radius  = SQR(sphere[isphere].radius + safety_radius_);
    detector_info->propagator             = propagator;
    detector_info->propagator_iscope      = isphere;
    detector_info->global_iscope          = detector_.size();
    detector_info->pe_processor           = pe_processor;
    detector_info->bandwidth_manager      = bandwidth_manager;
    detector_info->pe_generator           = pe_generator;
    detector_info->nrays_to_refract       = 0;
    detector_info->nrays_to_propagate     = 0;
    detector_info->nrays_refracted        = 0;
    detector_info->nrays_propagated       = 0;
    detector_.emplace_back(detector_info);
    propagator_info->detector_infos.emplace_back(detector_info);
  }

  bandwidth_manager_.emplace_back(bandwidth_manager);

  if(propagator_set_.empty()) {
    add_propagator_set(config_.scattering_radius(), "default_propagator_set");
  }
  PropagatorSet* propagator_set = propagator_set_.back();
  propagator_set->propagators.emplace_back(propagator_info);
  propagator_info->propagator_set = propagator_set;

  update_detector_efficiencies();
}

template<typename VCLArchitecture>
calin::simulation::vcl_ray_propagator::DaviesCottonVCLFocalPlaneRayPropagator<VCLArchitecture>*
VCLIACTArray<VCLArchitecture>::add_davies_cotton_propagator(
  calin::simulation::vs_optics::VSOArray* array, PEProcessor* pe_processor,
  const DetectionEfficiency& detector_efficiency, const AngularEfficiency& fp_angular_efficiency,
  SplinePEAmplitudeGenerator* pe_generator, const std::string& propagator_name,
  bool adopt_array, bool adopt_pe_processor, bool adopt_pe_generator)
{
  auto* propagator = new calin::simulation::vcl_ray_propagator::DaviesCottonVCLFocalPlaneRayPropagator<VCLArchitecture>(
    array, this->rng_, ref_index_, adopt_array, /* adopt_rng= */ false);

  auto* bandwidth_manager = new VCLDCBandwidthManager<VCLArchitecture>(
    &atm_abs_, detector_efficiency, fp_angular_efficiency, zobs_,
    config_.detector_energy_lo(), config_.detector_energy_hi(),
    config_.detector_energy_bin_width(), propagator_name);

  add_propagator(propagator, pe_processor, bandwidth_manager, pe_generator, propagator_name,
    /* adopt_propagator = */ true, adopt_pe_processor, adopt_pe_generator);

  return propagator;
}

template<typename VCLArchitecture>
calin::simulation::vcl_ray_propagator::DaviesCottonVCLFocalPlaneRayPropagator<VCLArchitecture>*
VCLIACTArray<VCLArchitecture>::add_davies_cotton_propagator(
  calin::simulation::vs_optics::VSOArray* array, PEProcessor* pe_processor,
  const DetectionEfficiency& detector_efficiency, const AngularEfficiency& fp_angular_efficiency,
  const std::string& propagator_name, SplinePEAmplitudeGenerator* pe_generator,
  bool adopt_array, bool adopt_pe_processor, bool adopt_pe_generator)
{
  return add_davies_cotton_propagator(array, pe_processor, detector_efficiency,
    fp_angular_efficiency, pe_generator, propagator_name, adopt_array, adopt_pe_processor,
    adopt_pe_generator);
}

template<typename VCLArchitecture>
calin::simulation::vcl_ray_propagator::DaviesCottonVCLFocalPlaneRayPropagator<VCLArchitecture>*
VCLIACTArray<VCLArchitecture>::add_davies_cotton_propagator(
  const ix::simulation::vs_optics::IsotropicDCArrayParameters& param, PEProcessor* pe_processor,
  const DetectionEfficiency& detector_efficiency, const AngularEfficiency& fp_angular_efficiency,
  SplinePEAmplitudeGenerator* pe_generator, const std::string& propagator_name,
  bool adopt_pe_processor, bool adopt_pe_generator)
{
  auto* array = new calin::simulation::vs_optics::VSOArray;
  calin::math::rng::VCLToScalarRNGCore scalar_core(this->rng_->core());
  calin::math::rng::RNG scalar_rng(&scalar_core);
  array->generateFromArrayParameters(param, scalar_rng);
  return add_davies_cotton_propagator(array, pe_processor, detector_efficiency,
    fp_angular_efficiency, pe_generator, propagator_name,
    /* adopt_array= */ true, adopt_pe_processor, adopt_pe_generator);
}

template<typename VCLArchitecture>
calin::simulation::vcl_ray_propagator::DaviesCottonVCLFocalPlaneRayPropagator<VCLArchitecture>*
VCLIACTArray<VCLArchitecture>::add_davies_cotton_propagator(
  const ix::simulation::vs_optics::IsotropicDCArrayParameters& param, PEProcessor* pe_processor,
  const DetectionEfficiency& detector_efficiency, const AngularEfficiency& fp_angular_efficiency,
  const std::string& propagator_name, SplinePEAmplitudeGenerator* pe_generator,
  bool adopt_pe_processor, bool adopt_pe_generator)
{
  return add_davies_cotton_propagator(param, pe_processor, detector_efficiency,
    fp_angular_efficiency, pe_generator, propagator_name, adopt_pe_processor,
    adopt_pe_generator);
}

template<typename VCLArchitecture>
calin::simulation::vcl_ray_propagator::PerfectOpticsVCLFocalPlaneRayPropagator<VCLArchitecture>*
VCLIACTArray<VCLArchitecture>::add_perfect_optics_propagator(
  const Eigen::VectorXd& x, const Eigen::VectorXd& y, const Eigen::VectorXd& z, 
  double radius, double focal_length, double field_of_view_radius, PEProcessor* pe_processor,
  const DetectionEfficiency& detector_efficiency,
  const std::string& propagator_name, SplinePEAmplitudeGenerator* pe_generator,
  bool adopt_pe_processor, bool adopt_pe_generator)
{
  if(x.size() != y.size() or x.size() != z.size()) {
    throw std::runtime_error("Number of telescope x, y and z coordinates must be equal");
  }
  auto* propagator = new calin::simulation::vcl_ray_propagator::PerfectOpticsVCLFocalPlaneRayPropagator<VCLArchitecture>(
    ref_index_);
  for(unsigned iscope=0; iscope<x.size(); ++iscope) {
    Eigen::Vector3d r0;
    r0 << x[iscope], y[iscope], z[iscope];
    propagator->add_telescope(r0, radius, config_.observation_level(), focal_length, field_of_view_radius);
  }

  auto* bandwidth_manager = new VCLSimpleBandwidthManager<VCLArchitecture>(
    &atm_abs_, detector_efficiency, zobs_,
    config_.detector_energy_lo(), config_.detector_energy_hi(),
    config_.detector_energy_bin_width(), propagator_name);

  add_propagator(propagator, pe_processor, bandwidth_manager, pe_generator, propagator_name,
    /* adopt_propagator = */ true, adopt_pe_processor, adopt_pe_generator);

  return propagator;
}

template<typename VCLArchitecture>
calin::simulation::vcl_ray_propagator::AllSkyVCLFocalPlaneRayPropagator<VCLArchitecture>*
VCLIACTArray<VCLArchitecture>::add_all_sky_propagator(
  Eigen::VectorXd& r0, double radius, double field_of_view_radius, PEProcessor* pe_processor,
  const DetectionEfficiency& detector_efficiency,
  const std::string& propagator_name, SplinePEAmplitudeGenerator* pe_generator,
  bool adopt_pe_processor, bool adopt_pe_generator)
{
  auto* propagator = new calin::simulation::vcl_ray_propagator::AllSkyVCLFocalPlaneRayPropagator<VCLArchitecture>(
    config_.observation_level(), r0, radius, field_of_view_radius, ref_index_);

  auto* bandwidth_manager = new VCLSimpleBandwidthManager<VCLArchitecture>(
    &atm_abs_, detector_efficiency, zobs_,
    config_.detector_energy_lo(), config_.detector_energy_hi(),
    config_.detector_energy_bin_width(), propagator_name);

  add_propagator(propagator, pe_processor, bandwidth_manager, pe_generator, propagator_name,
    /* adopt_propagator = */ true, adopt_pe_processor, adopt_pe_generator);

  return propagator;
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
point_telescope_az_el_phi_deg(unsigned iscope,
  double az_deg, double el_deg, double phi_deg)
{
  if(iscope >= detector_.size()) {
    throw std::out_of_range("Telescope ID out of range");
  }

  DetectorInfo* detector(detector_[iscope]);
  PropagatorInfo* ipropagator(detector->propagator_info);
  unsigned propagator_isphere = iscope-ipropagator->detector0;

  ipropagator->propagator->point_telescope_az_el_phi_deg(
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
  for(auto* propagator : propagator_) {
    for(unsigned propagator_isphere=0; propagator_isphere<propagator->ndetector;
        ++propagator_isphere) {
      unsigned isphere = propagator->detector0 + propagator_isphere;
      propagator->propagator->point_telescope_az_el_phi_deg(
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
  for(auto* propagator : propagator_) {
    auto spheres = propagator->propagator->detector_spheres();
    if(spheres.size() != propagator->ndetector) {
      // this should never happen
      throw std::runtime_error("Number of detectors proposed by propagator must remain constant over events.");
    }
    for(unsigned propagator_isphere=0; propagator_isphere<propagator->ndetector;
        ++propagator_isphere) {
      const auto& sphere(spheres[propagator_isphere]);
      auto* detector(propagator->detector_infos[propagator_isphere]);
      if(sphere.iobs != config_.observation_level()) {
        throw std::runtime_error("Detector observation level does not match configured value.");
      }
      detector->sphere = sphere;
      detector->unscattered_sphere = sphere;
      double zn = std::atan2(std::sqrt(SQR(sphere.obs_dir.x())+SQR(sphere.obs_dir.y())), sphere.obs_dir.z());
      znmin = std::min(znmin, std::max(zn - sphere.field_of_view_radius, 0.0));
      znmax = std::max(znmax, std::min(zn + sphere.field_of_view_radius, M_PI_2));
    }
    wmax_ = std::cos(znmin);
    wmin_ = std::cos(znmax);
    safety_radius_ = this->atm_->refraction_safety_radius(znmax, config_.observation_level());
    for(unsigned propagator_isphere=0; propagator_isphere<propagator->ndetector;
        ++propagator_isphere) {
      auto* detector(propagator->detector_infos[propagator_isphere]);
      const auto& sphere(detector->sphere);
      detector->squared_radius         = SQR(sphere.radius);
      detector->squared_safety_radius  = SQR(sphere.radius + safety_radius_);
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
  for(const auto* ibandwidth_manager : bandwidth_manager_) {
    bandwidth = std::max(bandwidth, ibandwidth_manager->bandwidth());
  }
  return bandwidth;
}

template<typename VCLArchitecture> calin::math::spline_interpolation::CubicSpline*
VCLIACTArray<VCLArchitecture>::new_height_dependent_pe_bandwidth_spline() const
{
  if(bandwidth_manager_.empty()) {
    return nullptr;
  }
  std::vector<double> heights = bandwidth_manager_.front()->
    detector_bandwidth_spline()->xknot_as_stdvec();
  std::vector<double> bandwidths(heights.size(), 0.0);
  for(const auto* ibandwidth_manager : bandwidth_manager_) {
    std::vector<double> detector_bandwidths =
      ibandwidth_manager->bandwidth_vs_height(heights, wmax_);
    std::transform(bandwidths.begin(), bandwidths.end(),
      detector_bandwidths.begin(), bandwidths.begin(),
      [](double a, double b) { return std::max(a,b); });
  }
  return new calin::math::spline_interpolation::CubicSpline(heights, bandwidths);
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{
  calin::math::rng::VCLToScalarRNGCore scalar_core(this->rng_->core());
  calin::math::rng::RNG scalar_rng(&scalar_core);
  Eigen::Vector3d e1(1.0, 0.0, 0.0);
  Eigen::Vector3d e2(0.0, 1.0, 0.0);
  calin::math::geometry::rotate_in_place_z_to_u_Rzy(e1, event.u0);
  calin::math::geometry::rotate_in_place_z_to_u_Rzy(e2, event.u0);
  for(auto* propagator : propagator_) {
    propagator->pe_processor->start_processing();
  }
  for(auto* detector : detector_) {
    detector->nrays_to_refract = 0;
    detector->nrays_to_propagate = 0;
    detector->nrays_refracted = 0;
    detector->nrays_propagated = 0;
    detector->sphere = detector->unscattered_sphere;
  }
  for(auto* propagator_set : propagator_set_) {
    if(propagator_set->scattering_radius > 0.0) {
      double b = propagator_set->scattering_radius*std::sqrt(scalar_rng.uniform());
      propagator_set->scattered_distance = b;
      double theta = scalar_rng.uniform()*M_PI*2.0;
      double bx = b*std::cos(theta);
      double by = b*std::sin(theta);
      Eigen::Vector3d bvec = bx*e1 + by*e2;
      bvec -= (bvec.z()/event.u0.z())*event.u0;
      propagator_set->scattered_offset = bvec;
      for(auto* propagator_info : propagator_set->propagators) {
        for(auto* detector_info : propagator_info->detector_infos) {
          detector_info->sphere.r0 += bvec;
        }
      }
    } else {
      propagator_set->scattered_offset = Eigen::Vector3d::Zero();
    }
  }
  make_detector_grid();

  return VCLIACTTrackVisitor<VCLArchitecture>::visit_event(event, kill_event);
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::leave_event()
{
  VCLIACTTrackVisitor<VCLArchitecture>::leave_event();

  for(auto* detector : detector_) {
    if(detector->nrays_to_refract) {
      do_refract_rays_for_detector(detector);
    }
    if(detector->nrays_to_propagate) {
      do_propagate_rays_for_detector(detector);
    }
  }
  for(auto* propagator : propagator_) {
    propagator->pe_processor->finish_processing();
  }
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

  int64_vt grid_hexid;
  if(grid_ncells_ == 1) {
    // If there is only one cell, we can skip the grid calculation
    grid_hexid = 0;
  } else {
    calin::math::ray::VCLRay<double_real> ray_copy(ray);
    ray_copy.propagate_to_z_plane_with_mask(ray_mask, -zobs_);
    double_vt grid_x = ray_copy.x() * grid_sep_inv_;
    double_vt grid_y = ray_copy.y() * grid_sep_inv_;
    int64_vt u;
    int64_vt v;
    calin::math::hex_array::VCL<VCLArchitecture>::xy_to_uv(grid_x, grid_y, u, v);
    grid_hexid = calin::math::hex_array::VCL<VCLArchitecture>::uv_to_hexid_ccw(u, v);
  }

  for(unsigned icell_detector=0;icell_detector<grid_ndetector_per_cell_; ++icell_detector) {
    int64_vt igrid_array = 
      vcl::min(grid_hexid,grid_ncells_)*grid_ndetector_per_cell_ + icell_detector;
    double_vt detector_x = vcl::lookup<0x40000000>(igrid_array, grid_x);
    double_vt detector_y = vcl::lookup<0x40000000>(igrid_array, grid_y);
    double_vt detector_z = vcl::lookup<0x40000000>(igrid_array, grid_z);
    double_vt detector_squared_safety_radius = vcl::lookup<0x40000000>(igrid_array, grid_ssr);
    int64_vt idetector = vcl::lookup<0x40000000>(igrid_array, grid_idetector);
    Vector3d_vt detector_pos { detector_x, detector_y, detector_z };

    auto intersecting_rays = ray_mask 
      & (ray.squared_distance_at_closest_approach(detector_pos) < detector_squared_safety_radius)
      & double_bvt(idetector != -1);
    unsigned intersecting_rays_bitmask = vcl::to_bits(intersecting_rays);

    int64_at idetector_array;
    idetector.store(idetector_array);

    for(unsigned iray=0; iray<VCLArchitecture::num_double; ++iray) {
      if(intersecting_rays_bitmask & 1) {
        auto* detector = detector_[idetector_array[iray]];
        switch(config_.refraction_mode()) {
        case calin::ix::simulation::vcl_iact::REFRACT_NO_RAYS:
        case calin::ix::simulation::vcl_iact::REFRACT_ALL_RAYS:
          detector->rays_to_propagate.insert_one_ray(detector->nrays_to_propagate, ray_array.extract_one_ray(iray));
          detector->bandwidths_to_propagate[detector->nrays_to_propagate] = bandwidth_array[iray];
          detector->ray_weights_to_propagate[detector->nrays_to_propagate] = ray_weight_array[iray];
          ++detector->nrays_to_propagate;
          if(detector->nrays_to_propagate == VCLArchitecture::num_double) {
            do_propagate_rays_for_detector(detector);
          }
          break;
        case calin::ix::simulation::vcl_iact::REFRACT_ONLY_CLOSE_RAYS:
        default:
          detector->rays_to_refract.insert_one_ray(detector->nrays_to_refract, ray_array.extract_one_ray(iray));
          detector->bandwidths_to_refract[detector->nrays_to_refract] = bandwidth_array[iray];
          detector->ray_weights_to_refract[detector->nrays_to_refract] = ray_weight_array[iray];
          ++detector->nrays_to_refract;
          if(detector->nrays_to_refract == VCLArchitecture::num_double) {
            do_refract_rays_for_detector(detector);
          }
          break;
        }
      }
      intersecting_rays_bitmask >>= 1;
    }
  }
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
do_refract_rays_for_detector(DetectorInfo* detector)
{
  Ray ray;
  detector->rays_to_refract.get_rays(ray);
  double_bvt ray_mask = VCLArchitecture::double_iota()<detector->nrays_to_refract;
  detector->nrays_refracted += detector->nrays_to_refract;
  detector->nrays_to_refract = 0;

  double_vt dz = ray.z()-zobs_;
  this->atm_->template vcl_propagate_ray_with_refraction_and_mask<VCLArchitecture>(ray, ray_mask, config_.observation_level());
  // Note propagation backwards by distance since uz() is negative
  ray.propagate_dist_with_mask(ray_mask, dz/ray.uz(), ref_index_);

  detector->rays_to_refract.set_rays(ray);

  auto intersecting_rays = ray_mask & (ray.squared_distance_at_closest_approach(detector->sphere.r0.template cast<double_vt>()) < detector->squared_radius);
  unsigned intersecting_rays_bitmask = vcl::to_bits(intersecting_rays);
  if(intersecting_rays_bitmask) {
    for(unsigned iray=0; iray<VCLArchitecture::num_double; ++iray) {
      if(intersecting_rays_bitmask & 1) {
        detector->rays_to_propagate.insert_one_ray(detector->nrays_to_propagate, detector->rays_to_refract.extract_one_ray(iray));
        detector->bandwidths_to_propagate[detector->nrays_to_propagate] = detector->bandwidths_to_refract[iray];
        detector->ray_weights_to_propagate[detector->nrays_to_propagate] = detector->ray_weights_to_refract[iray];
        ++detector->nrays_to_propagate;
        if(detector->nrays_to_propagate == VCLArchitecture::num_double) {
          do_propagate_rays_for_detector(detector);
        }
      }
      intersecting_rays_bitmask >>= 1;
    }
  }
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::
do_propagate_rays_for_detector(DetectorInfo* detector)
{
  Ray ray;
  detector->rays_to_propagate.get_rays(ray);
  double_bvt ray_mask = VCLArchitecture::double_iota()<detector->nrays_to_propagate;
  detector->nrays_propagated += detector->nrays_to_propagate;
  detector->nrays_to_propagate = 0;

  double_vt emission_z = ray.z();
  double_vt emission_uz = ray.uz();

  if(detector->propagator_info->propagator_set->scattering_radius > 0.0) {
    // Translate the ray to unscattered frame since this is how the propagator works
    ray.translate_origin(detector->propagator_info->propagator_set->scattered_offset.template cast<double_vt>());
  }

  FocalPlaneParameters fp_parameters;
  ray_mask = detector->propagator->propagate_rays_to_focal_plane(
    detector->propagator_iscope, ray, ray_mask, fp_parameters);

  if(not this->do_color_photons_) {
    double_vt emission_bandwidth;
    emission_bandwidth.load(detector->bandwidths_to_propagate);
    double_vt ground_bandwidth =
      detector->bandwidth_manager->bandwidth_for_pe(emission_z, emission_uz,
        fp_parameters.fplane_ux, fp_parameters.fplane_uy, fp_parameters.fplane_uz);
    double_vt bw_scaled_uniform_rand = emission_bandwidth * this->rng_->uniform_double();
    double_vt bw_scaled_detection_prob = fp_parameters.detection_prob * ground_bandwidth;
    ray_mask &= bw_scaled_uniform_rand < bw_scaled_detection_prob;
  }

  unsigned fp_rays_bitmask = vcl::to_bits(ray_mask);
  if(fp_rays_bitmask) {
    double_at fplane_x;
    double_at fplane_y;
    double_at fplane_ux;
    double_at fplane_uy;
    double_at fplane_t;
    int64_at pixel_id;

    if(this->do_color_photons_) {

    }

    fp_parameters.fplane_x.store(fplane_x);
    fp_parameters.fplane_z.store(fplane_y);
    fp_parameters.fplane_ux.store(fplane_ux);
    fp_parameters.fplane_uz.store(fplane_uy);
    fp_parameters.fplane_t.store(fplane_t);
    fp_parameters.pixel_id.store(pixel_id);

    if(detector->pe_generator != nullptr) {
      double_vt weight;
      weight.load(detector->ray_weights_to_propagate);
      weight *= detector->pe_generator->
        template vcl_generate_amplitude<VCLArchitecture>(*this->rng_);
      weight.store(detector->ray_weights_to_propagate);
    }

    for(unsigned iray=0; iray<VCLArchitecture::num_double; ++iray) {
      if(fp_rays_bitmask & 1) {
        detector->pe_processor->process_focal_plane_hit(detector->propagator_iscope,
          pixel_id[iray], fplane_x[iray], fplane_y[iray], fplane_ux[iray], fplane_uy[iray],
          fplane_t[iray], detector->ray_weights_to_propagate[iray]);
      }
      fp_rays_bitmask >>= 1;
    }
  }
}

namespace {
  inline double projected_radius(double r, double zcenter, double zobs, double sin_zn_max, double tan_zn_max)
  {
    double dz = std::abs(zcenter - zobs) + r/sin_zn_max;
    return dz*tan_zn_max;
  }
}

template<typename VCLArchitecture> void VCLIACTArray<VCLArchitecture>::make_detector_grid()
{
  // Create a hexagonal grid of detectors based on the detector positions
  // and the observation zenith angle. The grid is used to efficiently
  // find detectors that are close to a given ray at the observation level.

  if(detector_.size() <= std::max(1U, config_.grid_theshold())) 
  {
    grid_sep_ = std::numeric_limits<double>::infinity();
    grid_sep_inv_ = 0.0;

    grid_ncells_ = 1;
    grid_ndetector_per_cell_ = detector_.size();
    unsigned array_size = (grid_ncells_ + 1) * grid_ndetector_per_cell_;
    if(array_size > grid_array_size_) {
      grid_array_size_ = 2*array_size;
      calin::util::memory::safe_aligned_recalloc(grid_x, grid_array_size_);
      calin::util::memory::safe_aligned_recalloc(grid_y, grid_array_size_);
      calin::util::memory::safe_aligned_recalloc(grid_z, grid_array_size_);
      calin::util::memory::safe_aligned_recalloc(grid_ssr, grid_array_size_);
      calin::util::memory::safe_aligned_recalloc(grid_idetector, grid_array_size_);
    }

    std::fill(grid_x, grid_x+grid_array_size_, 0.0);
    std::fill(grid_y, grid_y+grid_array_size_, 0.0);
    std::fill(grid_z, grid_z+grid_array_size_, 0.0);
    std::fill(grid_ssr, grid_ssr+grid_array_size_, 0.0);
    std::fill(grid_idetector, grid_idetector+grid_array_size_, -1);

    for(int idetector=0; idetector<int(detector_.size()); ++idetector) {
      int ii = idetector;
      grid_x[ii] = detector_[idetector]->sphere.r0.x();
      grid_y[ii] = detector_[idetector]->sphere.r0.y();
      grid_z[ii] = detector_[idetector]->sphere.r0.z();
      grid_ssr[ii] = detector_[idetector]->squared_safety_radius;
      grid_idetector[ii] = detector_[idetector]->global_iscope;
    }

    return;    
  }

  using namespace calin::util::log;

  const double cos60 = 0.5;
  const double sin60 = 0.5*CALIN_HEX_ARRAY_SQRT3;

  double zn = std::acos(wmin_);
  double sin_zn_max = std::sin(zn);
  double tan_zn_max = std::tan(zn);

  if(config_.grid_separation() > 0.0) {
    grid_sep_ = config_.grid_separation();
  } else {
    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymin = std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();
    double rmax = 0.0;
    for(auto* detector : detector_) {
      double r = projected_radius(detector->sphere.radius + safety_radius_, 
        detector->sphere.r0.z(), zobs_, sin_zn_max, tan_zn_max);
      xmin = std::min(xmin, detector->sphere.r0.x()-r);
      xmax = std::max(xmax, detector->sphere.r0.x()+r);
      ymin = std::min(ymin, detector->sphere.r0.y()-r);
      ymax = std::max(ymax, detector->sphere.r0.y()+r);
      rmax = std::max(rmax, r);
    }
    double area = (xmax-xmin)*(ymax-ymin);
    grid_sep_ = std::max(4*rmax, std::sqrt(area/std::max(config_.grid_area_divisor(),1.0)));
  }
  grid_sep_inv_ = 1.0/grid_sep_;

  std::map<unsigned, std::vector<DetectorInfo*>> detector_grid;
  for(auto* detector : detector_) {
    double x = detector->sphere.r0.x()*grid_sep_inv_;
    double y = detector->sphere.r0.y()*grid_sep_inv_;
    double r = projected_radius(detector->sphere.radius + safety_radius_, 
      detector->sphere.r0.z(), zobs_, sin_zn_max, tan_zn_max)*grid_sep_inv_;
    double dx = x;
    double dy = y;
    int u;
    int v;
    calin::math::hex_array::xy_to_uv_with_remainder(dx,dy,u,v);
    unsigned hexid = calin::math::hex_array::uv_to_hexid(u,v);
    detector_grid[hexid].emplace_back(detector);

    double dx_cos60 = dx * cos60;
    double dy_sin60 = dy * sin60;

    double dx_pos60 = std::abs(dx_cos60 - dy_sin60);
    double dx_neg60 = std::abs(dx_cos60 + dy_sin60);

    if(dx+r > 0.5) {
      // Add the detector to the next hexagon in the x direction
      unsigned hexid_next = calin::math::hex_array::uv_to_hexid(u+1,v);
      detector_grid[hexid_next].emplace_back(detector);
    }
    if(dx_neg60+r > 0.5) {
      // Add the detector to the next hexagon in the x-60 direction
      unsigned hexid_next = calin::math::hex_array::uv_to_hexid(u,v+1);
      detector_grid[hexid_next].emplace_back(detector);
    }
    if(dx_pos60+r > 0.5) {
      // Add the detector to the next hexagon in the x+60 direction
      unsigned hexid_next = calin::math::hex_array::uv_to_hexid(u+1,v-1);
      detector_grid[hexid_next].emplace_back(detector);
    }
    if(dx-r < -0.5) {
      // Add the detector to the next hexagon in the -x direction
      unsigned hexid_next = calin::math::hex_array::uv_to_hexid(u-1,v);
      detector_grid[hexid_next].emplace_back(detector);
    }  
    if(dx_neg60-r < -0.5) {
      // Add the detector to the next hexagon in the -x-60 direction
      unsigned hexid_next = calin::math::hex_array::uv_to_hexid(u,v-1);
      detector_grid[hexid_next].emplace_back(detector);
    }
    if(dx_pos60-r < -0.5) {
      // Add the detector to the next hexagon in the -x+60 direction
      unsigned hexid_next = calin::math::hex_array::uv_to_hexid(u-1,v+1);
      detector_grid[hexid_next].emplace_back(detector);
    }
  }
  unsigned max_detectors_per_cell = 0;
  unsigned hexid_max = 0;
  for(auto& [hexid, detectors] : detector_grid) {
    hexid_max = std::max(hexid_max, hexid);
    max_detectors_per_cell = std::max(max_detectors_per_cell, unsigned(detectors.size()));
  }

  grid_ncells_ = hexid_max+1;
  grid_ndetector_per_cell_ = max_detectors_per_cell;
  unsigned array_size = (grid_ncells_ + 1) * grid_ndetector_per_cell_;
  if(array_size > grid_array_size_) {
    grid_array_size_ = 2*array_size;
    calin::util::memory::safe_aligned_recalloc(grid_x, grid_array_size_);
    calin::util::memory::safe_aligned_recalloc(grid_y, grid_array_size_);
    calin::util::memory::safe_aligned_recalloc(grid_z, grid_array_size_);
    calin::util::memory::safe_aligned_recalloc(grid_ssr, grid_array_size_);
    calin::util::memory::safe_aligned_recalloc(grid_idetector, grid_array_size_);
  }

  std::fill(grid_x, grid_x+grid_array_size_, 0.0);
  std::fill(grid_y, grid_y+grid_array_size_, 0.0);
  std::fill(grid_z, grid_z+grid_array_size_, 0.0);
  std::fill(grid_ssr, grid_ssr+grid_array_size_, 0.0);
  std::fill(grid_idetector, grid_idetector+grid_array_size_, -1);

  for(auto& [hexid, detectors] : detector_grid) {
    for(int idetector=0; idetector<int(detectors.size()); ++idetector) {
      int ii = hexid*grid_ndetector_per_cell_ + idetector;
      grid_x[ii] = detectors[idetector]->sphere.r0.x();
      grid_y[ii] = detectors[idetector]->sphere.r0.y();
      grid_z[ii] = detectors[idetector]->sphere.r0.z();
      grid_ssr[ii] = detectors[idetector]->squared_safety_radius;
      grid_idetector[ii] = detectors[idetector]->global_iscope;
    }
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
  config.set_grid_theshold(4);
  config.set_grid_area_divisor(256.0);
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

template<typename VCLArchitecture> std::string VCLIACTArray<VCLArchitecture>::banner() const
{
  constexpr double EV_NM = 1239.84193009239; // gunits: c/(ev/h) -> nm
  using calin::util::string::double_to_string_with_commas;
  using calin::math::special::SQR;
  std::ostringstream stream;
  double prop_delay_ct =
    this->atm_->propagation_ct_correction(this->atm_->top_of_atmosphere()) -
    this->atm_->propagation_ct_correction(zobs_);
  double thetac_zobs = std::acos(1/(this->atm_->n_minus_one(zobs_)+1));
  double thetac_5 = std::acos(1/(this->atm_->n_minus_one(5e5)+1));
  double thetac_10 = std::acos(1/(this->atm_->n_minus_one(10e5)+1));
  double thetac_15 = std::acos(1/(this->atm_->n_minus_one(15e5)+1));
  double yield_zobs = calin::simulation::air_cherenkov_tracker::YIELD_CONST*SQR(std::sin(thetac_zobs))*100;
  double yield_5 = calin::simulation::air_cherenkov_tracker::YIELD_CONST*SQR(std::sin(thetac_5))*100;
  double yield_10 = calin::simulation::air_cherenkov_tracker::YIELD_CONST*SQR(std::sin(thetac_10))*100;
  double yield_15 = calin::simulation::air_cherenkov_tracker::YIELD_CONST*SQR(std::sin(thetac_15))*100;
  stream
    << "Class : " << calin::util::vcl::templated_class_name<VCLArchitecture>("VCLIACTArray") << '\n'
    << "Number of focal-plane propagators : " << propagator_.size() << ", with "
    << detector_.size() << " detectors";
  if(propagator_set_.size()) {
    stream << ", in " << propagator_set_.size() << " sets.\n";
  } else {
    stream << ".\n";
  }
  
  std::vector<std::pair<std::string, unsigned>> message_counts;
  for(const auto* ipropagator : propagator_) {
    auto banner = ipropagator->propagator->banner("- "+ipropagator->name+": ", "  ");
    bool banner_found = false;
    for(auto& message : message_counts) {
      if(message.first == banner) {
        message.second++;
        banner_found = true;
        break;
      }
    }
    if(not banner_found) {
      message_counts.emplace_back(banner, 1);
    }
  }
  for(const auto& message : message_counts) {
    if(message.second == 1) {
      stream << message.first << '\n';
    } else if(message.second > 1) {
      stream << message.first << " (x" << message.second << ")\n";
    }
  }
  stream
    << "Detector zenith range : " << double_to_string_with_commas(std::acos(wmax_)/M_PI*180.0,1)
    << " to " << double_to_string_with_commas(std::acos(wmin_)/M_PI*180.0,1) << " degrees.\n"
    << "Observation level : " << double_to_string_with_commas(zobs_/1e5,3) << " km, "
    << "thickness " << double_to_string_with_commas(this->atm_->thickness(zobs_),1) << " g/cm^2\n"
    << "- Cherenkov angle at " << double_to_string_with_commas(zobs_/1e5,3) << ", 5, 10, 15 km : "
    << double_to_string_with_commas(thetac_zobs/M_PI*180,3) << ", "
    << double_to_string_with_commas(thetac_5/M_PI*180,3) << ", "
    << double_to_string_with_commas(thetac_10/M_PI*180,3) << ", "
    << double_to_string_with_commas(thetac_15/M_PI*180,3) << " deg\n"
    << "- Cherenkov yield at " << double_to_string_with_commas(zobs_/1e5,3) << ", 5, 10, 15 km : "
    << double_to_string_with_commas(yield_zobs,2) << ", "
    << double_to_string_with_commas(yield_5,2) << ", "
    << double_to_string_with_commas(yield_10,2) << ", "
    << double_to_string_with_commas(yield_15,2) <<  " ph/m/eV\n"
    << "- Vertical propagation delay from "
    << double_to_string_with_commas(this->atm_->top_of_atmosphere()/1e5,0)
    << " to " << double_to_string_with_commas(zobs_/1e5,3) << " km : "
    << double_to_string_with_commas(prop_delay_ct*0.03335641,2) <<  "ns ("
    << double_to_string_with_commas(prop_delay_ct,1) << " cm)\n";
  if(this->variable_bandwidth_spline_) {
    double bw_zobs = this->variable_bandwidth_spline_->value(zobs_);
    double bw_5 = this->variable_bandwidth_spline_->value(5e5);
    double bw_10 = this->variable_bandwidth_spline_->value(10e5);
    double bw_15 = this->variable_bandwidth_spline_->value(15e5);
    double bw_toa = this->variable_bandwidth_spline_->value(this->atm_->top_of_atmosphere());
    stream << "Cherenkov ray mode : PEs WITH HEIGHT-DEPENDENT BANDWIDTH\n"
      << "- " << double_to_string_with_commas(zobs_/1e5,3) << ", 5, 10, 15, "
      << double_to_string_with_commas(this->atm_->top_of_atmosphere()/1e5, 0) << " km : "
      << double_to_string_with_commas(bw_zobs,3) << ", "
      << double_to_string_with_commas(bw_5,3) << ", "
      << double_to_string_with_commas(bw_10,3) << ", "
      << double_to_string_with_commas(bw_15,3) << ", "
      << double_to_string_with_commas(bw_toa,3) << " eV\n"
      << "- PE yield at " << double_to_string_with_commas(zobs_/1e5,3) << ", 5, 10, 15 km : "
      << double_to_string_with_commas(yield_zobs*bw_zobs,2) << ", "
      << double_to_string_with_commas(yield_5*bw_5,2) << ", "
      << double_to_string_with_commas(yield_10*bw_10,2) << ", "
      << double_to_string_with_commas(yield_15*bw_15,2) <<  " PE/m\n";
  } else if (this->do_color_photons_) {
    stream << "Cherenkov ray mode : PHOTONS WITH FIXED BANDWIDTH "
      << double_to_string_with_commas(this->fixed_bandwidth_,3) << " eV\n"
      << "- Energy range "
      << double_to_string_with_commas(this->min_cherenkov_energy_,3) << " - "
      << double_to_string_with_commas(this->min_cherenkov_energy_+this->fixed_bandwidth_,3) << " eV ("
      << double_to_string_with_commas(EV_NM/(this->min_cherenkov_energy_+this->fixed_bandwidth_),0) << " - "
      << double_to_string_with_commas(EV_NM/this->min_cherenkov_energy_,0) << " nm)\n";
  } else {
    stream << "Cherenkov ray mode : PEs WITH FIXED BANDWIDTH "
      << double_to_string_with_commas(this->fixed_bandwidth_,3) << " eV\n"
      << "- PE yield at " << double_to_string_with_commas(zobs_/1e5,3) << ", 5, 10, 15 km : "
      << double_to_string_with_commas(yield_zobs*this->fixed_bandwidth_,2) << ", "
      << double_to_string_with_commas(yield_5*this->fixed_bandwidth_,2) << ", "
      << double_to_string_with_commas(yield_10*this->fixed_bandwidth_,2) << ", "
      << double_to_string_with_commas(yield_15*this->fixed_bandwidth_,2) <<  " PE/cm\n";
  }
  if(config_.refraction_mode() == calin::ix::simulation::vcl_iact::REFRACT_NO_RAYS) {
    stream << "Refraction mode : DISABLED (NO REFRACTION APPLIED)\n";
  } else if(config_.refraction_mode() == calin::ix::simulation::vcl_iact::REFRACT_ALL_RAYS) {
    stream << "Refraction mode : REFRACT ALL RAYS\n";
  } else { // config_.refraction_mode() == calin::ix::simulation::vcl_iact::REFRACT_ONLY_CLOSE_RAYS
    stream << "Refraction mode : REFRACT ONLY CLOSE RAYS\n"
      << "- Refraction safety radius : " << double_to_string_with_commas(safety_radius_*0.01,2) << " m\n";
  }
  if(config_.refraction_mode() != calin::ix::simulation::vcl_iact::REFRACT_NO_RAYS) {
    stream
      << "- Displacement from 5, 10, 15 km at Zn="
      << double_to_string_with_commas(std::acos(wmin_)/M_PI*180,1) << " deg : "
      << double_to_string_with_commas(this->atm_->refraction_displacement(5e5, std::acos(wmin_), config_.observation_level())*0.01,2) << ", "
      << double_to_string_with_commas(this->atm_->refraction_displacement(10e5, std::acos(wmin_), config_.observation_level())*0.01,2) << ", "
      << double_to_string_with_commas(this->atm_->refraction_displacement(15e5, std::acos(wmin_), config_.observation_level())*0.01,2) << " m\n"
      << "- Displacement from 5, 10, 15 km at Zn="
      << double_to_string_with_commas(std::acos(wmax_)/M_PI*180,1) << " deg : "
      << double_to_string_with_commas(this->atm_->refraction_displacement(5e5, std::acos(wmax_), config_.observation_level())*0.01,2) << ", "
      << double_to_string_with_commas(this->atm_->refraction_displacement(10e5, std::acos(wmax_), config_.observation_level())*0.01,2) << ", "
      << double_to_string_with_commas(this->atm_->refraction_displacement(15e5, std::acos(wmax_), config_.observation_level())*0.01,2) << " m\n"
      << "- Bending from 5, 10, 15 km at Zn="
      << double_to_string_with_commas(std::acos(wmin_)/M_PI*180,1) << " deg : "
      << double_to_string_with_commas(this->atm_->refraction_bending(5e5, std::acos(wmin_), config_.observation_level())/M_PI*180*3600,1) << ", "
      << double_to_string_with_commas(this->atm_->refraction_bending(10e5, std::acos(wmin_), config_.observation_level())/M_PI*180*3600,1) << ", "
      << double_to_string_with_commas(this->atm_->refraction_bending(15e5, std::acos(wmin_), config_.observation_level())/M_PI*180*3600,1) << " arcsec\n"
      << "- Bending from 5, 10, 15 km at Zn="
      << double_to_string_with_commas(std::acos(wmax_)/M_PI*180,1) << " deg : "
      << double_to_string_with_commas(this->atm_->refraction_bending(5e5, std::acos(wmax_), config_.observation_level())/M_PI*180*3600,1) << ", "
      << double_to_string_with_commas(this->atm_->refraction_bending(10e5, std::acos(wmax_), config_.observation_level())/M_PI*180*3600,1) << ", "
      << double_to_string_with_commas(this->atm_->refraction_bending(15e5, std::acos(wmax_), config_.observation_level())/M_PI*180*3600,1) << " arcsec\n";
  }

  stream << "Detector efficiency bandwidths :\n";
  message_counts.clear();
  for(const auto* ibwm : bandwidth_manager_) {
    auto banner = ibwm->banner(wmin_, wmax_, "- ", "  ");
    bool banner_found = false;
    for(auto& message : message_counts) {
      if(message.first == banner) {
        message.second++;
        banner_found = true;
        break;
      }
    }
    if(not banner_found) {
      message_counts.emplace_back(banner, 1);
    }
  }
  for(const auto& message : message_counts) {
    if(message.second == 1) {
      stream << message.first << '\n';
    } else if(message.second > 1) {
      stream << message.first << " (x" << message.second << ")\n";
    }
  }

  message_counts.clear();
  for(const auto* ipropagator : propagator_) {
    if(ipropagator->pe_generator != nullptr) {
      auto banner = ipropagator->pe_generator->banner("- "+ipropagator->name+": ", "  ");
      bool banner_found = false;
      for(auto& message : message_counts) {
        if(message.first == banner) {
          message.second++;
          banner_found = true;
          break;
        }
      }
      if(not banner_found) {
        message_counts.emplace_back(banner, 1);
      }
    }
  }

  if(message_counts.empty()) {
    stream << "No photo-electron spectra configured.\n";
  } else {
    stream << "Single photo-electron spectra :\n";
    for(const auto& message : message_counts) {
      if(message.second == 1) {
        stream << message.first << '\n';
      } else if(message.second > 1) {
        stream << message.first << " (x" << message.second << ")\n";
      }
    }
  }
  return stream.str();
}

#endif // not defined SWIG

} } } // namespace calin::simulation::vcl_iact
