/*

   calin/simulation/vcl_detector_efficiency.hpp -- Stephen Fegan -- 2022-07-26

   Classes for implementing detector efficiecy curves with VCL vectors.

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

#include <sstream>
#include <cmath>
#include <stdexcept>

#include <util/vcl.hpp>
#include <util/string.hpp>
#include <simulation/detector_efficiency.hpp>

namespace calin { namespace simulation { namespace detector_efficiency {

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCLPEAmplitudeGenerator
{
public:
#ifndef SWIG
  using double_vt = typename VCLArchitecture::double_vt;
  using double_at = typename VCLArchitecture::double_at;
  using VCLRNG = calin::math::rng::VCLRNG<VCLArchitecture>;

  virtual double_vt vcl_generate_amplitude(VCLRNG& rng) const = 0;
#endif

  Eigen::VectorXd bulk_generate_amplitude(unsigned n,
      calin::math::rng::VCLRNG<VCLArchitecture>& rng) const
#ifndef SWIG
  {
    Eigen::VectorXd amplitudes(n);
    constexpr unsigned num_double = VCLArchitecture::num_double;
    unsigned i = 0;
    for(; i+num_double<=n; i+=num_double) {
      vcl_generate_amplitude(rng).store(amplitudes.data()+i);
    }
    if(i<n) {
      double_at last_amplitudes;
      vcl_generate_amplitude(rng).store(last_amplitudes);
      for(unsigned ilane=0; i<n; ++i, ++ilane) {
        amplitudes[i] = last_amplitudes[ilane];
      }
    }
    return amplitudes;
  }
#else
  ;
#endif

  virtual ~VCLPEAmplitudeGenerator() = default;
  virtual double mean_amplitude() = 0;
  virtual std::string banner(const std::string& indent0="", const std::string& indentN="") const = 0;
};

template<typename VCLArchitecture>
class alignas(VCLArchitecture::vec_bytes) VCLSplinePEAmplitudeGenerator:
  public VCLPEAmplitudeGenerator<VCLArchitecture>
{
public:
#ifndef SWIG
  using VCLRNG = calin::math::rng::VCLRNG<VCLArchitecture>;
  using double_vt = typename VCLArchitecture::double_vt;
#endif

  VCLSplinePEAmplitudeGenerator(SplinePEAmplitudeGenerator* pe_generator, bool adopt_pe_generator):
    pe_generator_(pe_generator), adopt_pe_generator_(adopt_pe_generator)
  {
    if(pe_generator_ == nullptr) {
      throw std::invalid_argument("Spline PE amplitude generator must not be null");
    }
  }

  ~VCLSplinePEAmplitudeGenerator() override {
    if(adopt_pe_generator_) delete pe_generator_;
  }

  double mean_amplitude() final {
    return pe_generator_->mean_amplitude();
  }

  std::string banner(const std::string& indent0="", const std::string& indentN="") const final {
    return pe_generator_->banner(indent0, indentN);
  }

#ifndef SWIG
  double_vt vcl_generate_amplitude(VCLRNG& rng) const final {
    return pe_generator_->template vcl_generate_amplitude<VCLArchitecture>(rng);
  }
#endif

private:
  SplinePEAmplitudeGenerator* pe_generator_ = nullptr;
  bool adopt_pe_generator_ = false;
};

template<typename VCLArchitecture>
class alignas(VCLArchitecture::vec_bytes) VCLSimpleSiPMPEAmplitudeGenerator:
  public VCLPEAmplitudeGenerator<VCLArchitecture>
{
public:
#ifndef SWIG
  using VCLRNG = calin::math::rng::VCLRNG<VCLArchitecture>;
  using double_vt = typename VCLArchitecture::double_vt;
#endif

  VCLSimpleSiPMPEAmplitudeGenerator(double crosstalk_mean, double spe_resolution):
    crosstalk_mean_(crosstalk_mean), exp_neg_crosstalk_mean_(std::exp(-crosstalk_mean)),
    spe_resolution_(spe_resolution)
  {
    if(not (crosstalk_mean_ >= 0.0 and crosstalk_mean_ < 1.0)) {
      throw std::domain_error("SiPM crosstalk mean must be in [0, 1)");
    }
    if(not (spe_resolution_ >= 0.0 and std::isfinite(spe_resolution_))) {
      throw std::domain_error("SiPM SPE resolution must be finite and nonnegative");
    }
  }

  double mean_amplitude() final {
    return 1.0 / (1.0 - crosstalk_mean_);
  }

  std::string banner(const std::string& indent0="", const std::string& indentN="") const final {
    using calin::util::string::double_to_string_with_commas;
    std::ostringstream stream;
    stream << indent0 << "Simple SiPM PE amplitude model"
      << '\n' << indentN << "Crosstalk mean : "
      << double_to_string_with_commas(crosstalk_mean_, 3)
      << '\n' << indentN << "SPE resolution : "
      << double_to_string_with_commas(spe_resolution_, 3)
      << '\n' << indentN << "Mean charge : "
      << double_to_string_with_commas(1.0 / (1.0 - crosstalk_mean_), 3);
    return stream.str();
  }

#ifndef SWIG
  double_vt vcl_generate_amplitude(VCLRNG& rng) const final {
    const double_vt n = rng.borel_tanner_k1_exp_neg_mu_double(
      double_vt(exp_neg_crosstalk_mean_));
    return n + vcl::sqrt(n) * spe_resolution_ * rng.normal_double();
  }
#endif

  double crosstalk_mean() const { return crosstalk_mean_; }
  double spe_resolution() const { return spe_resolution_; }

private:
  double crosstalk_mean_;
  double exp_neg_crosstalk_mean_;
  double spe_resolution_;
};

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCLDirectionResponse
{
public:
#ifndef SWIG
  using double_vt   = typename VCLArchitecture::double_vt;
#endif

  virtual ~VCLDirectionResponse() {
    // nothing to see here
  }

  virtual double scale() const {
    return 1.0;
  }

  virtual std::string banner() const {
    return "";
  }

#ifndef SWIG
  virtual double_vt detection_probability(const double_vt h_emission,
    const double_vt ux, const double_vt uy, const double_vt uz) const
  {
    return 1.0;
  }
#endif
};

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCLUY1DSplineDirectionResponse:
  public VCLDirectionResponse<VCLArchitecture>
{
public:
#ifndef SWIG
  using double_vt   = typename VCLArchitecture::double_vt;
  using Vector3d_vt = typename VCLArchitecture::Vector3d_vt;
#endif

  VCLUY1DSplineDirectionResponse(calin::math::spline_interpolation::CubicSpline* spline,
    bool adopt_spline): VCLDirectionResponse<VCLArchitecture>(),
    spline_(adopt_spline ? spline : new calin::math::spline_interpolation::CubicSpline(*spline)) { }

  VCLUY1DSplineDirectionResponse(const AngularEfficiency& aeff, bool rescale_to_unity = true,
      double dx_multiplier = 1.0) {
    calin::math::spline_interpolation::CubicSpline spline(aeff.all_xi(), aeff.all_yi());
    spline_ = spline.new_regularized_spline_points_multiplier(dx_multiplier);
    if(rescale_to_unity) {
      scale_ = spline_->ymax();
      spline_->rescale(1.0/scale_);
    }
  }

  virtual ~VCLUY1DSplineDirectionResponse() {
    delete spline_;
  }

  double scale() const final {
    return scale_;
  }

  std::string banner() const final {
    using calin::util::string::double_to_string_with_commas;
    double w0 = spline_->value(1.0);
    double wmax = spline_->x_at_ymax();
    double emax = spline_->value(wmax);
    double whalf = spline_->find(0.5*emax);
    std::ostringstream stream;
    stream << double_to_string_with_commas(w0*100,1)
      << "% at 0 deg; 100% at " << double_to_string_with_commas(std::acos(wmax)/M_PI*180,1)
      << " deg; 50% at " << double_to_string_with_commas(std::acos(whalf)/M_PI*180,1)
      << " deg (scale " << double_to_string_with_commas(scale_,3) << ')';
    return stream.str();
  }

#ifndef SWIG
  double_vt detection_probability(double_vt h_emission,
      const double_vt ux, const double_vt uy, const double_vt uz) const final {
    return spline_->vcl_value<VCLArchitecture>(vcl::abs(uy));
  }
#endif

  const calin::math::spline_interpolation::CubicSpline* spline() const { return spline_; }

private:
  calin::math::spline_interpolation::CubicSpline* spline_ = nullptr;
  double scale_ = 1.0;
};

} } } // namespace calin::simulation::detector_efficiency
