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

#include <util/vcl.hpp>
#include <util/string.hpp>
#include <simulation/detector_efficiency.hpp>

namespace calin { namespace simulation { namespace detector_efficiency {

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
