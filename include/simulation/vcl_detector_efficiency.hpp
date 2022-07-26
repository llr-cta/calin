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

#include <util/vcl.hpp>
#include <simulation/detector_efficiency.hpp>

namespace calin { namespace simulation { namespace detector_efficiency {

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCLDirectionResponse
{
public:
#ifndef SWIG
  using double_vt   = typename VCLArchitecture::double_vt;
  using Vector3d_vt = typename VCLArchitecture::Vector3d_vt;
#endif

  virtual ~VCLDirectionResponse() {
    // nothing to see here
  }

  virtual double scale() const {
    return 1.0;
  }

#ifndef SWIG
  virtual double_vt detection_probability(double_vt h_emission, const Vector3d_vt& u) const {
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

#ifndef SWIG
  double_vt detection_probability(double_vt h_emission, const Vector3d_vt& u) const final {
    return spline_->vcl_value<VCLArchitecture>(vcl::abs(u.y()));
  }
#endif

  const calin::math::spline_interpolation::CubicSpline* spline() const { return spline_; }

private:
  calin::math::spline_interpolation::CubicSpline* spline_ = nullptr;
  double scale_ = 1.0;
};

} } } // namespace calin::simulation::detector_efficiency
