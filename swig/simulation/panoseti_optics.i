//-*-mode:swig;-*-

/*

   calin/simulation/panoseti_optics.i -- Stephen Fegan -- 2036-08-04

   SWIG interface file for calin.simulation.panoseti_optics

   Copyright 2026, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.simulation") panoseti_optics
%feature(autodoc,2);

%{
#include "simulation/vcl_panoseti_raytracer.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "math/ray.i"
%import "math/rng.i"
%import "util/vcl.i"

%import "simulation/panoseti_optics.pb.i"

%apply Eigen::VectorXd& OUTPUT { Eigen::VectorXd& x_out, Eigen::VectorXd& y_out, Eigen::VectorXd& t_out };

namespace calin { namespace simulation { namespace vcl_raytracer {

  template<typename VCLReal> class alignas(VCLReal::vec_bytes) VCLPanosetiThinLensScopeRayTracer: public VCLReal
  {
  public:
    using real_t = typename VCLReal::real_t;
    using RNG = calin::math::rng::VCLRealRNG<VCLReal>;
    VCLPanosetiThinLensScopeRayTracer(const calin::ix::simulation::panoseti_optics::ArrayParameters& array_params,
        unsigned scope_id, const calin::math::spline_interpolation::CubicSpline& lens_refractive_index_spline,
        real_t air_refractive_index = 1.0,
        RNG* rng = nullptr, bool adopt_rng = false);
    ~VCLPanosetiThinLensScopeRayTracer();
    bool point_telescope(const Eigen::Vector3d& v);
    bool point_telescope_az_el(const double az_rad, const double el_rad);
    bool point_telescope_az_el_phi(double az_rad, double el_rad, double phi_rad);
    unsigned monochromatic_psf(Eigen::VectorXd& x_out, Eigen::VectorXd& y_out, Eigen::VectorXd& t_out, unsigned nray,
      double photon_energy_ev,
      double theta = 0, double phi = 0, double distance = std::numeric_limits<double>::infinity(), double radius = 0);
  };

} } } // namespace calin::simulation::vcl_raytracer

%include "simulation/vcl_panoseti_raytracer.hpp" // nothing is defined here for SWIG

// %template (VCLPanosetiThinLensScopeRayTracerDouble128) 
//   calin::simulation::vcl_raytracer::VCLPanosetiThinLensScopeRayTracer<calin::util::vcl::VCL128DoubleReal>;
%template (VCLPanosetiThinLensScopeRayTracerDouble256) 
  calin::simulation::vcl_raytracer::VCLPanosetiThinLensScopeRayTracer<calin::util::vcl::VCL256DoubleReal>;
// %template (VCLPanosetiThinLensScopeRayTracerDouble512) 
//   calin::simulation::vcl_raytracer::VCLPanosetiThinLensScopeRayTracer<calin::util::vcl::VCL512DoubleReal>;

// %template (VCLPanosetiThinLensScopeRayTracerFloat128) 
//   calin::simulation::vcl_raytracer::VCLPanosetiThinLensScopeRayTracer<calin::util::vcl::VCL128FloatReal>;
// %template (VCLPanosetiThinLensScopeRayTracerFloat256) 
//   calin::simulation::vcl_raytracer::VCLPanosetiThinLensScopeRayTracer<calin::util::vcl::VCL256FloatReal>;
// %template (VCLPanosetiThinLensScopeRayTracerFloat512) 
//   calin::simulation::vcl_raytracer::VCLPanosetiThinLensScopeRayTracer<calin::util::vcl::VCL512FloatReal>;
