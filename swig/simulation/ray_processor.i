//-*-mode:swig;-*-

/*

   calin/simulation/ray_processor.i -- Stephen Fegan -- 2017-01-19

   SWIG interface file for calin.simulation.ray_processor

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.simulation") ray_processor
%feature(autodoc,2);

%{
#include "simulation/pe_processor.hpp"
#include "simulation/ray_processor.hpp"
#include "simulation/vso_ray_processor.hpp"
#include "simulation/vcl_ray_processor.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "math/ray.i"
%import "math/ray_generator.i"
%import "math/moments_calc.i"

%include "simulation/pe_processor.hpp"
%newobject *::add_processor_and_pe_visitor;
%include "simulation/ray_processor.hpp"
%include "simulation/vso_ray_processor.hpp"

%template(StdVectorRayProcessorDetectorSphere)
  std::vector<calin::simulation::ray_processor::RayProcessorDetectorSphere>;

%include "simulation/vcl_ray_processor.hpp"

%template (VCLRayTracerRayProcessorDouble128) calin::simulation::vcl_ray_processor::VCLRayTracerRayProcessorDouble<calin::util::vcl::VCL128Architecture>;
%template (VCLRayTracerRayProcessorDouble256) calin::simulation::vcl_ray_processor::VCLRayTracerRayProcessorDouble<calin::util::vcl::VCL256Architecture>;
%template (VCLRayTracerRayProcessorDouble512) calin::simulation::vcl_ray_processor::VCLRayTracerRayProcessorDouble<calin::util::vcl::VCL512Architecture>;

%template (VCLRayTracerRayProcessorFloat128) calin::simulation::vcl_ray_processor::VCLRayTracerRayProcessorFloat<calin::util::vcl::VCL128Architecture>;
%template (VCLRayTracerRayProcessorFloat256) calin::simulation::vcl_ray_processor::VCLRayTracerRayProcessorFloat<calin::util::vcl::VCL256Architecture>;
%template (VCLRayTracerRayProcessorFloat512) calin::simulation::vcl_ray_processor::VCLRayTracerRayProcessorFloat<calin::util::vcl::VCL512Architecture>;
