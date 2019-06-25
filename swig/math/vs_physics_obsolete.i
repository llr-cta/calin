//-*-mode:swig;-*-

/*

   calin/math/vs_physics.i -- Stephen Fegan -- 2015-12-17

   SWIG interface file for calin.simulation.vs_physics - Vec3D, Vec4D, Particle

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.math") vs_physics
%feature(autodoc,2);

%{
#include "math/vs_vec3d.hpp"
#include "math/vs_vec4d.hpp"
#include "math/vs_particle.hpp"
#define SWIG_FILE_WITH_INIT
%}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "common_types.pb.i"
%import "math/rng.i"

%ignore calin::math::vs_physics::Vec3D::operator =;
%ignore calin::math::vs_physics::Vec4D::operator =;

%include "math/vs_vec3d.hpp"
%include "math/vs_vec4d.hpp"
%include "math/vs_particle.hpp"
