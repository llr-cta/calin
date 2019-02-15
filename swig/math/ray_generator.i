/*

   calin/math/ray_generator.i -- Stephen Fegan -- 2017-01-19

   SWIG interface file for calin.math.position_generator,
                           calin.math.direction_generator, and
                           calin.math.ray_generator

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

%module (package="calin.math") ray_generator
%feature(autodoc,2);

%{
#include "math/position_generator.hpp"
#include "math/direction_generator.hpp"
#include "math/ray_generator.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "math/rng.i"
%import "math/ray.i"

%apply Eigen::Vector3d &OUTPUT { Eigen::Vector3d& pos };
%apply double &OUTPUT { double& weight };
%include "math/position_generator.hpp"

%apply double &OUTPUT { double& theta };
%apply double &OUTPUT { double& phi };
%apply Eigen::Vector3d &OUTPUT { Eigen::Vector3d& dir };
%apply Eigen::Matrix3d &OUTPUT { Eigen::Matrix3d& trans_mat };
%include "math/direction_generator.hpp"

%include "math/ray_generator.hpp"
