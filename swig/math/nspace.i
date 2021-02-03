/*

   calin/math/nspace.i -- Stephen Fegan -- 2020-11-27

   SWIG interface file for calin.math.nspace

   Copyright 2020, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.math") nspace
%feature(autodoc,2);

%{
#include "math/nspace.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%apply double &OUTPUT { double& w0 };
%apply Eigen::VectorXd &OUTPUT { Eigen::VectorXd& w1 };
%apply Eigen::VectorXd &OUTPUT { Eigen::VectorXd& x_out };
%apply int64_t &OUTPUT { int64_t& array_index, int64_t& block_index };

%include "math/nspace.hpp"
