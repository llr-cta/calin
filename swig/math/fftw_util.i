/*

   calin/math/fftw_util.i -- Stephen Fegan -- 2017-05-19

   SWIG interface file for calin.math.fftw_util

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

%module (package="calin.math") fftw_util
%feature(autodoc,2);

%{
#include "math/fftw_util.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "math/fftw_util.pb.i"

%apply Eigen::VectorXd &OUTPUT { Eigen::VectorXd& oivec };
%apply Eigen::VectorXd &OUTPUT { Eigen::VectorXd& oqvec };

%include "math/fftw_util.hpp"
