/*

   calin/iact_data/algorithms.i -- Stephen Fegan -- 2021-07-23

   SWIG interface file for IACT algorithms

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

%module (package="calin.iact_data") algorithms
%feature(autodoc,2);

%{
#include "iact_data/algorithms.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%apply Eigen::VectorXi& OUTPUT { Eigen::VectorXi& channel_island_id };
%apply Eigen::VectorXi& OUTPUT { Eigen::VectorXi& island_count };

/* %import "iact_data/algorithms.pb.i" */
%include "iact_data/algorithms.hpp"
