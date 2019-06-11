/*

   calin/simulation/pmt.i -- Stephen Fegan -- 2016-03-21

   SWIG interface file for PMT simulation

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

%module (package="calin.simulation") pmt
%feature(autodoc,2);

%{
#include "simulation/pmt.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "math/rng.i"
%import "simulation/pmt.pb.i"

%apply double &OUTPUT { double& x };
%apply double &OUTPUT { double& m };
%apply int &OUTPUT { int& n };
%apply Eigen::VectorXd &OUTPUT { Eigen::VectorXd& x };
%apply Eigen::VectorXd &OUTPUT { Eigen::VectorXd& m };
%apply Eigen::VectorXi &OUTPUT { Eigen::VectorXi& n };

%include "simulation/pmt.hpp"

%clear double& x;
%clear double& m;
%clear int& n;
%clear Eigen::VectorXd& x;
%clear Eigen::VectorXd& m;
%clear Eigen::VectorXi& n;
