/*

  calin/calib/pmt_ses_models.i -- Stephen Fegan -- 2017-04-24

  PMT single-electron spectrum models

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

%module (package="calin.calib") pmt_ses_models

%{
#include "math/log_quadratic_spline_pdf_1d.hpp"
#include "calib/pmt_ses_models.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "math/function.i"
%import "math/histogram.i"

%import "math/function.i"
%import "math/pdf_1d.i"

%include "calib/pmt_ses_models.hpp"
