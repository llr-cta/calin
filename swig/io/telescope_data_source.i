/*

   calin/io/telescope_data_source.i -- Stephen Fegan -- 2016-01-14

   SWIG interface file for calin.io.telescope_data_source

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/


%module (package="calin.io") telescope_data_source

%{
#include "io/telescope_data_source.hpp"
#include "io/nectarcam_data_source.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "numpy.i"
%include "stdint.i"
%include "calin_typemaps.i"

%import "calin_global_definitions.i"

%newobject calin::io::telescope_data_source::TelescopeDataSource::getNextEvent();
%newobject calin::io::nectarcam_data_source::NectarCamZFITSDataSource::getNextEvent();

%include "io/telescope_data_source.hpp"
%include "io/nectarcam_data_source.hpp"
