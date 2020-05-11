/*

   calin/iact_data/nectarcam_ancillary_data.i -- Stephen Fegan -- 2020-05-11

   Classes to extract NectarCAM ancillary data from SQL database

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

%module (package="calin.iact_data") nectarcam_ancillary_data
%feature(autodoc,2);

%{
#include "iact_data/nectarcam_ancillary_data.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%newobject retrieve_nectarcam_ancillary_data(const std::string& db_file,
  int camera_id, int64_t start_time_sec, int64_t end_time_sec);
%newobject retrieve_nectarcam_ancillary_data(const std::string& db_file,
  int camera_id, const std::string& start_time_sec, const std::string& end_time_sec);


%import "iact_data/nectarcam_ancillary_data.pb.i"
%include "iact_data/nectarcam_ancillary_data.hpp"
