/*

   calin/iact_data/actl_event_decoder.i -- Stephen Fegan -- 2018-11-20

   SWIG interface file for calin.iact_data.acada_event_decoder

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.iact_data") acada_event_decoder
%feature(autodoc,2);

%{
#include <calin_global_definitions.hpp>
#include <calin_global_config.hpp>
#include <iact_data/acada_event_decoder.hpp>
#include <iact_data/instrument_layout.pb.h>
#include <iact_data/nectarcam_acada_event_decoder.hpp>
#include <iact_data/lstcam_acada_event_decoder.hpp>
#include <iact_data/unified_acada_event_decoder.hpp>
#include <iact_data/cta_acada_event_decoder.hpp>
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import <iact_data/telescope_event.pb.i>

%import "raw_acada_event_data_source.i"

%include <iact_data/acada_event_decoder.hpp>

%template(ACADACameraEventDecoder_L0)
  calin::iact_data::acada_event_decoder::ACADACameraEventDecoder<
    calin::iact_data::acada_data_source::ACADA_MessageSet_L0>;

%template(ACADACameraEventDecoder_R1v0) 
  calin::iact_data::acada_event_decoder::ACADACameraEventDecoder<
    calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0>;

%template(ACADACameraEventDecoder_R1v1) 
  calin::iact_data::acada_event_decoder::ACADACameraEventDecoder<
    calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1>;

%import <iact_data/nectarcam_configuration.pb.i>
%include <iact_data/nectarcam_acada_event_decoder.hpp>
%include <iact_data/lstcam_acada_event_decoder.hpp>
%include <iact_data/unified_acada_event_decoder.hpp>
%include <iact_data/cta_acada_event_decoder.hpp>
