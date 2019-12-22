/*

   calin/iact_data/actl_event_decoder.i -- Stephen Fegan -- 2018-11-20

   SWIG interface file for calin.iact_data.actl_event_decoder

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

%module (package="calin.iact_data") actl_event_decoder
%feature(autodoc,2);

%{
#include <ProtobufIFits.h>
#include <L0.pb.h>
#include <R1.pb.h>
#include <calin_global_definitions.hpp>
#include <calin_global_config.hpp>
#include <common_types.pb.h>
#include <iact_data/nectarcam_data_source.pb.h>
#include <iact_data/lstcam_data_source.pb.h>
#include <iact_data/actl_event_decoder.hpp>
#include <iact_data/nectarcam_actl_event_decoder.hpp>
#include <iact_data/lstcam_actl_event_decoder.hpp>
#include <iact_data/cta_actl_event_decoder.hpp>
#include <iact_data/cta_data_source.hpp>
#include <iact_data/actl_event_decoder.hpp>
#include <iact_data/nectarcam_data_source.hpp>
#include <iact_data/lstcam_data_source.hpp>
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%include "calin_global_config.hpp"
%import "calin_global_definitions.i"

%import "iact_data/telescope_data_source.i"

%ignore get_next(uint64_t& seq_index_out);
%ignore get_next(uint64_t& seq_index_out, google::protobuf::Arena** arena);

%include <iact_data/actl_event_decoder.hpp>

%import <pattern/delegation.hpp>
%template (DelegatorACTL_R1_CameraEventDecoder)
  calin::pattern::delegation::Delegator<
    calin::iact_data::actl_event_decoder::ACTL_R1_CameraEventDecoder>;
%template (DelegatorACTL_L0_CameraEventDecoder)
  calin::pattern::delegation::Delegator<
    calin::iact_data::actl_event_decoder::ACTL_L0_CameraEventDecoder>;

%include <iact_data/nectarcam_actl_event_decoder.hpp>
%include <iact_data/lstcam_actl_event_decoder.hpp>
%include <iact_data/cta_actl_event_decoder.hpp>
