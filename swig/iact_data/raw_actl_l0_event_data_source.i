/*

   calin/iact_data/raw_actl_event_data_source.i -- Stephen Fegan -- 2016-12-12

   SWIG interface file for calin.io.telescope_data_source

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.iact_data") raw_actl_l0_event_data_source
%feature(autodoc,2);

#define CALIN_NO_ACTL_R1

%{
#include <ProtobufIFits.h>
#include <L0.pb.h>
#include <R1.pb.h>
#include <calin_global_definitions.hpp>
#include <calin_global_config.hpp>
#include <common_types.pb.h>
#include <iact_data/telescope_event.pb.h>
#include <iact_data/instrument_layout.pb.h>
#include <iact_data/nectarcam_data_source.pb.h>
#include <iact_data/lstcam_data_source.pb.h>
#include <iact_data/zfits_actl_data_source.hpp>
#include <iact_data/telescope_data_source.hpp>
#include <iact_data/nectarcam_data_source.hpp>
#include <iact_data/lstcam_data_source.hpp>
#include <iact_data/cta_data_source.hpp>
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%include "calin_global_config.hpp"
%import "calin_global_definitions.i"

%import "google_protobuf.i"
%import "io/data_source.i"

%ignore get_next(uint64_t& seq_index_out);
%ignore get_next(uint64_t& seq_index_out, google::protobuf::Arena** arena);

namespace DataModel {

class CameraEvent: public google::protobuf::Message { };
class CameraRunHeader: public google::protobuf::Message { };
class CTAMessage: public google::protobuf::Message { };

}

%typemap(in, numinputs=0) DataModel::CameraEvent** CALIN_PROTOBUF_OUTPUT
  (DataModel::CameraEvent* temp = nullptr) {
    // typemap(in) DataModel::CameraEvent** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    $1 = &temp;
}

%typemap(argout) DataModel::CameraEvent** CALIN_PROTOBUF_OUTPUT {
    // typemap(argout) DataModel::CameraEvent** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(*$1), $*1_descriptor, SWIG_POINTER_OWN));
}

%apply DataModel::CameraEvent** CALIN_PROTOBUF_OUTPUT {
  DataModel::CameraEvent** event_out };
%apply uint64_t& OUTPUT { uint64_t& seq_index_out };
%apply google::protobuf::Arena** CALIN_ARENA_OUTPUT {
  google::protobuf::Arena** arena_out };

%import "iact_data/telescope_data_source.i"
%import "iact_data/zfits_data_source.pb.i"

%template(ACTL_L0_CameraEventDataSource)
  calin::io::data_source::DataSource<DataModel::CameraEvent>;
%template(ConstACTL_L0_CameraEventDataSource)
  calin::io::data_source::DataSource<const DataModel::CameraEvent>;
%template(ConstACTL_L0_CameraEventDataSink)
  calin::io::data_source::DataSink<const DataModel::CameraEvent>;

%newobject simple_get_next();

%extend calin::io::data_source::DataSource<DataModel::CameraEvent> {

  DataModel::CameraEvent* simple_get_next()
  {
    uint64_t unused_seq_index = 0;
    return $self->get_next(unused_seq_index);
  }

#if 0
  void get_next(uint64_t& seq_index_out, DataModel::CameraEvent** event_out,
    google::protobuf::Arena** arena_out)
  {
    seq_index_out = 0;
    *event_out = $self->get_next(seq_index_out, arena_out);
    if(*event_out != nullptr and *arena_out == nullptr)
      throw std::runtime_error("Memory allocation error: no Arena returned");
  }
#else
  void get_next(uint64_t& seq_index_out, DataModel::CameraEvent** event_out)
  {
    seq_index_out = 0;
    *event_out = $self->get_next(seq_index_out);
  }
#endif

}

%newobject get_run_header();

%include "iact_data/zfits_actl_data_source.hpp"
