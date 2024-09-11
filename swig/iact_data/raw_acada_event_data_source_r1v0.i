/*

   calin/iact_data/raw_actl_r1_event_data_source.i -- Stephen Fegan -- 2024-09-09

   SWIG interface file for calin.io.raw_acada_event_data_source

   Copyright 2024, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.iact_data") raw_acada_event_data_source_r1v0
%feature(autodoc,2);

%{
#include <iact_data/telescope_data_source.hpp>
#include <iact_data/zfits_acada_data_source.hpp>
// #include <iact_data/nectarcam_data_source.hpp>
// #include <iact_data/lstcam_data_source.hpp>
// #include <iact_data/cta_data_source.hpp>
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
%import "iact_data/zfits_data_source.pb.i"

%ignore get_next(uint64_t& seq_index_out);
%ignore get_next(uint64_t& seq_index_out, google::protobuf::Arena** arena);

%newobject simple_get_next();
%newobject get_run_header();

%include "io/chained_data_source.hpp"
// %import "iact_data/telescope_data_source.i"
%include "iact_data/acada_data_source.hpp"
%include "iact_data/zfits_acada_data_source.hpp"

%apply uint64_t& OUTPUT { uint64_t& seq_index_out };
%apply google::protobuf::Arena** CALIN_ARENA_OUTPUT {
  google::protobuf::Arena** arena_out };

/*

    RRRRRRRRRRRRRRRRR     1111111                              000000000     
    R::::::::::::::::R   1::::::1                            00:::::::::00   
    R::::::RRRRRR:::::R 1:::::::1                          00:::::::::::::00 
    RR:::::R     R:::::R111:::::1                         0:::::::000:::::::0
      R::::R     R:::::R   1::::1vvvvvvv           vvvvvvv0::::::0   0::::::0
      R::::R     R:::::R   1::::1 v:::::v         v:::::v 0:::::0     0:::::0
      R::::RRRRRR:::::R    1::::1  v:::::v       v:::::v  0:::::0     0:::::0
      R:::::::::::::RR     1::::l   v:::::v     v:::::v   0:::::0 000 0:::::0
      R::::RRRRRR:::::R    1::::l    v:::::v   v:::::v    0:::::0 000 0:::::0
      R::::R     R:::::R   1::::l     v:::::v v:::::v     0:::::0     0:::::0
      R::::R     R:::::R   1::::l      v:::::v:::::v      0:::::0     0:::::0
      R::::R     R:::::R   1::::l       v:::::::::v       0::::::0   0::::::0
    RR:::::R     R:::::R111::::::111     v:::::::v        0:::::::000:::::::0
    R::::::R     R:::::R1::::::::::1      v:::::v          00:::::::::::::00 
    R::::::R     R:::::R1::::::::::1       v:::v             00:::::::::00   
    RRRRRRRR     RRRRRRR111111111111        vvv                000000000     

*/

namespace ProtoR1 {
  class CameraEvent: public google::protobuf::Message { };
  class CameraConfiguration: public google::protobuf::Message { };
}

%typemap(in, numinputs=0) ProtoR1::CameraEvent** CALIN_PROTOBUF_OUTPUT
  (ProtoR1::CameraEvent* temp = nullptr) {
    // typemap(in) ProtoR1::CameraEvent** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    $1 = &temp;
}

%typemap(argout) ProtoR1::CameraEvent** CALIN_PROTOBUF_OUTPUT {
    // typemap(argout) ProtoR1::CameraEvent** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(*$1), $*1_descriptor, SWIG_POINTER_OWN));
}

%apply ProtoR1::CameraEvent** CALIN_PROTOBUF_OUTPUT {
  ProtoR1::CameraEvent** event_out };

%template(DataSource_R1v0) 
  calin::io::data_source::DataSource<ProtoR1::CameraEvent>;
%template(RandomAccessDataSource_R1v0)
  calin::io::data_source::RandomAccessDataSource<ProtoR1::CameraEvent>;

%extend calin::io::data_source::DataSource<ProtoR1::CameraEvent> {
  ProtoR1::CameraEvent* simple_get_next()
  {
    uint64_t unused_seq_index = 0;
    return $self->get_next(unused_seq_index);
  }

  void get_next(uint64_t& seq_index_out, ProtoR1::CameraEvent** event_out)
  {
    seq_index_out = 0;
    *event_out = $self->get_next(seq_index_out);
  }
}

%template(ACADACameraEventDataSource_R1v0) 
  calin::iact_data::acada_data_source::ACADACameraEventDataSource<ProtoR1::CameraEvent>;
%template(ACADACameraEventDataSourceWithRunHeader_R1v0)
  calin::iact_data::acada_data_source::ACADACameraEventDataSourceWithRunHeader<
    ProtoR1::CameraEvent,ProtoR1::CameraConfiguration>;
%template(ACADACameraEventRandomAccessDataSourceWithRunHeader_R1v0)
  calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
    ProtoR1::CameraEvent,ProtoR1::CameraConfiguration>;

%template(DataSourceOpener_R1v0)
  calin::io::data_source::DataSourceOpener<
    calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
      ProtoR1::CameraEvent,ProtoR1::CameraConfiguration> >;

%template(BasicChainedDataSource_ACADACameraEventRandomAccessDataSourceWithRunHeader_R1v0)
  calin::io::data_source::BasicChainedDataSource<
    calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
      ProtoR1::CameraEvent,ProtoR1::CameraConfiguration> >;

%template(BasicChainedRandomAccessDataSource_ACADACameraEventRandomAccessDataSourceWithRunHeader_R1v0)
  calin::io::data_source::BasicChainedRandomAccessDataSource<
    calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
      ProtoR1::CameraEvent,ProtoR1::CameraConfiguration> >;

%template(ZFITSACADACameraEventDataSourceOpener_R1v0) 
  calin::iact_data::zfits_acada_data_source::ZFITSACADACameraEventDataSourceOpener<
    ProtoR1::CameraEvent,ProtoR1::CameraConfiguration>;

%template(ZFITSSingleFileACADACameraEventDataSource_R1v0) 
  calin::iact_data::zfits_acada_data_source::ZFITSSingleFileACADACameraEventDataSource<
    ProtoR1::CameraEvent,ProtoR1::CameraConfiguration>;
%template(ZFITSACADACameraEventDataSource_R1v0) 
  calin::iact_data::zfits_acada_data_source::ZFITSACADACameraEventDataSource<
    ProtoR1::CameraEvent,ProtoR1::CameraConfiguration>;
