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

%module (package="calin.iact_data") raw_acada_event_data_source_l0
%feature(autodoc,2);

%{
#include <iact_data/telescope_data_source.hpp>
#include <iact_data/zfits_acada_data_source.hpp>
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
%include "iact_data/acada_data_source.hpp"
%include "iact_data/zfits_acada_data_source.hpp"

%apply uint64_t& OUTPUT { uint64_t& seq_index_out };
%apply google::protobuf::Arena** CALIN_ARENA_OUTPUT {
  google::protobuf::Arena** arena_out };

/*

              LLLLLLLLLLL                       000000000
              L:::::::::L                     00:::::::::00
              L:::::::::L                   00:::::::::::::00
              LL:::::::LL                  0:::::::000:::::::0
                L:::::L                    0::::::0   0::::::0
                L:::::L                    0:::::0     0:::::0
                L:::::L                    0:::::0     0:::::0
                L:::::L                    0:::::0 000 0:::::0
                L:::::L                    0:::::0 000 0:::::0
                L:::::L                    0:::::0     0:::::0
                L:::::L                    0:::::0     0:::::0
                L:::::L         LLLLLL     0::::::0   0::::::0
              LL:::::::LLLLLLLLL:::::L     0:::::::000:::::::0
              L::::::::::::::::::::::L      00:::::::::::::00
              L::::::::::::::::::::::L        00:::::::::00
              LLLLLLLLLLLLLLLLLLLLLLLL          000000000

*/

namespace ProtoDataModel {
  class CameraEvent: public google::protobuf::Message { };
  class CameraRunHeader: public google::protobuf::Message { };
}

%typemap(in, numinputs=0) ProtoDataModel::CameraEvent** CALIN_PROTOBUF_OUTPUT
  (ProtoDataModel::CameraEvent* temp = nullptr) {
    // typemap(in) ProtoDataModel::CameraEvent** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    $1 = &temp;
}

%typemap(argout) ProtoDataModel::CameraEvent** CALIN_PROTOBUF_OUTPUT {
    // typemap(argout) ProtoDataModel::CameraEvent** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(*$1), $*1_descriptor, SWIG_POINTER_OWN));
}

%apply ProtoDataModel::CameraEvent** CALIN_PROTOBUF_OUTPUT {
  ProtoDataModel::CameraEvent** event_out };

%template(DataSource_L0) 
  calin::io::data_source::DataSource<ProtoDataModel::CameraEvent>;
%template(RandomAccessDataSource_L0)
  calin::io::data_source::RandomAccessDataSource<ProtoDataModel::CameraEvent>;

%extend calin::io::data_source::DataSource<ProtoDataModel::CameraEvent> {
  ProtoDataModel::CameraEvent* simple_get_next()
  {
    uint64_t unused_seq_index = 0;
    return $self->get_next(unused_seq_index);
  }

  void get_next(uint64_t& seq_index_out, ProtoDataModel::CameraEvent** event_out)
  {
    seq_index_out = 0;
    *event_out = $self->get_next(seq_index_out);
  }
}

%template(ACADACameraEventDataSource_L0) 
  calin::iact_data::acada_data_source::ACADACameraEventDataSource<ProtoDataModel::CameraEvent>;
%template(ACADACameraEventDataSourceWithRunHeader_L0)
  calin::iact_data::acada_data_source::ACADACameraEventDataSourceWithRunHeader<
    ProtoDataModel::CameraEvent,ProtoDataModel::CameraRunHeader>;
%template(ACADACameraEventRandomAccessDataSourceWithRunHeader_L0)
  calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
    ProtoDataModel::CameraEvent,ProtoDataModel::CameraRunHeader>;

%template(DataSourceOpener_L0)
  calin::io::data_source::DataSourceOpener<
    calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
      ProtoDataModel::CameraEvent,ProtoDataModel::CameraRunHeader> >;

%template(BasicChainedDataSource_ACADACameraEventRandomAccessDataSourceWithRunHeader_L0)
  calin::io::data_source::BasicChainedDataSource<
    calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
      ProtoDataModel::CameraEvent,ProtoDataModel::CameraRunHeader> >;

%template(BasicChainedRandomAccessDataSource_ACADACameraEventRandomAccessDataSourceWithRunHeader_L0)
  calin::io::data_source::BasicChainedRandomAccessDataSource<
    calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
      ProtoDataModel::CameraEvent,ProtoDataModel::CameraRunHeader> >;

%template(ZFITSACADACameraEventDataSourceOpener_L0) 
  calin::iact_data::zfits_acada_data_source::ZFITSACADACameraEventDataSourceOpener<
    ProtoDataModel::CameraEvent,ProtoDataModel::CameraRunHeader>;

%template(ZFITSSingleFileACADACameraEventDataSource_L0) 
  calin::iact_data::zfits_acada_data_source::ZFITSSingleFileACADACameraEventDataSource<
    ProtoDataModel::CameraEvent,ProtoDataModel::CameraRunHeader>;
%template(ZFITSACADACameraEventDataSource_L0) 
  calin::iact_data::zfits_acada_data_source::ZFITSACADACameraEventDataSource<
    ProtoDataModel::CameraEvent,ProtoDataModel::CameraRunHeader>;
