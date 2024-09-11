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

%module (package="calin.iact_data") raw_acada_event_data_source_r1v1
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

        RRRRRRRRRRRRRRRRR     1111111                        1111111   
        R::::::::::::::::R   1::::::1                       1::::::1   
        R::::::RRRRRR:::::R 1:::::::1                      1:::::::1   
        RR:::::R     R:::::R111:::::1                      111:::::1   
          R::::R     R:::::R   1::::1vvvvvvv           vvvvvvv1::::1   
          R::::R     R:::::R   1::::1 v:::::v         v:::::v 1::::1   
          R::::RRRRRR:::::R    1::::1  v:::::v       v:::::v  1::::1   
          R:::::::::::::RR     1::::l   v:::::v     v:::::v   1::::l   
          R::::RRRRRR:::::R    1::::l    v:::::v   v:::::v    1::::l   
          R::::R     R:::::R   1::::l     v:::::v v:::::v     1::::l   
          R::::R     R:::::R   1::::l      v:::::v:::::v      1::::l   
          R::::R     R:::::R   1::::l       v:::::::::v       1::::l   
        RR:::::R     R:::::R111::::::111     v:::::::v     111::::::111
        R::::::R     R:::::R1::::::::::1      v:::::v      1::::::::::1
        R::::::R     R:::::R1::::::::::1       v:::v       1::::::::::1
        RRRRRRRR     RRRRRRR111111111111        vvv        111111111111

*/

namespace R1v1 {
  class Event: public google::protobuf::Message { };
  class CameraConfiguration: public google::protobuf::Message { };
}

%typemap(in, numinputs=0) R1v1::Event** CALIN_PROTOBUF_OUTPUT
  (R1v1::Event* temp = nullptr) {
    // typemap(in) R1v1::Event** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    $1 = &temp;
}

%typemap(argout) R1v1::Event** CALIN_PROTOBUF_OUTPUT {
    // typemap(argout) R1v1::Event** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(*$1), $*1_descriptor, SWIG_POINTER_OWN));
}

%apply R1v1::Event** CALIN_PROTOBUF_OUTPUT {
  R1v1::Event** event_out };

%template(DataSource_R1v1) 
  calin::io::data_source::DataSource<R1v1::Event>;
%template(RandomAccessDataSource_R1v1)
  calin::io::data_source::RandomAccessDataSource<R1v1::Event>;

%extend calin::io::data_source::DataSource<R1v1::Event> {
  R1v1::Event* simple_get_next()
  {
    uint64_t unused_seq_index = 0;
    return $self->get_next(unused_seq_index);
  }

  void get_next(uint64_t& seq_index_out, R1v1::Event** event_out)
  {
    seq_index_out = 0;
    *event_out = $self->get_next(seq_index_out);
  }
}

%template(ACADACameraEventDataSource_R1v1) 
  calin::iact_data::acada_data_source::ACADACameraEventDataSource<R1v1::Event>;
%template(ACADACameraEventDataSourceWithRunHeader_R1v1)
  calin::iact_data::acada_data_source::ACADACameraEventDataSourceWithRunHeader<
    R1v1::Event,R1v1::CameraConfiguration>;
%template(ACADACameraEventRandomAccessDataSourceWithRunHeader_R1v1)
  calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
    R1v1::Event,R1v1::CameraConfiguration>;

%template(DataSourceOpener_R1v1)
  calin::io::data_source::DataSourceOpener<
    calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
      R1v1::Event,R1v1::CameraConfiguration> >;

%template(BasicChainedDataSource_ACADACameraEventRandomAccessDataSourceWithRunHeader_R1v1)
  calin::io::data_source::BasicChainedDataSource<
    calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
      R1v1::Event,R1v1::CameraConfiguration> >;

%template(BasicChainedRandomAccessDataSource_ACADACameraEventRandomAccessDataSourceWithRunHeader_R1v1)
  calin::io::data_source::BasicChainedRandomAccessDataSource<
    calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
      R1v1::Event,R1v1::CameraConfiguration> >;

%template(ZFITSACADACameraEventDataSourceOpener_R1v1) 
  calin::iact_data::zfits_acada_data_source::ZFITSACADACameraEventDataSourceOpener<
    R1v1::Event,R1v1::CameraConfiguration>;

%template(ZFITSSingleFileACADACameraEventDataSource_R1v1) 
  calin::iact_data::zfits_acada_data_source::ZFITSSingleFileACADACameraEventDataSource<
    R1v1::Event,R1v1::CameraConfiguration>;
%template(ZFITSACADACameraEventDataSource_R1v1) 
  calin::iact_data::zfits_acada_data_source::ZFITSACADACameraEventDataSource<
    R1v1::Event,R1v1::CameraConfiguration>;
