/*

   calin/iact_data/raw_acada_event_data_source.i -- Stephen Fegan -- 2024-09-09

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

%module (package="calin.iact_data") raw_acada_event_data_source
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
%newobject get_data_stream();

%apply std::vector<std::string>& OUTPUT { std::vector<std::string>& fragment_filenames };
%apply unsigned& OUTPUT { unsigned& num_missing_fragments };

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

namespace calin::iact_data::acada_data_source {
  class ACADA_EventMessage_L0: public google::protobuf::Message { };
  class ACADA_HeaderMessage_L0: public google::protobuf::Message { };
}

%typemap(in, numinputs=0) calin::iact_data::acada_data_source::ACADA_EventMessage_L0** CALIN_PROTOBUF_OUTPUT
  (calin::iact_data::acada_data_source::ACADA_EventMessage_L0* temp = nullptr) {
    // typemap(in) calin::iact_data::acada_data_source::ACADA_EventMessage_L0** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    $1 = &temp;
}

%typemap(argout) calin::iact_data::acada_data_source::ACADA_EventMessage_L0** CALIN_PROTOBUF_OUTPUT {
    // typemap(argout) calin::iact_data::acada_data_source::ACADA_EventMessage_L0** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(*$1), $*1_descriptor, SWIG_POINTER_OWN));
}

%typemap(in, numinputs=0) calin::iact_data::acada_data_source::ACADA_HeaderMessage_L0** CALIN_PROTOBUF_OUTPUT
  (calin::iact_data::acada_data_source::ACADA_HeaderMessage_L0* temp = nullptr) {
    // typemap(in) calin::iact_data::acada_data_source::ACADA_HeaderMessage_L0** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    $1 = &temp;
}

%typemap(argout) calin::iact_data::acada_data_source::ACADA_HeaderMessage_L0** CALIN_PROTOBUF_OUTPUT {
    // typemap(argout) calin::iact_data::acada_data_source::ACADA_HeaderMessage_L0** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(*$1), $*1_descriptor, SWIG_POINTER_OWN));
}

%apply calin::iact_data::acada_data_source::ACADA_EventMessage_L0** CALIN_PROTOBUF_OUTPUT {
  calin::iact_data::acada_data_source::ACADA_EventMessage_L0** message_out };
%apply calin::iact_data::acada_data_source::ACADA_HeaderMessage_L0** CALIN_PROTOBUF_OUTPUT {
  calin::iact_data::acada_data_source::ACADA_HeaderMessage_L0** message_out };

%template(DataSource_ACADAEventMessage_L0)
  calin::io::data_source::DataSource<calin::iact_data::acada_data_source::ACADA_EventMessage_L0>;
%template(DataSource_ACADAHeaderMessage_L0)
  calin::io::data_source::DataSource<calin::iact_data::acada_data_source::ACADA_HeaderMessage_L0>;

%template(RandomAccessDataSource_ACADAEventMessage_L0)
  calin::io::data_source::RandomAccessDataSource<calin::iact_data::acada_data_source::ACADA_EventMessage_L0>;
%template(RandomAccessDataSource_ACADAHeaderMessage_L0)
  calin::io::data_source::RandomAccessDataSource<calin::iact_data::acada_data_source::ACADA_HeaderMessage_L0>;

%extend calin::io::data_source::DataSource<calin::iact_data::acada_data_source::ACADA_EventMessage_L0> {
  calin::iact_data::acada_data_source::ACADA_EventMessage_L0* simple_get_next()
  {
    uint64_t unused_seq_index = 0;
    return $self->get_next(unused_seq_index);
  }

  void get_next(uint64_t& seq_index_out, calin::iact_data::acada_data_source::ACADA_EventMessage_L0** message_out)
  {
    seq_index_out = 0;
    *message_out = $self->get_next(seq_index_out);
  }
}

%extend calin::io::data_source::DataSource<calin::iact_data::acada_data_source::ACADA_HeaderMessage_L0> {
  calin::iact_data::acada_data_source::ACADA_HeaderMessage_L0* simple_get_next()
  {
    uint64_t unused_seq_index = 0;
    return $self->get_next(unused_seq_index);
  }

  void get_next(uint64_t& seq_index_out, calin::iact_data::acada_data_source::ACADA_HeaderMessage_L0** message_out)
  {
    seq_index_out = 0;
    *message_out = $self->get_next(seq_index_out);
  }
}

%template(ACADACameraEventDataSource_L0) 
  calin::iact_data::acada_data_source::ACADACameraEventDataSource<
    calin::iact_data::acada_data_source::ACADA_EventMessage_L0>;
%template(ACADACameraEventDataSourceWithRunHeader_L0)
  calin::iact_data::acada_data_source::ACADACameraEventDataSourceWithRunHeader<
    calin::iact_data::acada_data_source::ACADA_MessageSet_L0>;
%template(ACADACameraEventRandomAccessDataSourceWithRunHeader_L0)
  calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
    calin::iact_data::acada_data_source::ACADA_MessageSet_L0>;

%template(DataSourceOpener_L0)
  calin::io::data_source::DataSourceOpener<
    calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
      calin::iact_data::acada_data_source::ACADA_MessageSet_L0> >;

%template(BasicChainedDataSource_ACADACameraEventRandomAccessDataSourceWithRunHeader_L0)
  calin::io::data_source::BasicChainedDataSource<
    calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
      calin::iact_data::acada_data_source::ACADA_MessageSet_L0> >;

%template(BasicChainedRandomAccessDataSource_ACADACameraEventRandomAccessDataSourceWithRunHeader_L0)
  calin::io::data_source::BasicChainedRandomAccessDataSource<
    calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
      calin::iact_data::acada_data_source::ACADA_MessageSet_L0> >;

%template(ZFITSACADACameraEventDataSourceOpener_L0) 
  calin::iact_data::zfits_acada_data_source::ZFITSACADACameraEventDataSourceOpener<
    calin::iact_data::acada_data_source::ACADA_MessageSet_L0>;

%template(ZFITSSingleFileACADACameraEventDataSource_L0) 
  calin::iact_data::zfits_acada_data_source::ZFITSSingleFileACADACameraEventDataSource<
    calin::iact_data::acada_data_source::ACADA_MessageSet_L0>;
%template(ZFITSACADACameraEventDataSource_L0) 
  calin::iact_data::zfits_acada_data_source::ZFITSACADACameraEventDataSource<
    calin::iact_data::acada_data_source::ACADA_MessageSet_L0>;

%template(ZFITSSingleFileSingleMessageDataSource_ACADAEventMessage_L0)
  calin::iact_data::zfits_acada_data_source::
    ZFITSSingleFileSingleMessageDataSource<calin::iact_data::acada_data_source::ACADA_EventMessage_L0>;
%template(ZFITSSingleFileSingleMessageDataSource_ACADAHeaderMessage_L0)
  calin::iact_data::zfits_acada_data_source::
    ZFITSSingleFileSingleMessageDataSource<calin::iact_data::acada_data_source::ACADA_HeaderMessage_L0>;

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

namespace calin::iact_data::acada_data_source {
  class ACADA_EventMessage_R1v0: public google::protobuf::Message { };
  class ACADA_HeaderMessage_R1v0: public google::protobuf::Message { };
}

%typemap(in, numinputs=0) calin::iact_data::acada_data_source::ACADA_EventMessage_R1v0** CALIN_PROTOBUF_OUTPUT
  (calin::iact_data::acada_data_source::ACADA_EventMessage_R1v0* temp = nullptr) {
    // typemap(in) calin::iact_data::acada_data_source::ACADA_EventMessage_R1v0** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    $1 = &temp;
}

%typemap(argout) calin::iact_data::acada_data_source::ACADA_EventMessage_R1v0** CALIN_PROTOBUF_OUTPUT {
    // typemap(argout) calin::iact_data::acada_data_source::ACADA_EventMessage_R1v0** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(*$1), $*1_descriptor, SWIG_POINTER_OWN));
}

%typemap(in, numinputs=0) calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v0** CALIN_PROTOBUF_OUTPUT
  (calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v0* temp = nullptr) {
    // typemap(in) calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v0** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    $1 = &temp;
}

%typemap(argout) calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v0** CALIN_PROTOBUF_OUTPUT {
    // typemap(argout) calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v0** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(*$1), $*1_descriptor, SWIG_POINTER_OWN));
}

%apply calin::iact_data::acada_data_source::ACADA_EventMessage_R1v0** CALIN_PROTOBUF_OUTPUT {
  calin::iact_data::acada_data_source::ACADA_EventMessage_R1v0** message_out };
%apply calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v0** CALIN_PROTOBUF_OUTPUT {
  calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v0** message_out };

%template(DataSource_ACADAEventMessage_R1v0)
  calin::io::data_source::DataSource<calin::iact_data::acada_data_source::ACADA_EventMessage_R1v0>;
%template(DataSource_ACADAHeaderMessage_R1v0)
  calin::io::data_source::DataSource<calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v0>;

%template(RandomAccessDataSource_ACADAEventMessage_R1v0)
  calin::io::data_source::RandomAccessDataSource<calin::iact_data::acada_data_source::ACADA_EventMessage_R1v0>;
%template(RandomAccessDataSource_ACADAHeaderMessage_R1v0)
  calin::io::data_source::RandomAccessDataSource<calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v0>;

%extend calin::io::data_source::DataSource<calin::iact_data::acada_data_source::ACADA_EventMessage_R1v0> {
  calin::iact_data::acada_data_source::ACADA_EventMessage_R1v0* simple_get_next()
  {
    uint64_t unused_seq_index = 0;
    return $self->get_next(unused_seq_index);
  }

  void get_next(uint64_t& seq_index_out, calin::iact_data::acada_data_source::ACADA_EventMessage_R1v0** message_out)
  {
    seq_index_out = 0;
    *message_out = $self->get_next(seq_index_out);
  }
}

%extend calin::io::data_source::DataSource<calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v0> {
  calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v0* simple_get_next()
  {
    uint64_t unused_seq_index = 0;
    return $self->get_next(unused_seq_index);
  }

  void get_next(uint64_t& seq_index_out, calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v0** message_out)
  {
    seq_index_out = 0;
    *message_out = $self->get_next(seq_index_out);
  }
}

%template(ACADACameraEventDataSource_R1v0) 
  calin::iact_data::acada_data_source::ACADACameraEventDataSource<
    calin::iact_data::acada_data_source::ACADA_EventMessage_R1v0>;
%template(ACADACameraEventDataSourceWithRunHeader_R1v0)
  calin::iact_data::acada_data_source::ACADACameraEventDataSourceWithRunHeader<
    calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0>;
%template(ACADACameraEventRandomAccessDataSourceWithRunHeader_R1v0)
  calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
    calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0>;

%template(DataSourceOpener_R1v0)
  calin::io::data_source::DataSourceOpener<
    calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
      calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0> >;

%template(BasicChainedDataSource_ACADACameraEventRandomAccessDataSourceWithRunHeader_R1v0)
  calin::io::data_source::BasicChainedDataSource<
    calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
      calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0> >;

%template(BasicChainedRandomAccessDataSource_ACADACameraEventRandomAccessDataSourceWithRunHeader_R1v0)
  calin::io::data_source::BasicChainedRandomAccessDataSource<
    calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
      calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0> >;

%template(ZFITSACADACameraEventDataSourceOpener_R1v0) 
  calin::iact_data::zfits_acada_data_source::ZFITSACADACameraEventDataSourceOpener<
    calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0>;

%template(ZFITSSingleFileACADACameraEventDataSource_R1v0) 
  calin::iact_data::zfits_acada_data_source::ZFITSSingleFileACADACameraEventDataSource<
    calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0>;
%template(ZFITSACADACameraEventDataSource_R1v0) 
  calin::iact_data::zfits_acada_data_source::ZFITSACADACameraEventDataSource<
    calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0>;

%template(ZFITSSingleFileSingleMessageDataSource_ACADAEventMessage_R1v0)
  calin::iact_data::zfits_acada_data_source::
    ZFITSSingleFileSingleMessageDataSource<calin::iact_data::acada_data_source::ACADA_EventMessage_R1v0>;
%template(ZFITSSingleFileSingleMessageDataSource_ACADAHeaderMessage_R1v0)
  calin::iact_data::zfits_acada_data_source::
    ZFITSSingleFileSingleMessageDataSource<calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v0>;

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

namespace calin::iact_data::acada_data_source {
  class ACADA_EventMessage_R1v1: public google::protobuf::Message { };
  class ACADA_HeaderMessage_R1v1: public google::protobuf::Message { };
  class ACADA_DataStreamMessage_R1v1: public google::protobuf::Message { };
}

%typemap(in, numinputs=0) calin::iact_data::acada_data_source::ACADA_EventMessage_R1v1** CALIN_PROTOBUF_OUTPUT
  (calin::iact_data::acada_data_source::ACADA_EventMessage_R1v1* temp = nullptr) {
    // typemap(in) calin::iact_data::acada_data_source::ACADA_EventMessage_R1v1** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    $1 = &temp;
}

%typemap(argout) calin::iact_data::acada_data_source::ACADA_EventMessage_R1v1** CALIN_PROTOBUF_OUTPUT {
    // typemap(argout) calin::iact_data::acada_data_source::ACADA_EventMessage_R1v1** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(*$1), $*1_descriptor, SWIG_POINTER_OWN));
}

%typemap(in, numinputs=0) calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v1** CALIN_PROTOBUF_OUTPUT
  (calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v1* temp = nullptr) {
    // typemap(in) calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v1** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    $1 = &temp;
}

%typemap(argout) calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v1** CALIN_PROTOBUF_OUTPUT {
    // typemap(argout) calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v1** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(*$1), $*1_descriptor, SWIG_POINTER_OWN));
}

%typemap(in, numinputs=0) calin::iact_data::acada_data_source::ACADA_DataStreamMessage_R1v1** CALIN_PROTOBUF_OUTPUT
  (calin::iact_data::acada_data_source::ACADA_DataStreamMessage_R1v1* temp = nullptr) {
    // typemap(in) calin::iact_data::acada_data_source::ACADA_DataStreamMessage_R1v1** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    $1 = &temp;
}

%typemap(argout) calin::iact_data::acada_data_source::ACADA_DataStreamMessage_R1v1** CALIN_PROTOBUF_OUTPUT {
    // typemap(argout) calin::iact_data::acada_data_source::ACADA_DataStreamMessage_R1v1** CALIN_PROTOBUF_OUTPUT - raw_actl_event_data_source.i
    %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(*$1), $*1_descriptor, SWIG_POINTER_OWN));
}

%apply calin::iact_data::acada_data_source::ACADA_EventMessage_R1v1** CALIN_PROTOBUF_OUTPUT {
  calin::iact_data::acada_data_source::ACADA_EventMessage_R1v1** message_out };
%apply calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v1** CALIN_PROTOBUF_OUTPUT {
  calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v1** message_out };
%apply calin::iact_data::acada_data_source::ACADA_DataStreamMessage_R1v1** CALIN_PROTOBUF_OUTPUT {
  calin::iact_data::acada_data_source::ACADA_DataStreamMessage_R1v1** message_out };

%template(DataSource_ACADAEventMessage_R1v1)
  calin::io::data_source::DataSource<calin::iact_data::acada_data_source::ACADA_EventMessage_R1v1>;
%template(DataSource_ACADAHeaderMessage_R1v1)
  calin::io::data_source::DataSource<calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v1>;
%template(DataSource_ACADADataStreamMessage_R1v1)
  calin::io::data_source::DataSource<calin::iact_data::acada_data_source::ACADA_DataStreamMessage_R1v1>;

%template(RandomAccessDataSource_ACADAEventMessage_R1v1)
  calin::io::data_source::RandomAccessDataSource<calin::iact_data::acada_data_source::ACADA_EventMessage_R1v1>;
%template(RandomAccessDataSource_ACADAHeaderMessage_R1v1)
  calin::io::data_source::RandomAccessDataSource<calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v1>;
%template(RandomAccessDataSource_ACADADataStreamMessage_R1v1)
  calin::io::data_source::RandomAccessDataSource<calin::iact_data::acada_data_source::ACADA_DataStreamMessage_R1v1>;

%extend calin::io::data_source::DataSource<calin::iact_data::acada_data_source::ACADA_EventMessage_R1v1> {
  calin::iact_data::acada_data_source::ACADA_EventMessage_R1v1* simple_get_next()
  {
    uint64_t unused_seq_index = 0;
    return $self->get_next(unused_seq_index);
  }

  void get_next(uint64_t& seq_index_out, calin::iact_data::acada_data_source::ACADA_EventMessage_R1v1** message_out)
  {
    seq_index_out = 0;
    *message_out = $self->get_next(seq_index_out);
  }
}

%extend calin::io::data_source::DataSource<calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v1> {
  calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v1* simple_get_next()
  {
    uint64_t unused_seq_index = 0;
    return $self->get_next(unused_seq_index);
  }

  void get_next(uint64_t& seq_index_out, calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v1** message_out)
  {
    seq_index_out = 0;
    *message_out = $self->get_next(seq_index_out);
  }
}

%extend calin::io::data_source::DataSource<calin::iact_data::acada_data_source::ACADA_DataStreamMessage_R1v1> {
  calin::iact_data::acada_data_source::ACADA_DataStreamMessage_R1v1* simple_get_next()
  {
    uint64_t unused_seq_index = 0;
    return $self->get_next(unused_seq_index);
  }

  void get_next(uint64_t& seq_index_out, calin::iact_data::acada_data_source::ACADA_DataStreamMessage_R1v1** message_out)
  {
    seq_index_out = 0;
    *message_out = $self->get_next(seq_index_out);
  }
}

%template(ACADACameraEventDataSource_R1v1) 
  calin::iact_data::acada_data_source::ACADACameraEventDataSource<
    calin::iact_data::acada_data_source::ACADA_EventMessage_R1v1>;
%template(ACADACameraEventDataSourceWithRunHeader_R1v1)
  calin::iact_data::acada_data_source::ACADACameraEventDataSourceWithRunHeader<
    calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1>;
%template(ACADACameraEventRandomAccessDataSourceWithRunHeader_R1v1)
  calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
    calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1>;

%template(DataSourceOpener_R1v1)
  calin::io::data_source::DataSourceOpener<
    calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
      calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1> >;

%template(BasicChainedDataSource_ACADACameraEventRandomAccessDataSourceWithRunHeader_R1v1)
  calin::io::data_source::BasicChainedDataSource<
    calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
      calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1> >;

%template(BasicChainedRandomAccessDataSource_ACADACameraEventRandomAccessDataSourceWithRunHeader_R1v1)
  calin::io::data_source::BasicChainedRandomAccessDataSource<
    calin::iact_data::acada_data_source::ACADACameraEventRandomAccessDataSourceWithRunHeader<
      calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1> >;

%template(ZFITSACADACameraEventDataSourceOpener_R1v1) 
  calin::iact_data::zfits_acada_data_source::ZFITSACADACameraEventDataSourceOpener<
    calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1>;

%template(ZFITSSingleFileACADACameraEventDataSource_R1v1) 
  calin::iact_data::zfits_acada_data_source::ZFITSSingleFileACADACameraEventDataSource<
    calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1>;
%template(ZFITSACADACameraEventDataSource_R1v1) 
  calin::iact_data::zfits_acada_data_source::ZFITSACADACameraEventDataSource<
    calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1>;

%template(ZFITSSingleFileSingleMessageDataSource_ACADAEventMessage_R1v1)
  calin::iact_data::zfits_acada_data_source::
    ZFITSSingleFileSingleMessageDataSource<calin::iact_data::acada_data_source::ACADA_EventMessage_R1v1>;
%template(ZFITSSingleFileSingleMessageDataSource_ACADAHeaderMessage_R1v1)
  calin::iact_data::zfits_acada_data_source::
    ZFITSSingleFileSingleMessageDataSource<calin::iact_data::acada_data_source::ACADA_HeaderMessage_R1v1>;
%template(ZFITSSingleFileSingleMessageDataSource_ACADADataStreamMessage_R1v1)
  calin::iact_data::zfits_acada_data_source::
    ZFITSSingleFileSingleMessageDataSource<calin::iact_data::acada_data_source::ACADA_DataStreamMessage_R1v1>;
