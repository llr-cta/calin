/*

   calin/iact_data/telescope_data_source.i -- Stephen Fegan -- 2016-01-14

   SWIG interface file for calin.io.telescope_data_source

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.iact_data") telescope_data_source
%feature(autodoc,2);

%{
#include <calin_global_config.hpp>
#include <iact_data/telescope_data_source.hpp>
#include <iact_data/nectarcam_data_source.hpp>
#include <iact_data/nectarcam_configuration.hpp>
#include <iact_data/lstcam_data_source.hpp>
#include <iact_data/cta_data_source.hpp>
using namespace calin::io;
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%include "calin_global_config.hpp"
%import "calin_global_definitions.i"

%import "iact_data/telescope_event.pb.i"
%import "iact_data/nectarcam_configuration.pb.i"
%import "iact_data/lstcam_configuration.pb.i"
%import "iact_data/telescope_run_configuration.pb.i"
%import "iact_data/zfits_data_source.pb.i"
%import "iact_data/nectarcam_data_source.pb.i"
%import "iact_data/lstcam_data_source.pb.i"
%import "iact_data/cta_data_source.pb.i"

%newobject decode_nmc_xml_file(const std::string& filename);
%include "iact_data/nectarcam_configuration.hpp"

%newobject get_run_configuration();

%ignore get_next(uint64_t& seq_index_out);
%ignore get_next(uint64_t& seq_index_out, google::protobuf::Arena** arena);

%typemap(in, numinputs=0) calin::ix::iact_data::telescope_event::TelescopeEvent** CALIN_PROTOBUF_OUTPUT
  (calin::ix::iact_data::telescope_event::TelescopeEvent* temp = nullptr) {
    // typemap(in) calin::ix::iact_data::telescope_event::TelescopeEvent** CALIN_PROTOBUF_OUTPUT - telescope_data_source.i
    $1 = &temp;
}

%typemap(argout)
  calin::ix::iact_data::telescope_event::TelescopeEvent** CALIN_PROTOBUF_OUTPUT {
    // typemap(argout) calin::ix::iact_data::telescope_event::TelescopeEvent** CALIN_PROTOBUF_OUTPUT - telescope_data_source.i
    %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(*$1), $*1_descriptor, 0));
}

%apply calin::ix::iact_data::telescope_event::TelescopeEvent** CALIN_PROTOBUF_OUTPUT {
  calin::ix::iact_data::telescope_event::TelescopeEvent** event_out };
%apply uint64_t& OUTPUT { uint64_t& seq_index_out };
%apply google::protobuf::Arena** CALIN_ARENA_OUTPUT {
  google::protobuf::Arena** arena_out };

%newobject new_data_source();
%newobject new_data_sink();

%import "io/data_source.i"

/* %import "io/data_source.hpp"
%import "io/chained_data_source.hpp" */

%template(TelescopeDataSource)
  calin::io::data_source::DataSource<
    calin::ix::iact_data::telescope_event::TelescopeEvent>;

%newobject simple_get_next();

%extend calin::io::data_source::DataSource<
    calin::ix::iact_data::telescope_event::TelescopeEvent> {

  calin::ix::iact_data::telescope_event::TelescopeEvent* simple_get_next()
  {
    uint64_t unused_seq_index = 0;
    return $self->get_next(unused_seq_index);
  }

  void get_next(uint64_t& seq_index_out,
    calin::ix::iact_data::telescope_event::TelescopeEvent** event_out,
    google::protobuf::Arena** arena_out)
  {
    seq_index_out = 0;
    *event_out = $self->get_next(seq_index_out, arena_out);
    if(*event_out != nullptr and *arena_out == nullptr)
      throw std::runtime_error("Memory allocation error: no Arena returned");
  }
}

%template(TelescopeRandomAccessDataSource)
  calin::io::data_source::RandomAccessDataSource<
    calin::ix::iact_data::telescope_event::TelescopeEvent>;

%template(VectorTelescopeDataSource)
  std::vector<calin::io::data_source::DataSource<
    calin::ix::iact_data::telescope_event::TelescopeEvent>*>;

%template(VectorTelescopeRandomAccessDataSource)
  std::vector<calin::io::data_source::RandomAccessDataSource<
    calin::ix::iact_data::telescope_event::TelescopeEvent>*>;

%template(FileTelescopeDataSource)
  calin::io::data_source::ProtobufFileDataSource<
    calin::ix::iact_data::telescope_event::TelescopeEvent>;

%template(TelescopeDataSink)
  calin::io::data_source::DataSink<
    calin::ix::iact_data::telescope_event::TelescopeEvent>;

%template(FileTelescopeDataSink)
  calin::io::data_source::ProtobufFileDataSink<
    calin::ix::iact_data::telescope_event::TelescopeEvent>;

%template(TelescopeDataSourceFactory)
  calin::io::data_source::DataSourceFactory<
    calin::ix::iact_data::telescope_event::TelescopeEvent>;

%template(TelescopeDataSinkFactory)
  calin::io::data_source::DataSinkFactory<
    calin::ix::iact_data::telescope_event::TelescopeEvent>;

%include "io/buffered_data_source.hpp"


%template(BufferedTelescopeDataSource)
  calin::io::data_source::BufferedDataSource<
    calin::ix::iact_data::telescope_event::TelescopeEvent>;

%template(MultiThreadTelescopeDataSourceBuffer)
  calin::io::data_source::UnidirectionalBufferedDataSourcePump<
    calin::ix::iact_data::telescope_event::TelescopeEvent>;

%include "iact_data/telescope_data_source.hpp"

%template(DataSourceOpenerTelescopeRandomAccessDataSource)
  calin::io::data_source::DataSourceOpener<
    calin::io::data_source::RandomAccessDataSource<
      calin::ix::iact_data::telescope_event::TelescopeEvent> >;
//    calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSource>;

%template(BasicChainedDataSourceTelescopeRandomAccessDataSource)
  calin::io::data_source::BasicChainedDataSource<
    calin::io::data_source::RandomAccessDataSource<
      calin::ix::iact_data::telescope_event::TelescopeEvent> >;

%template(BasicChainedRandomAccessDataSourceTelescopeRandomAccessDataSource)
  calin::io::data_source::BasicChainedRandomAccessDataSource<
    calin::io::data_source::RandomAccessDataSource<
      calin::ix::iact_data::telescope_event::TelescopeEvent> >;

%template(DataSourceOpenerTelescopeRandomAccessDataSourceWithRunConfig)
  calin::io::data_source::DataSourceOpener<
    calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig>;

%template(BasicChainedDataSourceTelescopeRandomAccessDataSourceWithRunConfig)
  calin::io::data_source::BasicChainedDataSource<
    calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig>;

%template(BasicChainedRandomAccessDataSourceTelescopeRandomAccessDataSourceWithRunConfig)
  calin::io::data_source::BasicChainedRandomAccessDataSource<
    calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig>;

%template(VectorTelescopeRandomAccessDataSourceWithRunConfig)
  std::vector<calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig*>;

%include "pattern/delegation.hpp"

%template(DelegatorTelescopeRandomAccessDataSourceWithRunConfig)
  calin::pattern::delegation::Delegator<
    calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig>;

%include "iact_data/zfits_data_source.hpp"

%include "iact_data/nectarcam_data_source.hpp"
%include "iact_data/lstcam_data_source.hpp"
%include "iact_data/cta_data_source.hpp"
