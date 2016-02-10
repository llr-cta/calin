/*

   calin/iact_data/telescope_data_source.i -- Stephen Fegan -- 2016-01-14

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

%module (package="calin.iact_data") telescope_data_source

%{
#include "iact_data/telescope_data_source.hpp"
#include "iact_data/nectarcam_data_source.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "numpy.i"
%include "stdint.i"
%include "calin_typemaps.i"

%import "calin_global_definitions.i"
%import "iact_data/telescope_event.pb.i"
%import "iact_data/zfits_data_source.pb.i"
%import "iact_data/nectarcam_data_source.pb.i"

%newobject get_next();
//%newobject calin::iact_data::telescope_data_source::TelescopeDataSource::get_next();
//%newobject calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSource::get_next();
//%newobject calin::iact_data::zfits_data_source::ZFITSSingleFileDataSource::get_next();
//%newobject calin::iact_data::nectarcam_data_source::NectarCamZFITSDataSource::get_next();

%import "io/data_source.hpp"

%template(TelescopeDataSource)
  calin::io::data_source::DataSource<
    calin::ix::iact_data::telescope_event::TelescopeEvent>;

%template(TelescopeRandomAccessDataSource)
  calin::io::data_source::RandomAccessDataSource<
    calin::ix::iact_data::telescope_event::TelescopeEvent>;

%template(FileTelescopeDataSource)
  calin::io::data_source::ProtobufFileDataSource<
    calin::ix::iact_data::telescope_event::TelescopeEvent>;

%template(TelescopeDataSink)
  calin::io::data_source::DataSink<
    calin::ix::iact_data::telescope_event::TelescopeEvent>;

%template(FileTelescopeDataSink)
  calin::io::data_source::ProtobufFileDataSink<
    calin::ix::iact_data::telescope_event::TelescopeEvent>;

%include "iact_data/telescope_data_source.hpp"

%template(DataSourceOpenerTelescopeRandomAccessDataSource)
  calin::io::data_source::DataSourceOpener<
    calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSource>;

%template(BasicChainedDataSourceTelescopeRandomAccessDataSource)
  calin::io::data_source::BasicChainedDataSource<
    calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSource>;

%template(BasicChaninedRandomAccessDataSourceTelescopeRandomAccessDataSource)
  calin::io::data_source::BasicChaninedRandomAccessDataSource<
    calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSource>;

%include "iact_data/zfits_data_source.hpp"
%include "iact_data/nectarcam_data_source.hpp"
