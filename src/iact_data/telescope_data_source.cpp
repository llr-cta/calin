/*

   calin/iact_data/telescope_data_source.cpp -- Stephen Fegan -- 2016-01-08

   A supplier of single telescope data, for example:
   ix::iact_data::telescope_event::TelescopeEvent

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

//#define CALIN_TELESCOPE_DATA_SOURCE_NO_EXTERN
#include <iact_data/telescope_data_source.hpp>

using namespace calin::iact_data::telescope_data_source;

template class calin::io::data_source::DataSource<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;
template class calin::io::data_source::RandomAccessDataSource<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;
template class calin::io::data_source::ProtobufFileDataSource<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;
template class calin::io::data_source::BufferedDataSource<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;
template class calin::io::data_source::MultiThreadDataSourceBuffer<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;

template class calin::io::data_source::DataSink<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;
template class calin::io::data_source::ProtobufFileDataSink<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;

TelescopeRandomAccessDataSourceWithRunConfig::
~TelescopeRandomAccessDataSourceWithRunConfig()
{
  // nothing to see here
}
