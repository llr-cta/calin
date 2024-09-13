/*

   calin/iact_data/zfits_acada_data_source.cpp -- Stephen Fegan -- 2024-09-05

   Base class for source of "raw" ACADA data types

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

#include <iact_data/acada_data_source.hpp>

using namespace calin::iact_data::acada_data_source;

template<typename EventMessage>
ACADACameraEventDataSource<EventMessage>::
~ACADACameraEventDataSource()
{
  // nothing to see here
}

template<typename EventMessage, typename HeaderMessage, typename DataStreamMessage>
ACADACameraEventDataSourceWithRunHeader<EventMessage,HeaderMessage,DataStreamMessage>::
~ACADACameraEventDataSourceWithRunHeader()
{
  // nothing to see here
}

template<typename EventMessage, typename HeaderMessage, typename DataStreamMessage>
ACADACameraEventRandomAccessDataSourceWithRunHeader<EventMessage,HeaderMessage,DataStreamMessage>::
~ACADACameraEventRandomAccessDataSourceWithRunHeader()
{
  // nothing to see here
}

namespace calin { namespace iact_data { namespace acada_data_source {

template class ACADACameraEventDataSource<ACADA_EventMessage_L0>;
template class ACADACameraEventDataSourceWithRunHeader<ACADA_EventMessage_L0, ACADA_HeaderMessage_L0>;
template class ACADACameraEventRandomAccessDataSourceWithRunHeader<ACADA_EventMessage_L0, ACADA_HeaderMessage_L0>;

template class ACADACameraEventDataSource<ACADA_EventMessage_R1v0>;
template class ACADACameraEventDataSourceWithRunHeader<ACADA_EventMessage_R1v0, ACADA_HeaderMessage_R1v0>;
template class ACADACameraEventRandomAccessDataSourceWithRunHeader<ACADA_EventMessage_R1v0, ACADA_HeaderMessage_R1v0>;

template class ACADACameraEventDataSource<ACADA_EventMessage_R1v1>;
template class ACADACameraEventDataSourceWithRunHeader<
  ACADA_EventMessage_R1v1, ACADA_HeaderMessage_R1v1, ACADA_DataStreamMessage_R1v1>;
template class ACADACameraEventRandomAccessDataSourceWithRunHeader<
  ACADA_EventMessage_R1v1, ACADA_HeaderMessage_R1v1, ACADA_DataStreamMessage_R1v1>;

} } } // namespace calin::iact_data::acada_data_source
