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

template<typename EventMessage, typename HeaderMessage>
ACADACameraEventRandomAccessDataSourceWithRunHeader<EventMessage,HeaderMessage>::
ACADACameraEventRandomAccessDataSourceWithRunHeader():
  calin::io::data_source::RandomAccessDataSource<EventMessage>()
{
  // nothing to see here
}

template<typename EventMessage, typename HeaderMessage>
ACADACameraEventRandomAccessDataSourceWithRunHeader<EventMessage,HeaderMessage>::
~ACADACameraEventRandomAccessDataSourceWithRunHeader()
{
  // nothing to see here
}

namespace calin { namespace iact_data { namespace acada_data_source {

template class ACADACameraEventRandomAccessDataSourceWithRunHeader<ACADA_L0_EventMessage, ACADA_L0_HeaderMessage>;
template class ACADACameraEventRandomAccessDataSourceWithRunHeader<ACADA_R1v0_EventMessage, ACADA_R1v0_HeaderMessage>;

} } } // namespace calin::iact_data::acada_data_source
