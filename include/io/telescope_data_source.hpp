/*

   calin/io/telescope_data_source.hpp -- Stephen Fegan -- 2016-01-08

   A supplier of single telescope data, for example:
   ix::iact::telescope_event::TelescopeEvent

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

#pragma once

#include <io/data_source.hpp>
#include <iact/telescope_event.pb.h>

namespace calin { namespace io { namespace telescope_data_source {

class TelescopeDataSource:
  public calin::io::data_source::DataSource<
    calin::ix::iact::telescope_event::TelescopeEvent>
{
public:
  virtual ~TelescopeDataSource();
  virtual calin::ix::iact::telescope_event::TelescopeEvent* getNext() = 0;
};

class RawFileTelescopeDataSource final public TelescopeDataSource 
{
public:
  RawFileTelescopeDataSource(const std::string& filename);
  ~RawFileTelescopeDataSource();
  calin::ix::iact::telescope_event::TelescopeEvent* getNext() override;
};

} } } // namespace calin::io::telescope_data_source
