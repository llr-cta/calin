/*

   calin/io/nectarcam_data_source.hpp -- Stephen Fegan -- 2016-01-11

   A supplier of single telescope data from NectarCam DAQ data files

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

#include <string>

#include <calin_global_config.hpp>
#include <io/telescope_data_source.hpp>

namespace calin { namespace io { namespace nectarcam_data_source {

#ifdef CALIN_HAVE_CTA_CAMERASTOACTL

class NectarCamZFITSDataSource:
  public calin::io::telescope_data_source::TelescopeDataSource
{
public:
  NectarCamZFITSDataSource(const std::string& filename);
  virtual ~NectarCamZFITSDataSource();
  calin::ix::iact::telescope_event::TelescopeEvent* getNextEvent() override;
private:
  std::string filename_;
  unsigned next_event_index_ = 0;
};

#endif // #ifdef CALIN_HAVE_CTA_CAMERASTOACTL

} } } // namespace calin::io::nectarcam_data_source
