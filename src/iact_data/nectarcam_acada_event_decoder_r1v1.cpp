/*

   calin/iact_data/nectarcam_acada_event_decoder_r1v0.cpp -- Stephen Fegan -- 2018-09-21

   A decoder of NectarCAM ACADA data in prototype R1 format

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <stdexcept>
#include <string>
#include <memory>
#include <numeric>

#include <util/log.hpp>
#include <util/file.hpp>
#include <iact_data/acada_data_source.hpp>
#include <iact_data/nectarcam_acada_event_decoder.hpp>
#include <iact_data/nectarcam_layout.hpp>
#include <iact_data/nectarcam_configuration.hpp>
#include <provenance/system_info.hpp>

using namespace calin::iact_data::acada_data_source;
using namespace calin::iact_data::nectarcam_acada_event_decoder;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::util::log;

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

NectarCam_ACADACameraEventDecoder_R1v1::
NectarCam_ACADACameraEventDecoder_R1v1(const std::string& filename, 
    const config_type& config):
  calin::iact_data::unified_acada_event_decoder::Unified_ACADACameraEventDecoder_R1v1(filename, config)
{
  // nothing to see here
}

NectarCam_ACADACameraEventDecoder_R1v1::
~NectarCam_ACADACameraEventDecoder_R1v1()
{
  // nothing to see here
}

bool NectarCam_ACADACameraEventDecoder_R1v1::decode(
  calin::ix::iact_data::telescope_event::TelescopeEvent* event,
  const calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1& cta_messages)
{
  if(not Unified_ACADACameraEventDecoder_R1v1::decode(event, cta_messages)) {
    return false;
  }

  // do more work

  return true;
}    

bool NectarCam_ACADACameraEventDecoder_R1v1::decode_run_config(
  calin::ix::iact_data::telescope_run_configuration:: TelescopeRunConfiguration* run_config,
  const calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1& cta_messages)
{
  if(not Unified_ACADACameraEventDecoder_R1v1::decode_run_config(run_config, cta_messages)) {
    return false;
  }

  // do more work

  return true;
}

NectarCam_ACADACameraEventDecoder_R1v1* NectarCam_ACADACameraEventDecoder_R1v1::clone() const
{
  return new NectarCam_ACADACameraEventDecoder_R1v1(*this);
}
