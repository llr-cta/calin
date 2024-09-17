/*

   calin/iact_data/cta_acada_event_decoder.hpp -- Stephen Fegan -- 2018-11-22

   A decoder of raw ACADA data from CTA DAQ data files

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#pragma once

#include <string>

#include <calin_global_definitions.hpp>
#include <calin_global_config.hpp>
#include <iact_data/cta_data_source.pb.h>
#include <iact_data/acada_event_decoder.hpp>
#include <pattern/delegation.hpp>

namespace calin { namespace iact_data { namespace cta_acada_event_decoder {

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

class CTA_ACADACameraEventDecoder_R1v0:
  public calin::iact_data::acada_event_decoder::ACADACameraEventDecoder_R1v0,
  public calin::pattern::delegation::Delegator<
    calin::iact_data::acada_event_decoder::ACADACameraEventDecoder_R1v0>
{
public:
  CALIN_TYPEALIAS(message_set_type, calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0);
  CALIN_TYPEALIAS(event_type, calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0::event_type);
  CALIN_TYPEALIAS(header_type, calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0::header_type);
  CALIN_TYPEALIAS(data_stream_type, calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0::data_stream_type);

  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig);

  CTA_ACADACameraEventDecoder_R1v0(
    calin::iact_data::acada_event_decoder::ACADACameraEventDecoder_R1v0* decoder,
    bool adopt_decoder = false);

  CTA_ACADACameraEventDecoder_R1v0(const std::string& filename, unsigned run_number = 0,
    const calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig& config = default_config());

  ~CTA_ACADACameraEventDecoder_R1v0();

  bool decode(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    const calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0& cta_messages) override;

  bool decode_run_config(
    calin::ix::iact_data::telescope_run_configuration:: TelescopeRunConfiguration* run_config,
    const calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0& cta_messages) override;

  CTA_ACADACameraEventDecoder_R1v0  * clone() const override;

  calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig config() const { return config_; }

  static calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig default_config();

protected:
  void ensure_deligate(const calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0& cta_messages);
  std::string filename_;
  unsigned run_number_ = 0;
  calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig config_ = default_config();
};

} } } // namespace calin::iact_data::cta_acada_event_decoder
