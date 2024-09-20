/*

   calin/iact_data/unified_acada_event_decoder.hpp -- Stephen Fegan -- 2024-09-17

   Decoder of ACADA event types unified across all telescope types

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

#pragma once

#include <string>

#include <calin_global_definitions.hpp>
#include <calin_global_config.hpp>
#include <iact_data/acada_event_decoder.hpp>
#include <iact_data/cta_data_source.pb.h>

namespace calin { namespace iact_data { namespace unified_acada_event_decoder {

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

class Unified_ACADACameraEventDecoder_R1v1:
  public calin::iact_data::acada_event_decoder::ACADACameraEventDecoder<
    calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1>
{
public:
  CALIN_TYPEALIAS(config_type, calin::ix::iact_data::cta_data_source::UnifiedCameraEventDecoderConfig);

  CALIN_TYPEALIAS(message_set_type, calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1);
  CALIN_TYPEALIAS(event_type, calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1::event_type);
  CALIN_TYPEALIAS(header_type, calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1::header_type);
  CALIN_TYPEALIAS(data_stream_type, calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1::data_stream_type);

  Unified_ACADACameraEventDecoder_R1v1(const std::string& filename,
    const config_type& config = default_config());

  ~Unified_ACADACameraEventDecoder_R1v1();

  bool decode(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    const calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1& cta_messages) override;

  bool decode_run_config(
    calin::ix::iact_data::telescope_run_configuration:: TelescopeRunConfiguration* run_config,
    const calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1& cta_messages) override;

  Unified_ACADACameraEventDecoder_R1v1* clone() const override;

  calin::ix::iact_data::cta_data_source::UnifiedCameraEventDecoderConfig config() const { return config_; }

  static calin::ix::iact_data::cta_data_source::UnifiedCameraEventDecoderConfig default_config() {
    config_type config;
    return config;
  }

protected:
  void copy_single_gain_waveforms(
    const calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
    const int16_t* cta_waveforms, const uint8_t* cta_pixel_mask, unsigned nsample,
    calin::ix::iact_data::telescope_event::Waveforms* calin_waveforms,
    uint8_t has_gain_mask, const std::string& which_gain) const;


  config_type config_;
  std::string filename_;
  int64_t run_start_time_ = 0;
  int32_t nmod_configured_ = 0;
  int32_t nchan_configured_ = 0;
  int32_t ncamera_clock_ = 0;
  int32_t nmodule_clock_ = 0;
};

} } } // namespace calin::iact_data::unified_event_decoder
