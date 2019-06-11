/*

   calin/iact_data/cta_data_source.hpp -- Stephen Fegan -- 2018-11-22

   A supplier of single telescope data from CTA DAQ data files

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

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
#include <iact_data/nectarcam_actl_event_decoder.hpp>
#include <iact_data/lstcam_actl_event_decoder.hpp>
#include <pattern/delegation.hpp>

namespace calin { namespace iact_data { namespace cta_actl_event_decoder {

#ifdef CALIN_HAVE_CTA_CAMERASTOACTL

/*

                      RRRRRRRRRRRRRRRRR          1111111
                      R::::::::::::::::R        1::::::1
                      R::::::RRRRRR:::::R      1:::::::1
                      RR:::::R     R:::::R     111:::::1
                        R::::R     R:::::R        1::::1
                        R::::R     R:::::R        1::::1
                        R::::RRRRRR:::::R         1::::1
                        R:::::::::::::RR          1::::l
                        R::::RRRRRR:::::R         1::::l
                        R::::R     R:::::R        1::::l
                        R::::R     R:::::R        1::::l
                        R::::R     R:::::R        1::::l
                      RR:::::R     R:::::R     111::::::111
                      R::::::R     R:::::R     1::::::::::1
                      R::::::R     R:::::R     1::::::::::1
                      RRRRRRRR     RRRRRRR     111111111111

*/


class CTA_ACTL_R1_CameraEventDecoder:
  public calin::iact_data::actl_event_decoder::ACTL_R1_CameraEventDecoder,
  public calin::pattern::delegation::Delegator<
    calin::iact_data::actl_event_decoder::ACTL_R1_CameraEventDecoder>
{
public:
  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig);

  CTA_ACTL_R1_CameraEventDecoder(const std::string& filename, unsigned run_number = 0,
    const calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig& config = default_config());

  ~CTA_ACTL_R1_CameraEventDecoder();

  virtual bool decode(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    const R1::CameraEvent* cta_event) override;

  virtual bool decode_run_config(
    calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config,
    const R1::CameraConfiguration* cta_run_header,
    const R1::CameraEvent* cta_event) override;

  calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig config() const { return config_; }

  static calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig default_config();

protected:
  std::string filename_;
  unsigned run_number_ = 0;
  calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig config_ = default_config();
};

#endif

} } } // namespace calin::iact_data::cta_actl_event_decoder
