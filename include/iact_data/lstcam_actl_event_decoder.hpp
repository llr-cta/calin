/*

   calin/iact_data/lstcam_data_source.hpp -- Stephen Fegan -- 2018-10-15

   A supplier of single telescope data from LSTCam DAQ data files

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

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
#include <iact_data/lstcam_data_source.pb.h>
#include <iact_data/zfits_data_source.hpp>

namespace calin { namespace iact_data { namespace lstcam_actl_event_decoder {

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

class LSTCam_ACTL_R1_CameraEventDecoder:
  public actl_event_decoder::ACTL_R1_CameraEventDecoder
{
public:
  CALIN_TYPEALIAS(config_type, calin::ix::iact_data::
    lstcam_data_source::LSTCamCameraEventDecoderConfig);

  LSTCam_ACTL_R1_CameraEventDecoder(const std::string& filename, unsigned run_number = 0,
    const calin::ix::iact_data::lstcam_data_source::LSTCamCameraEventDecoderConfig& config = default_config());

  ~LSTCam_ACTL_R1_CameraEventDecoder();

  virtual bool decode(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    const R1::CameraEvent* cta_event) override;

  virtual bool decode_run_config(
    calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config,
    const R1::CameraConfiguration* cta_run_header,
    const R1::CameraEvent* cta_event) override;

  calin::ix::iact_data::lstcam_data_source::LSTCamCameraEventDecoderConfig config() const { return config_; }

  static calin::ix::iact_data::lstcam_data_source::LSTCamCameraEventDecoderConfig default_config() {
    config_type config = config_type::default_instance();
    config.set_separate_channel_waveforms(true);
    config.set_counts_to_time_133megahertz(30797ULL);
    return config;
  }

protected:
  void copy_single_gain_waveforms(
    calin::ix::iact_data::telescope_event::DigitizedSkyImage* calin_image,
    const int16_t* cta_waveforms, const uint8_t* cta_pixel_mask,
    uint8 has_gain_mask, const std::string& which_gain) const;

  config_type config_;
  std::string filename_;
  unsigned nmod_ = 0;
  unsigned nsample_ = 0;
  unsigned run_number_ = 0;
  unsigned telescope_id_ = 0;
  bool exchange_gain_channels_ = false;
  int64_t run_start_time_ = 0;
};

#endif

} } } // namespace calin::iact_data::lstcam_data_source
