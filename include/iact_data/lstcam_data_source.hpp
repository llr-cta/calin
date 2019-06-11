/*

   calin/iact_data/lstcam_data_source.hpp -- Stephen Fegan -- 2018-10-16

   A supplier of single telescope data from LSTCam DAQ data files

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
#include <iact_data/lstcam_actl_event_decoder.hpp>
#include <iact_data/lstcam_data_source.pb.h>
#include <iact_data/zfits_data_source.hpp>
#include <pattern/delegation.hpp>

namespace calin { namespace iact_data { namespace lstcam_data_source {

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

class LSTCamZFITSDataSource_R1:
  public calin::iact_data::zfits_data_source::ZFITSDataSource_R1
{
public:
  CALIN_TYPEALIAS(decoder_config_type,
    calin::ix::iact_data::lstcam_data_source::LSTCamCameraEventDecoderConfig);

#if 0
  void set_decoder_config(const decoder_config_type& config) {
    decoder_->set_config(config); }
#endif
  calin::ix::iact_data::lstcam_data_source::LSTCamCameraEventDecoderConfig decoder_config() const {
    return decoder_->config(); }
#if 0
  decoder_config_type* mutable_decoder_config() {
    return decoder_->mutable_config(); }
#endif
  static calin::ix::iact_data::lstcam_data_source::LSTCamCameraEventDecoderConfig default_decoder_config() {
    return lstcam_actl_event_decoder::LSTCam_ACTL_R1_CameraEventDecoder::default_config(); }

  LSTCamZFITSDataSource_R1(const std::string& filename,
    const config_type& config,
    const decoder_config_type& decoder_config = default_decoder_config());
  LSTCamZFITSDataSource_R1(const std::string& filename,
    const decoder_config_type& decoder_config = default_decoder_config(),
    const config_type& config = default_config());
  virtual ~LSTCamZFITSDataSource_R1();
private:
  lstcam_actl_event_decoder::LSTCam_ACTL_R1_CameraEventDecoder* decoder_;
};

#endif

} } } // namespace calin::iact_data::lstcam_data_source
