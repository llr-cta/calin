/*

   calin/iact_data/lstcam_data_source.hpp -- Stephen Fegan -- 2018-10-16

   A supplier of single telescope data from LSTCam DAQ data files

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
#include <iact_data/lstcam_acada_event_decoder.hpp>
#include <iact_data/lstcam_data_source.pb.h>
#include <iact_data/zfits_data_source.hpp>
#include <pattern/delegation.hpp>

namespace calin { namespace iact_data { namespace lstcam_data_source {

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

class LSTCamZFITSDataSource_R1v0:
  public calin::iact_data::zfits_data_source::ZFITSDataSource_R1v0
{
public:
  CALIN_TYPEALIAS(decoder_config_type,
    calin::ix::iact_data::lstcam_data_source::LSTCamCameraEventDecoderConfig);

  calin::ix::iact_data::lstcam_data_source::LSTCamCameraEventDecoderConfig decoder_config() const {
    return decoder_->config(); }

  static calin::ix::iact_data::lstcam_data_source::LSTCamCameraEventDecoderConfig default_decoder_config() {
    return lstcam_acada_event_decoder::LSTCam_ACADACameraEventDecoder_R1v0::default_config(); }

  LSTCamZFITSDataSource_R1v0(const std::string& filename,
    calin::iact_data::zfits_acada_data_source::
      ZFITSACADACameraEventDataSource<calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0>* acada_zfits,
    const decoder_config_type& decoder_config = default_decoder_config(),
    bool adopt_acada_zfits = false);
  LSTCamZFITSDataSource_R1v0(const std::string& filename,
    const config_type& config,
    const decoder_config_type& decoder_config = default_decoder_config());
  LSTCamZFITSDataSource_R1v0(const std::string& filename,
    const decoder_config_type& decoder_config = default_decoder_config(),
    const config_type& config = default_config());
  virtual ~LSTCamZFITSDataSource_R1v0();
private:
  lstcam_acada_event_decoder::LSTCam_ACADACameraEventDecoder_R1v0* decoder_;
};

} } } // namespace calin::iact_data::lstcam_data_source
