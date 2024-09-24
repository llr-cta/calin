/*

   calin/iact_data/nectarcam_acada_event_decoder.hpp -- Stephen Fegan -- 2016-01-11

   A supplier of single telescope data from NectarCam DAQ data files

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

#pragma once

#include <string>

#include <calin_global_definitions.hpp>
#include <calin_global_config.hpp>
#include <iact_data/nectarcam_data_source.pb.h>
#include <iact_data/acada_data_source.hpp>
#include <iact_data/acada_event_decoder.hpp>
#include <iact_data/unified_acada_event_decoder.hpp>

namespace calin { namespace iact_data { namespace nectarcam_acada_event_decoder {

/*

              LLLLLLLLLLL                       000000000
              L:::::::::L                     00:::::::::00
              L:::::::::L                   00:::::::::::::00
              LL:::::::LL                  0:::::::000:::::::0
                L:::::L                    0::::::0   0::::::0
                L:::::L                    0:::::0     0:::::0
                L:::::L                    0:::::0     0:::::0
                L:::::L                    0:::::0 000 0:::::0
                L:::::L                    0:::::0 000 0:::::0
                L:::::L                    0:::::0     0:::::0
                L:::::L                    0:::::0     0:::::0
                L:::::L         LLLLLL     0::::::0   0::::::0
              LL:::::::LLLLLLLLL:::::L     0:::::::000:::::::0
              L::::::::::::::::::::::L      00:::::::::::::00
              L::::::::::::::::::::::L        00:::::::::00
              LLLLLLLLLLLLLLLLLLLLLLLL          000000000

*/

class NectarCam_ACADACameraEventDecoder_L0:
  public calin::iact_data::acada_event_decoder::ACADACameraEventDecoder<
    calin::iact_data::acada_data_source::ACADA_MessageSet_L0>
{
public:
  CALIN_TYPEALIAS(config_type, calin::ix::iact_data::
    nectarcam_data_source::NectarCamCameraEventDecoderConfig);

  CALIN_TYPEALIAS(message_set_type, calin::iact_data::acada_data_source::ACADA_MessageSet_L0);
  CALIN_TYPEALIAS(event_type, calin::iact_data::acada_data_source::ACADA_MessageSet_L0::event_type);
  CALIN_TYPEALIAS(header_type, calin::iact_data::acada_data_source::ACADA_MessageSet_L0::header_type);
  CALIN_TYPEALIAS(data_stream_type, calin::iact_data::acada_data_source::ACADA_MessageSet_L0::data_stream_type);

  NectarCam_ACADACameraEventDecoder_L0(const std::string& filename, unsigned run_number = 0,
    const calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig& config = default_config());

  calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig config() const { return config_; }
  static calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig default_config() {
    config_type config = config_type::default_instance();
    config.set_nmc_xml_suffix(".NMC.xml");
    config.set_separate_channel_waveforms(true);
    return config;
  }

  virtual ~NectarCam_ACADACameraEventDecoder_L0();

  bool decode(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    const calin::iact_data::acada_data_source::ACADA_MessageSet_L0& cta_messages) override;

  bool decode_run_config(
    calin::ix::iact_data::telescope_run_configuration:: TelescopeRunConfiguration* run_config,
    const calin::iact_data::acada_data_source::ACADA_MessageSet_L0& cta_messages) override;

  NectarCam_ACADACameraEventDecoder_L0* clone() const override;

protected:
  void copy_single_gain_integrals(const event_type* cta_event,
    const calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
    const ProtoDataModel::PixelsChannel& cta_image,
    calin::ix::iact_data::telescope_event::DigitizedSkyImage* calin_image,
    const std::string& which_gain,
    calin::ix::iact_data::telescope_event::SignalType signal_type) const;

  virtual void copy_single_gain_waveforms(const event_type* cta_event,
    const calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
    const ProtoDataModel::PixelsChannel& cta_image,
    calin::ix::iact_data::telescope_event::DigitizedSkyImage* calin_image,
    const std::string& which_gain,
    calin::ix::iact_data::telescope_event::SignalType signal_type) const;

  unsigned get_nmod_from_event(const ProtoDataModel::CameraEvent* cta_event) const;

  config_type config_;
  std::string filename_;
  unsigned run_number_ = 0;
  bool exchange_gain_channels_ = false;
  int64_t run_start_time_ = 0;
};

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

class NectarCam_ACADACameraEventDecoder_R1v0:
  public calin::iact_data::acada_event_decoder::ACADACameraEventDecoder<
    calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0>
{
public:
  CALIN_TYPEALIAS(config_type, calin::ix::iact_data::
    nectarcam_data_source::NectarCamCameraEventDecoderConfig);

  CALIN_TYPEALIAS(message_set_type, calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0);
  CALIN_TYPEALIAS(event_type, calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0::event_type);
  CALIN_TYPEALIAS(header_type, calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0::header_type);
  CALIN_TYPEALIAS(data_stream_type, calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0::data_stream_type);

  NectarCam_ACADACameraEventDecoder_R1v0(const std::string& filename, unsigned run_number = 0,
    const config_type& config = default_config());

  ~NectarCam_ACADACameraEventDecoder_R1v0();

  virtual bool decode(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    const calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0& cta_messages) override;

  bool decode_run_config(
    calin::ix::iact_data::telescope_run_configuration:: TelescopeRunConfiguration* run_config,
    const calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0& cta_messages) override;

  NectarCam_ACADACameraEventDecoder_R1v0* clone() const override;

  calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig config() const { return config_; }

  static calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig default_config() {
    config_type config = config_type::default_instance();
    config.set_nmc_xml_suffix(".NMC.xml");
    config.set_separate_channel_waveforms(true);
    return config;
  }

protected:
  void copy_single_gain_integrals(const event_type* cta_event,
    const calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
    const int16_t* cta_charges, const uint8_t* cta_pixel_mask,
    calin::ix::iact_data::telescope_event::DigitizedSkyImage* calin_image,
    const std::string& which_gain) const;

  void copy_single_gain_waveforms(
    const calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
    const int16_t* cta_waveforms, const uint8_t* cta_pixel_mask,
    calin::ix::iact_data::telescope_event::Waveforms* calin_waveforms,
    uint8_t has_gain_mask, const std::string& which_gain) const;

  config_type config_;
  std::string filename_;
  unsigned nmod_ = 0;
  unsigned nsample_ = 0;
  unsigned run_number_ = 0;
  unsigned telescope_id_ = 0;
  bool exchange_gain_channels_ = false;
  int64_t run_start_time_ = 0;
};

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

class NectarCam_ACADACameraEventDecoder_R1v1:
  public calin::iact_data::unified_acada_event_decoder::Unified_ACADACameraEventDecoder_R1v1
{
public:
  CALIN_TYPEALIAS(config_type, calin::ix::iact_data::
    nectarcam_data_source::NectarCamCameraEventDecoderConfig);

  CALIN_TYPEALIAS(message_set_type, calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1);
  CALIN_TYPEALIAS(event_type, calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1::event_type);
  CALIN_TYPEALIAS(header_type, calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1::header_type);
  CALIN_TYPEALIAS(data_stream_type, calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1::data_stream_type);

  NectarCam_ACADACameraEventDecoder_R1v1(const std::string& filename, 
    const config_type& config = default_config());

  ~NectarCam_ACADACameraEventDecoder_R1v1();

  virtual bool decode(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    const calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1& cta_messages) override;

  bool decode_run_config(
    calin::ix::iact_data::telescope_run_configuration:: TelescopeRunConfiguration* run_config,
    const calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1& cta_messages) override;

  NectarCam_ACADACameraEventDecoder_R1v1* clone() const override;

  calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig config() const;

  static calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig default_config() {
    config_type config = config_type::default_instance();
    // config.set_nmc_xml_suffix(".NMC.xml");
    config.set_separate_channel_waveforms(true);
    return config;
  }

protected:
  static inline calin::ix::iact_data::cta_data_source::UnifiedCameraEventDecoderConfig unified_decoder_config(
      config_type config);

};


} } } // namespace calin::iact_data::nectarcam_acada_event_decoder
