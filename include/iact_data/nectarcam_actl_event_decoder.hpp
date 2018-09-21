/*

   calin/iact_data/nectarcam_data_source.hpp -- Stephen Fegan -- 2016-01-11

   A supplier of single telescope data from NectarCam DAQ data files

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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
#include <iact_data/nectarcam_data_source.pb.h>
#include <iact_data/zfits_data_source.hpp>

namespace calin { namespace iact_data { namespace nectarcam_actl_event_decoder {

#ifdef CALIN_HAVE_CTA_CAMERASTOACTL

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

class NectarCAM_ACTL_L0_CameraEventDecoder:
  public actl_event_decoder::ACTL_L0_CameraEventDecoder
{
public:
  CALIN_TYPEALIAS(config_type, calin::ix::iact_data::
    nectarcam_data_source::NectarCamCameraEventDecoderConfig);

  NectarCAM_ACTL_L0_CameraEventDecoder(const std::string& filename, unsigned run_number = 0,
    const config_type& config = default_config());

  //void set_config(const config_type& config) { config_.CopyFrom(config); }
  config_type config() const { return config_; }
  //config_type* mutable_config() { return &config_; }
  static config_type default_config() {
    config_type config = config_type::default_instance();
    config.set_nmc_xml_suffix(".NMC.xml");
    config.set_separate_channel_waveforms(true);
    return config;
  }

  virtual ~NectarCAM_ACTL_L0_CameraEventDecoder();

  bool decode(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    const DataModel::CameraEvent* cta_event) override;

  bool decode_run_config(
    calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config,
    const DataModel::CameraRunHeader* cta_run_header,
    const DataModel::CameraEvent* cta_event) override;

protected:
  void copy_single_gain_integrals(const DataModel::CameraEvent* cta_event,
    const calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
    const DataModel::PixelsChannel& cta_image,
    calin::ix::iact_data::telescope_event::DigitizedSkyImage* calin_image,
    const std::string& which_gain) const;

  virtual void copy_single_gain_waveforms(const DataModel::CameraEvent* cta_event,
    const calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
    const DataModel::PixelsChannel& cta_image,
    calin::ix::iact_data::telescope_event::DigitizedSkyImage* calin_image,
    const std::string& which_gain) const;

  unsigned get_nmod_from_event(const DataModel::CameraEvent* cta_event) const;

  config_type config_;
  std::string filename_;
  unsigned run_number_ = 0;
  bool exchange_gain_channels_ = false;
  int64_t run_start_time_ = 0;
};

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

class NectarCAM_ACTL_R1_CameraEventDecoder:
  public actl_event_decoder::ACTL_R1_CameraEventDecoder
{
public:
  CALIN_TYPEALIAS(config_type, calin::ix::iact_data::
    nectarcam_data_source::NectarCamCameraEventDecoderConfig);

  NectarCAM_ACTL_R1_CameraEventDecoder(const std::string& filename, unsigned run_number = 0,
    const config_type& config = default_config());

  ~NectarCAM_ACTL_R1_CameraEventDecoder();

  virtual bool decode(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    const R1::CameraEvent* cta_event) override;

  virtual bool decode_run_config(
    calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config,
    const R1::CameraConfiguration* cta_run_header,
    const R1::CameraEvent* cta_event) override;

  config_type config() const { return config_; }

  static config_type default_config() {
    config_type config = config_type::default_instance();
    config.set_nmc_xml_suffix(".NMC.xml");
    config.set_separate_channel_waveforms(true);
    return config;
  }

protected:
  void copy_single_gain_integrals(const DataModel::CameraEvent* cta_event,
    const calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
    const DataModel::PixelsChannel& cta_image,
    calin::ix::iact_data::telescope_event::DigitizedSkyImage* calin_image,
    const std::string& which_gain) const;

  virtual void copy_single_gain_waveforms(const DataModel::CameraEvent* cta_event,
    const calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
    const DataModel::PixelsChannel& cta_image,
    calin::ix::iact_data::telescope_event::DigitizedSkyImage* calin_image,
    const std::string& which_gain) const;

  config_type config_;
  std::string filename_;
  unsigned run_number_ = 0;
  bool exchange_gain_channels_ = false;
  int64_t run_start_time_ = 0;  
};

#endif

} } } // namespace calin::iact_data::nectarcam_data_source
