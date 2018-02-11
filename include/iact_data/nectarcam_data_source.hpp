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

namespace calin { namespace iact_data { namespace nectarcam_data_source {

#ifdef CALIN_HAVE_CTA_CAMERASTOACTL

class NectarCamCameraEventDecoder:
  public zfits_data_source::CTACameraEventDecoder
{
public:
  CALIN_TYPEALIAS(config_type, calin::ix::iact_data::
    nectarcam_data_source::NectarCamCameraEventDecoderConfig);

  NectarCamCameraEventDecoder(const std::string& filename,
    unsigned run_number = 0,
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

  virtual ~NectarCamCameraEventDecoder();

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

class NectarCamZFITSDataSource:
  public calin::iact_data::zfits_data_source::ZFITSDataSource
{
public:
  CALIN_TYPEALIAS(decoder_config_type,
    NectarCamCameraEventDecoder::config_type);

#if 0
  void set_decoder_config(const decoder_config_type& config) {
    decoder_->set_config(config); }
#endif
  decoder_config_type decoder_config() const {
    return decoder_->config(); }
#if 0
  decoder_config_type* mutable_decoder_config() {
    return decoder_->mutable_config(); }
#endif
  static decoder_config_type default_decoder_config() {
    return NectarCamCameraEventDecoder::default_config(); }

  NectarCamZFITSDataSource(const std::string& filename,
    const config_type& config,
    const decoder_config_type& decoder_config = default_decoder_config());
  NectarCamZFITSDataSource(const std::string& filename,
    const decoder_config_type& decoder_config = default_decoder_config(),
    const config_type& config = default_config());
  virtual ~NectarCamZFITSDataSource();
private:
  NectarCamCameraEventDecoder* decoder_;
};

#endif

} } } // namespace calin::iact_data::nectarcam_data_source
