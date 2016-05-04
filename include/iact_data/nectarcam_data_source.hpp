/*

   calin/iact_data/nectarcam_data_source.hpp -- Stephen Fegan -- 2016-01-11

   A supplier of single telescope data from NectarCam DAQ data files

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

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

  NectarCamCameraEventDecoder(const config_type& config = default_config()):
    zfits_data_source::CTACameraEventDecoder(), config_(config) { }

  void set_config(const config_type& config) { config_.CopyFrom(config); }
  const config_type& config() const { return config_; }
  config_type* mutable_config() { return &config_; }
  static config_type default_config() {
    return config_type::default_instance(); }

  virtual ~NectarCamCameraEventDecoder();
  calin::ix::iact_data::telescope_event::TelescopeEvent*
    decode(const DataModel::CameraEvent* cta_event) override;

  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* decode_run_config(
      const DataModel::CameraRunHeader* cta_run_header,
      const DataModel::CameraEvent* cta_event) override;

private:
  void copy_single_gain_image(const DataModel::CameraEvent* cta_event,
    const calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
    const DataModel::PixelsChannel& cta_image,
    calin::ix::iact_data::telescope_event::DigitizedSkyImage* calin_image,
    const std::string& which_gain) const;

  unsigned get_nmod_from_event(const DataModel::CameraEvent* cta_event) const;

  config_type config_;
};

class NectarCamZFITSDataSource:
  public calin::iact_data::zfits_data_source::ZFITSDataSource
{
public:
  CALIN_TYPEALIAS(decoder_config_type,
    NectarCamCameraEventDecoder::config_type);

  void set_decoder_config(const decoder_config_type& config) {
    decoder_->set_config(config); }
  const decoder_config_type& decoder_config() const {
    return decoder_->config(); }
  decoder_config_type* mutable_decoder_config() {
    return decoder_->mutable_config(); }
  static decoder_config_type default_decoder_config() {
    return decoder_config_type::default_instance(); }

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
