/*

   calin/iact_data/zfits_data_source.hpp -- Stephen Fegan -- 2016-01-21

   A supplier of single telescope data from CTA ACTL ZFits data files

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
#include <io/chained_data_source.hpp>
#include <iact_data/zfits_actl_data_source.hpp>
#include <iact_data/zfits_data_source.pb.h>
#include <iact_data/telescope_data_source.hpp>

namespace calin { namespace iact_data { namespace zfits_data_source {

#ifdef CALIN_HAVE_CTA_CAMERASTOACTL

class CTACameraEventDecoder
{
public:
  virtual ~CTACameraEventDecoder();
  virtual bool decode(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    const DataModel::CameraEvent* cta_event) = 0;
  virtual bool decode_run_config(
    calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config,
    const DataModel::CameraRunHeader* cta_run_header,
    const DataModel::CameraEvent* cta_event) = 0;
};

class ZFITSDataSource:
  public calin::iact_data::telescope_data_source::
    TelescopeRandomAccessDataSourceWithRunConfig
{
public:
  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig);

  ZFITSDataSource(const std::string& filename,
    CTACameraEventDecoder* decoder, bool adopt_decoder = false,
    const config_type& config = default_config());
  virtual ~ZFITSDataSource();

  calin::ix::iact_data::telescope_event::
  TelescopeEvent* get_next(uint64_t& seq_index_out,
      google::protobuf::Arena** arena = nullptr) override;

  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* get_run_configuration() override;

  uint64_t size() override;
  void set_next_index(uint64_t next_index) override;

  static config_type default_config() {
    return zfits_actl_data_source::ZFITSACTLDataSource::default_config(); }

private:
  CTACameraEventDecoder* decoder_;
  bool adopt_decoder_ = false;
  calin::iact_data::zfits_actl_data_source::
    ZFITSACTLDataSource* actl_zfits_ = nullptr;
  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* run_config_ = nullptr;
};

#endif
} } } // namespace calin::iact_data::nectarcam_data_source
