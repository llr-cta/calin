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

class CTACameraEventDecoder
{
public:
  virtual ~CTACameraEventDecoder();
  virtual calin::ix::iact_data::telescope_event::TelescopeEvent*
    decode(const DataModel::CameraEvent* cta_event) = 0;
  virtual calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* decode_run_config(
      const DataModel::CameraRunHeader* cta_run_header,
      const DataModel::CameraEvent* cta_event) = 0;
};

class ZFITSSingleFileDataSource:
  public telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig
{
public:
  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig);

  ZFITSSingleFileDataSource(const std::string& filename,
    CTACameraEventDecoder* decoder, bool adopt_decoder = false,
    const config_type& config = default_config());
  virtual ~ZFITSSingleFileDataSource();

  calin::ix::iact_data::telescope_event::TelescopeEvent* get_next() override;
  uint64_t size() override;
  uint64_t next_index() override;
  void set_next_index(uint64_t next_index) override;

  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* get_run_configuration() override;

  static config_type default_config();

private:
  std::string filename_;
  ACTL::IO::ProtobufIFits* zfits_ = nullptr;
  CTACameraEventDecoder* decoder_ = nullptr;
  bool adopt_decoder_ = false;
  uint64_t next_event_index_ = 0;
  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* run_config_ = nullptr;
};

class ZFITSDataSource:
  public calin::io::data_source::BasicChainedRandomAccessDataSource<
    calin::iact_data::telescope_data_source::
      TelescopeRandomAccessDataSourceWithRunConfig>
{
public:
  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig);

  ZFITSDataSource(const std::string& filename,
    CTACameraEventDecoder* decoder, bool adopt_decoder = false,
    const config_type& config = default_config());
  virtual ~ZFITSDataSource();

  calin::ix::iact_data::telescope_event::TelescopeEvent* get_next() override;

  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* get_run_configuration() override;

  static config_type default_config() {
    config_type config = config_type::default_instance();
    config.set_extension(".fits.fz");
    return config;
  }
  const config_type& config() const { return config_; }

private:
  config_type config_;
  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* run_config_ = nullptr;
};

class ZFitsDataSourceOpener:
  public calin::io::data_source::DataSourceOpener<
    calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig>
{
public:
  CALIN_TYPEALIAS(data_source_type,
    calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSource);
  ZFitsDataSourceOpener(std::string filename,
    CTACameraEventDecoder* decoder, bool adopt_decoder = false,
    const ZFITSDataSource::config_type& config =
      ZFITSDataSource::default_config());
  virtual ~ZFitsDataSourceOpener();
  unsigned num_sources() override;
  calin::iact_data::telescope_data_source::
    TelescopeRandomAccessDataSourceWithRunConfig* open(unsigned isource) override;
private:
  std::vector<std::string> filenames_;
  CTACameraEventDecoder* decoder_ = nullptr;
  bool adopt_decoder_ = false;
  ZFITSDataSource::config_type config_;
};


} } } // namespace calin::iact_data::nectarcam_data_source
