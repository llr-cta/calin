/*

   calin/iact_data/cta_data_source.hpp -- Stephen Fegan -- 2018-11-06

   A supplier of single telescope data from CTA DAQ data files

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
#include <iact_data/cta_data_source.pb.h>
#include <iact_data/zfits_data_source.hpp>
#include <iact_data/zfits_actl_data_source.hpp>
#include <pattern/delegation.hpp>

namespace calin { namespace iact_data { namespace cta_data_source {

#ifdef CALIN_HAVE_CTA_CAMERASTOACTL

class CTAZFITSDataSource:
  public calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig,
  public calin::pattern::delegation::Delegator<
    calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig>
{
public:
  CALIN_TYPEALIAS(decoder_config_type,
    calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig);
  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig);

  static calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig default_config() {
    config_type config = calin::iact_data::zfits_data_source::ZFITSDataSource_R1::default_config();
    config.set_data_model(calin::ix::iact_data::zfits_data_source::ACTL_DATA_MODEL_AUTO_DETECT);
    config.set_run_header_table_name(""); // Differs between L0 and R1 so let downstream decode
    return config;
  }

  static calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig default_decoder_config();

  CTAZFITSDataSource(const std::string& filename,
    const calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig& config,
    const calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig& decoder_config = default_decoder_config());
  CTAZFITSDataSource(const std::string& filename,
    const calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig& decoder_config = default_decoder_config(),
    const calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig& config = default_config());
  virtual ~CTAZFITSDataSource();

  calin::ix::iact_data::telescope_event::TelescopeEvent* get_next(
    uint64_t& seq_index_out, google::protobuf::Arena** arena = nullptr) override;
  uint64_t size() override;
  void set_next_index(uint64_t next_index) override;

  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* get_run_configuration() override;

private:
  static TelescopeRandomAccessDataSourceWithRunConfig* construct_delegate(
    const std::string& filename, config_type config, decoder_config_type decoder_config);
};


#endif

} } } // namespace calin::iact_data::cta_data_source
