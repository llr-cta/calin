/*

   calin/iact_data/cta_data_source.hpp -- Stephen Fegan -- 2018-11-06

   A supplier of single telescope data from CTA DAQ data files

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
#include <iact_data/cta_data_source.pb.h>
#include <iact_data/zfits_data_source.hpp>
#include <iact_data/zfits_actl_data_source.hpp>
#include <pattern/delegation.hpp>

namespace calin { namespace iact_data { namespace cta_data_source {

#ifdef CALIN_HAVE_CTA_CAMERASTOACTL

class CTAZFITSDataSource:
  public calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig,
  public calin::io::data_source::FragmentList,
  private calin::pattern::delegation::Delegator<
    calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig>
{
public:
  CALIN_TYPEALIAS(decoder_config_type,
    calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig);
  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig);

  static calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig default_config();
  static calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig default_decoder_config();

  CTAZFITSDataSource(const std::string& filename,
    const calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig& config,
    const calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig& decoder_config = default_decoder_config());
  CTAZFITSDataSource(const std::string& filename,
    const calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig& decoder_config = default_decoder_config(),
    const calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig& config = default_config());

  CTAZFITSDataSource(const std::string& filename, CTAZFITSDataSource* base_data_source,
    const calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig& config = default_config());
  CTAZFITSDataSource(const std::string& filename,
    const calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig& config,
    CTAZFITSDataSource* base_data_source);

  virtual ~CTAZFITSDataSource();

  calin::ix::iact_data::telescope_event::TelescopeEvent* get_next(
    uint64_t& seq_index_out, google::protobuf::Arena** arena = nullptr) override;
  uint64_t size() override;
  void set_next_index(uint64_t next_index) override;

  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* get_run_configuration() override;

  unsigned current_fragment_index() const override;
  unsigned num_fragments() const override;
  std::string fragment_name(unsigned index) const override;

private:
  static TelescopeRandomAccessDataSourceWithRunConfig* construct_delegate(
    std::string filename, config_type config, decoder_config_type decoder_config);
  static TelescopeRandomAccessDataSourceWithRunConfig* copy_base_data_source(
    std::string filename, config_type config, CTAZFITSDataSource* base_data_source);
};


#endif

} } } // namespace calin::iact_data::cta_data_source
