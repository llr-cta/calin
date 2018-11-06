/*

   calin/iact_data/cta_data_source.cpp -- Stephen Fegan -- 2018-11-06

   A supplier of single telescope data from CTA ZFITS data files

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

#include <stdexcept>
#include <string>
#include <memory>
#include <numeric>

#include <util/log.hpp>
#include <util/file.hpp>
#include <iact_data/cta_data_source.hpp>
#include <iact_data/nectarcam_actl_event_decoder.hpp>
#include <iact_data/nectarcam_data_source.hpp>
#include <iact_data/nectarcam_data_source.pb.h>
#include <iact_data/lstcam_actl_event_decoder.hpp>
#include <iact_data/lstcam_data_source.hpp>
#include <iact_data/lstcam_data_source.pb.h>

using namespace calin::iact_data::cta_data_source;
using namespace calin::iact_data::nectarcam_data_source;
using namespace calin::ix::iact_data::nectarcam_data_source;
using namespace calin::iact_data::nectarcam_actl_event_decoder;
using namespace calin::iact_data::lstcam_data_source;
using namespace calin::ix::iact_data::lstcam_data_source;
using namespace calin::iact_data::lstcam_actl_event_decoder;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::iact_data::zfits_actl_data_source;
using namespace calin::util::log;
using namespace calin::util::file;

#include <ProtobufIFits.h>
#include <IFits.h>
#include <L0.pb.h>
#include <R1.pb.h>

CTAZFITSDataSource::
CTAZFITSDataSource(const std::string& filename,
    const config_type& config, const decoder_config_type& decoder_config):
  TelescopeRandomAccessDataSourceWithRunConfig(),
  Delegator<calin::iact_data::telescope_data_source::
    TelescopeRandomAccessDataSourceWithRunConfig>(
      construct_delegate(filename, config, decoder_config), /* adopt_delegate=*/ true)
{
  // nothing to see here
}

CTAZFITSDataSource::
CTAZFITSDataSource(const std::string& filename,
    const decoder_config_type& decoder_config, const config_type& config):
  CTAZFITSDataSource(filename, config, decoder_config)
{
  // nothing to see here
}

CTAZFITSDataSource::~CTAZFITSDataSource()
{
  // nothing to see here
}

calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig
CTAZFITSDataSource::default_decoder_config()
{
  calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig config;
  config.mutable_nectarcam()->CopyFrom(
    NectarCamZFITSDataSource_R1::default_decoder_config());
  config.mutable_lstcam()->CopyFrom(
    LSTCamZFITSDataSource_R1::default_decoder_config());
  return config;
}

calin::ix::iact_data::telescope_event::TelescopeEvent*
CTAZFITSDataSource::get_next(uint64_t& seq_index_out, google::protobuf::Arena** arena)
{
  return delegate_->get_next(seq_index_out, arena);
}

uint64_t CTAZFITSDataSource::size()
{
  return delegate_->size();
}

void CTAZFITSDataSource::set_next_index(uint64_t next_index)
{
  return delegate_->set_next_index(next_index);
}

calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration*
CTAZFITSDataSource::get_run_configuration()
{
  return delegate_->get_run_configuration();
}

calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig*
CTAZFITSDataSource::construct_delegate(const std::string& filename,
  const config_type& config, decoder_config_type decoder_config)
{
  if(config.data_model() ==
      calin::ix::iact_data::zfits_data_source::ACTL_DATA_MODEL_L0) {
    throw std::runtime_error(
      "CTAZFITSDataSource::construct_delegate: L0 data format not supported");
  }

  if(decoder_config.camera_type()
      == calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig::AUTOMATIC)
  {
    if(!is_file(filename))
      throw std::runtime_error(
        "CTAZFITSDataSource::construct_delegate: File not found: " + filename);
    ZFITSSingleFileACTL_R1_CameraEventDataSource src(filename, config);
    auto* header = src.get_run_header();
    if(header->has_nectarcam())
      decoder_config.set_camera_type(
        calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig::NECTARCAM);
    else if(header->has_lstcam())
      decoder_config.set_camera_type(
        calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig::LSTCAM);
    else
      throw std::runtime_error(
        "CTAZFITSDataSource::construct_delegate: no known camera extension found");
    delete header;
  }

  if(decoder_config.camera_type() == calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig::NECTARCAM) {
    return new NectarCamZFITSDataSource_R1(filename, config, decoder_config.nectarcam());
  } else {
    return new LSTCamZFITSDataSource_R1(filename, config, decoder_config.lstcam());
  }
}
