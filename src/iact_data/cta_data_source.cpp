/*

   calin/iact_data/cta_data_source.cpp -- Stephen Fegan -- 2018-11-06

   A supplier of single telescope data from CTA ZFITS data files

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

CTAZFITSDataSource::
CTAZFITSDataSource(const std::string& filename,
    CTAZFITSDataSource* base_data_source, const config_type& config):
  TelescopeRandomAccessDataSourceWithRunConfig(),
  Delegator<calin::iact_data::telescope_data_source::
    TelescopeRandomAccessDataSourceWithRunConfig>(
      copy_base_data_source(filename, config, base_data_source), /* adopt_delegate=*/ true)
{
  // nothing to see here
}

CTAZFITSDataSource::
CTAZFITSDataSource(const std::string& filename,
    const config_type& config, CTAZFITSDataSource* base_data_source):
  CTAZFITSDataSource(filename, base_data_source, config)
{
  // nothing to see here
}

CTAZFITSDataSource::~CTAZFITSDataSource()
{
  // nothing to see here
}

calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig
CTAZFITSDataSource::default_config()
{
  config_type config = calin::iact_data::zfits_data_source::ZFITSDataSource_R1::default_config();
  config.set_data_model(calin::ix::iact_data::zfits_data_source::ACTL_DATA_MODEL_AUTO_DETECT);
  config.set_run_header_table_name(""); // Differs between L0 and R1 so let downstream decode
  return config;
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

unsigned CTAZFITSDataSource::current_fragment_index() const
{
  return dynamic_cast<const calin::io::data_source::FragmentList&>(
    *this->delegate()).current_fragment_index();
}

unsigned CTAZFITSDataSource::num_fragments() const
{
  return dynamic_cast<const calin::io::data_source::FragmentList&>(
    *this->delegate()).num_fragments();
}

std::string CTAZFITSDataSource::fragment_name(unsigned index) const
{
  return dynamic_cast<const calin::io::data_source::FragmentList&>(
    *this->delegate()).fragment_name(index);
}

calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig*
CTAZFITSDataSource::construct_delegate(std::string filename,
  config_type config, decoder_config_type decoder_config)
{
  expand_filename_in_place(filename);

  if(!is_file(filename))
    throw std::runtime_error(
      "CTAZFITSDataSource::construct_delegate: File not found: " + filename);

  if(config.data_model() ==
      calin::ix::iact_data::zfits_data_source::ACTL_DATA_MODEL_AUTO_DETECT) {
    if(zfits_actl_data_source::is_zfits_r1(filename, config.events_table_name())) {
      config.set_data_model(calin::ix::iact_data::zfits_data_source::ACTL_DATA_MODEL_R1);
    } else if (zfits_actl_data_source::is_zfits_l0(filename, config.events_table_name())) {
      config.set_data_model(calin::ix::iact_data::zfits_data_source::ACTL_DATA_MODEL_L0);
    }
  }

  if(config.data_model() == calin::ix::iact_data::zfits_data_source::ACTL_DATA_MODEL_L0) {
    if(decoder_config.camera_type() ==
          calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig::AUTO_DETECT
        or decoder_config.camera_type() ==
          calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig::NECTARCAM) {
      return new NectarCamZFITSDataSource_L0(filename, config, decoder_config.nectarcam());
    }
    throw std::runtime_error(
      "CTAZFITSDataSource::construct_delegate: L0 data format not supported with camera type: "
      + calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig::CameraType_Name(decoder_config.camera_type()));
  } else if(config.data_model() == calin::ix::iact_data::zfits_data_source::ACTL_DATA_MODEL_R1) {
    if(decoder_config.camera_type()
        == calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig::AUTO_DETECT)
    {
      ZFITSSingleFileACTL_R1_CameraEventDataSource src(filename, config);
      auto* header = src.get_run_header();
      if(header == nullptr) {
        throw std::runtime_error(
          "CTAZFITSDataSource::construct_delegate: no run header found: " + filename);
      }
      if(header->has_nectarcam()) {
        decoder_config.set_camera_type(
          calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig::NECTARCAM);
      } else if(header->has_lstcam()) {
        decoder_config.set_camera_type(
          calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig::LSTCAM);
      } else {
        delete header;
        throw std::runtime_error(
          "CTAZFITSDataSource::construct_delegate: no known camera extension found");
      }
      delete header;
    }

    if(decoder_config.camera_type() == calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig::NECTARCAM) {
      return new NectarCamZFITSDataSource_R1(filename, config, decoder_config.nectarcam());
    } else {
      return new LSTCamZFITSDataSource_R1(filename, config, decoder_config.lstcam());
    }
  } else {
    throw std::runtime_error(
      "CTAZFITSDataSource::construct_delegate: unsupported data format: "
      + calin::ix::iact_data::zfits_data_source::ACTLDataModel_Name(config.data_model()));
  }
  return nullptr;
}

calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig*
CTAZFITSDataSource::copy_base_data_source(
  std::string filename, config_type config, CTAZFITSDataSource* base_data_source)
{
  expand_filename_in_place(filename);

  if(!is_file(filename))
    throw std::runtime_error(
      "CTAZFITSDataSource::copy_base_data_source: File not found: " + filename);

  if(auto* zfits = dynamic_cast<zfits_data_source::ZFITSDataSource_R1*>(base_data_source->delegate())) {
    return new zfits_data_source::ZFITSDataSource_R1(filename, zfits, config);
  }else if(auto* zfits = dynamic_cast<zfits_data_source::ZFITSDataSource_L0*>(base_data_source->delegate())) {
    return new zfits_data_source::ZFITSDataSource_L0(filename, zfits, config);
  } else {
    throw std::runtime_error(
      "CTAZFITSDataSource::copy_base_data_source: unsupported data source");
  }
}
