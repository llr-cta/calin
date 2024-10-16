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
#include <iact_data/nectarcam_acada_event_decoder.hpp>
#include <iact_data/nectarcam_data_source.hpp>
#include <iact_data/nectarcam_data_source.pb.h>
#include <iact_data/lstcam_acada_event_decoder.hpp>
#include <iact_data/lstcam_data_source.hpp>
#include <iact_data/lstcam_data_source.pb.h>

using namespace calin::ix::iact_data::cta_data_source;
using namespace calin::iact_data::cta_data_source;
using namespace calin::iact_data::nectarcam_data_source;
using namespace calin::ix::iact_data::nectarcam_data_source;
using namespace calin::iact_data::nectarcam_acada_event_decoder;
using namespace calin::iact_data::lstcam_data_source;
using namespace calin::ix::iact_data::lstcam_data_source;
using namespace calin::iact_data::lstcam_acada_event_decoder;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::iact_data::zfits_acada_data_source;
using namespace calin::util::log;
using namespace calin::util::file;

CTAZFITSDataSource::
CTAZFITSDataSource(const std::string& filename,
    const config_type& config, const decoder_config_type& decoder_config):
  BasicZFITSDataSource(),
  Delegator<calin::iact_data::zfits_data_source::BasicZFITSDataSource>(
      construct_delegate(filename, config, decoder_config), /* adopt_delegate=*/ true),
  config_(config), decoder_config_(decoder_config)
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

calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig
CTAZFITSDataSource::default_config()
{
  config_type config = calin::iact_data::zfits_data_source::ZFITSDataSource_R1v1::default_config();
  config.set_data_model(calin::ix::iact_data::zfits_data_source::ACADA_DATA_MODEL_AUTO_DETECT);
  config.set_events_table_name(""); // Differs between R1v0 and R1v1 so let downstream decide
  config.set_run_header_table_name(""); // Differs between L0, R1v0 and R1v1 so let downstream decide
  config.set_data_stream_table_name(""); // Differs between L0, R1v0 and R1v1 so let downstream decide
  return config;
}

CTACameraEventDecoderConfig
CTAZFITSDataSource::default_decoder_config()
{
  CTACameraEventDecoderConfig config;
  config.mutable_nectarcam()->CopyFrom(
    NectarCamZFITSDataSource_R1v0::default_decoder_config());
  config.mutable_lstcam()->CopyFrom(
    LSTCamZFITSDataSource_R1v0::default_decoder_config());
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
  return delegate_->current_fragment_index();
}

unsigned CTAZFITSDataSource::num_fragments() const
{
  return delegate_->num_fragments();
}

unsigned CTAZFITSDataSource::num_missing_fragments() const
{
  return delegate_->num_missing_fragments();
}

std::string CTAZFITSDataSource::fragment_name(unsigned index) const
{
  return delegate_->fragment_name(index);
}

calin::iact_data::zfits_data_source::BasicZFITSDataSource* 
CTAZFITSDataSource::new_of_type(const std::string& filename, const config_type& config,
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* force_run_configuration)
{
  return delegate_->new_of_type(filename, config, force_run_configuration);
}

calin::iact_data::zfits_data_source::BasicZFITSDataSource*
CTAZFITSDataSource::construct_delegate(std::string filename,
  config_type config, decoder_config_type decoder_config)
{
  filename = expand_filename(filename);

  std::vector<std::string> fragment_filenames;
  unsigned unused_num_missing_fragments = 0;
  config_type fragment_find_config = config;
  fragment_find_config.set_max_file_fragments(1);
  calin::iact_data::zfits_acada_data_source::generate_fragment_list(filename, fragment_find_config,
    fragment_filenames, unused_num_missing_fragments, /* log_missing_fragments = */ false);

  if(fragment_filenames.empty()) {
    throw std::runtime_error(
      "CTAZFITSDataSource::construct_delegate: File not found: " + expand_filename(filename));
  }

  if(config.data_model() ==
      calin::ix::iact_data::zfits_data_source::ACADA_DATA_MODEL_AUTO_DETECT) {
    config.set_data_model(calin::iact_data::zfits_acada_data_source::get_zfits_data_model(fragment_filenames[0]));

    if(config.data_model() ==
        calin::ix::iact_data::zfits_data_source::ACADA_DATA_MODEL_AUTO_DETECT) {
      throw std::runtime_error(
        "CTAZFITSDataSource::construct_delegate: could not auto-detect data model: " + filename);
    }
  }

  if(config.data_model() == calin::ix::iact_data::zfits_data_source::ACADA_DATA_MODEL_L0) {
    if(decoder_config.camera_type() == AUTO_DETECT
        or decoder_config.camera_type() == NECTARCAM) {
      return new NectarCamZFITSDataSource_L0(filename, config, decoder_config.nectarcam());
    }
    throw std::runtime_error(
      "CTAZFITSDataSource::construct_delegate: L0 data format not supported with camera type: "
      + CameraType_Name(decoder_config.camera_type()));
  } else if(config.data_model() == calin::ix::iact_data::zfits_data_source::ACADA_DATA_MODEL_R1V0) {
    if(decoder_config.camera_type()
        == AUTO_DETECT)
    {
      auto* src = new ZFITSACADACameraEventDataSource_R1v0(filename, config);
      auto* header = src->get_run_header();
      if(header == nullptr) {
        delete src;
        throw std::runtime_error(
          "CTAZFITSDataSource::construct_delegate: no run header found: " + filename);
      }
      if(header->has_nectarcam()) {
        decoder_config.set_camera_type(
          NECTARCAM);
        delete header;
        return new NectarCamZFITSDataSource_R1v0(filename, src, decoder_config.nectarcam(), /* adopt_zfits_src= */ true);
      } else if(header->has_lstcam()) {
        decoder_config.set_camera_type(
          LSTCAM);
        delete header;
        return new LSTCamZFITSDataSource_R1v0(filename, src, decoder_config.lstcam(), /* adopt_zfits_src= */ true);
      } else {
        delete header;
        delete src;
        throw std::runtime_error(
          "CTAZFITSDataSource::construct_delegate: no known camera extension found");
      }
      delete header;
    } else if(decoder_config.camera_type() == NECTARCAM) {
      return new NectarCamZFITSDataSource_R1v0(filename, config, decoder_config.nectarcam());
    } else if(decoder_config.camera_type() == LSTCAM) {
      return new LSTCamZFITSDataSource_R1v0(filename, config, decoder_config.lstcam());
    } else {
      throw std::runtime_error(
        "CTAZFITSDataSource::construct_delegate: R1v0 data format not supported with camera type: "
        + CameraType_Name(decoder_config.camera_type()));
    }
  } else if(config.data_model() == calin::ix::iact_data::zfits_data_source::ACADA_DATA_MODEL_R1V1) {
      if(decoder_config.camera_type() == AUTO_DETECT) {
        return new NectarCamZFITSDataSource_R1v1(filename, config, decoder_config.nectarcam());
      } else if(decoder_config.camera_type() == NECTARCAM) {
        return new NectarCamZFITSDataSource_R1v1(filename, config, decoder_config.nectarcam());
      } else {
        throw std::runtime_error(
          "CTAZFITSDataSource::construct_delegate: R1v1 data format not supported with camera type: "
          + CameraType_Name(decoder_config.camera_type()));
      }
  } else {
    throw std::runtime_error(
      "CTAZFITSDataSource::construct_delegate: unsupported data format: "
      + calin::ix::iact_data::zfits_data_source::ACADADataModel_Name(config.data_model()));
  }
  return nullptr;
}

CTAZFITSDataSourceFactory::~CTAZFITSDataSourceFactory()
{
  if(adopt_base_data_source_)delete base_data_source_;
  if(adopt_run_config_)delete run_config_;
}

calin::io::data_source::DataSource< calin::ix::iact_data::telescope_event::TelescopeEvent>* 
CTAZFITSDataSourceFactory::new_data_source()
{
  unsigned isource = isource_.fetch_add(1,std::memory_order_relaxed);
  if(isource<base_data_source_->num_fragments()) {
    auto config = base_data_source_->config();
    config.clear_forced_file_fragments_list();
    config.add_forced_file_fragments_list(base_data_source_->fragment_name(isource));
    return base_data_source_->new_of_type(config.forced_file_fragments_list(0), config, run_config_);
  } else {
    return nullptr;
  }
}
