/*

   calin/iact_data/nectarcam_data_source.cpp -- Stephen Fegan -- 2016-01-11

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

#include <stdexcept>
#include <memory>
#include <cctype>

#include <io/log.hpp>
#include <util/file.hpp>
#include <iact_data/zfits_data_source.hpp>

using namespace calin::iact_data::zfits_data_source;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::io::log;
using calin::util::file::is_file;
using calin::util::file::is_readable;
using calin::util::file::expand_filename;

#include <ProtobufIFits.h>
#include <L0.pb.h>

CTACameraEventDecoder::~CTACameraEventDecoder()
{
  // nothing to see here
}

// =============================================================================
// ZFITSSingleFileDataSource - single ZFits file with decoder
// Uses ZFITSSingleFileACTLDataSource to read events and decoder to translate
// them.
// =============================================================================

ZFITSSingleFileDataSource::
ZFITSSingleFileDataSource(calin::iact_data::zfits_actl_data_source::
      ZFITSSingleFileACTLDataSource* actl_zfits,
    CTACameraEventDecoder* decoder, bool adopt_decoder,
    const config_type& config, bool adopt_actl_zfits):
  TelescopeRandomAccessDataSourceWithRunConfig(),
  decoder_(decoder), adopt_decoder_(adopt_decoder),
  actl_zfits_(actl_zfits), adopt_actl_zfits_(adopt_actl_zfits)
{
  if(!config.dont_read_run_header())
  {
    DataModel::CameraEvent* actl_sample_event = nullptr;
    DataModel::CameraRunHeader* actl_run_header = nullptr;
    try {
      actl_zfits_->set_next_index(0);
      uint64_t unused_seq_index = 0;
      actl_sample_event = actl_zfits_->borrow_next_event(unused_seq_index);
    } catch(...) {
      // ignore errors that occur reading sample event;
    }
    try {
      actl_run_header = actl_zfits_->get_run_header();
    } catch(...) {
      // ignore errors that occur reading run header
    }
    run_config_ = new TelescopeRunConfiguration;
    decoder_->decode_run_config(run_config_, actl_run_header, actl_sample_event);
    delete actl_run_header;
    if(actl_sample_event)actl_zfits_->release_borrowed_event(actl_sample_event);
    actl_zfits_->set_next_index(0);
  }
}

ZFITSSingleFileDataSource::
ZFITSSingleFileDataSource(const std::string& filename,
    CTACameraEventDecoder* decoder, bool adopt_decoder,
    const config_type& config):
  ZFITSSingleFileDataSource(new calin::iact_data::zfits_actl_data_source::
    ZFITSSingleFileACTLDataSource(filename, config), decoder, adopt_decoder,
    config, true)
{
  // nothing to see here
}

ZFITSSingleFileDataSource::~ZFITSSingleFileDataSource()
{
  delete run_config_;
  if(adopt_actl_zfits_)delete actl_zfits_;
  if(adopt_decoder_)delete decoder_;
}

calin::ix::iact_data::telescope_event::TelescopeEvent*
ZFITSSingleFileDataSource::get_next(
  uint64_t& seq_index_out, google::protobuf::Arena** arena)
{
  DataModel::CameraEvent* cta_event =
    actl_zfits_->borrow_next_event(seq_index_out);
  if(!cta_event)return nullptr;
  TelescopeEvent* event = nullptr;
  TelescopeEvent* delete_event = nullptr;
  google::protobuf::Arena* delete_arena = nullptr;
  if(arena) {
    if(!*arena)*arena = delete_arena = new google::protobuf::Arena;
    event = google::protobuf::Arena::CreateMessage<TelescopeEvent>(*arena);
  }
  else event = delete_event = new TelescopeEvent;
  if(!event)
  {
    delete delete_arena;
    actl_zfits_->release_borrowed_event(cta_event);
    throw std::runtime_error("Could not allocate telescpe event");
  }
  if(!decoder_->decode(event, cta_event))
  {
    delete delete_arena;
    delete delete_event;
    actl_zfits_->release_borrowed_event(cta_event);
    throw std::runtime_error("Could not decode ACTL event");
  }
  actl_zfits_->release_borrowed_event(cta_event);
  event->set_source_event_index(seq_index_out);
  return event;
}

uint64_t ZFITSSingleFileDataSource::size()
{
  return actl_zfits_->size();
}

void ZFITSSingleFileDataSource::set_next_index(uint64_t next_index)
{
  actl_zfits_->set_next_index(next_index);
}

calin::ix::iact_data::telescope_run_configuration::
  TelescopeRunConfiguration* ZFITSSingleFileDataSource::get_run_configuration()
{
  if(!run_config_)return nullptr;
  auto* run_config = new TelescopeRunConfiguration();
  run_config->CopyFrom(*run_config_);
  return run_config;
}

// =============================================================================
// ZFITSDataSource - chained ZFits files with decoder
// =============================================================================

ZFITSDataSource::ZFITSDataSource(const std::string& filename,
    CTACameraEventDecoder* decoder, bool adopt_decoder,
    const config_type& config):
  calin::io::data_source::BasicChainedRandomAccessDataSource<
    calin::iact_data::telescope_data_source::
      TelescopeRandomAccessDataSourceWithRunConfig>(
    new ZFITSDataSourceOpener(filename, decoder, config), true),
  decoder_(decoder), adopt_decoder_(adopt_decoder),
  run_config_(source_->get_run_configuration())
{
  // nothing to see here
}

ZFITSDataSource::~ZFITSDataSource()
{
  delete run_config_;
  if(adopt_decoder_)delete decoder_;
}

calin::ix::iact_data::telescope_run_configuration::
TelescopeRunConfiguration* ZFITSDataSource::get_run_configuration()
{
  if(!run_config_)return nullptr;
  auto* run_config = new TelescopeRunConfiguration();
  run_config->CopyFrom(*run_config_);
  return run_config;
}

// =============================================================================
// ZFITSDataSourceOpener - opener for ZFITSSingleFileDataSource
// Uses ZFITSACTLDataSourceOpener to open underlying ZFITSSingleFileACTLDataSource
// objects
// =============================================================================

ZFITSDataSourceOpener::
ZFITSDataSourceOpener(std::string filename, CTACameraEventDecoder* decoder,
    const ZFITSDataSource::config_type& config):
  calin::io::data_source::DataSourceOpener<
    calin::iact_data::telescope_data_source::
      TelescopeRandomAccessDataSourceWithRunConfig>(),
  zfits_actl_opener_(new calin::iact_data::zfits_actl_data_source::
    ZFITSACTLDataSourceOpener(filename, config)),
  decoder_(decoder), config_(config)
{
  // nothing to see here
}

ZFITSDataSourceOpener::~ZFITSDataSourceOpener()
{
  delete zfits_actl_opener_;
}

unsigned ZFITSDataSourceOpener::num_sources()
{
  return zfits_actl_opener_->num_sources();
}

ZFITSSingleFileDataSource* ZFITSDataSourceOpener::open(unsigned isource)
{
  auto* zfits_actl = zfits_actl_opener_->open(isource);
  if(zfits_actl == nullptr)return nullptr;
  auto config = config_;
  if(zfits_actl_opener_->has_opened_file())config.set_dont_read_run_header(true);
  return new ZFITSSingleFileDataSource(zfits_actl, decoder_, false, config, true);
}
