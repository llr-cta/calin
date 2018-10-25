/*

   calin/iact_data/nectarcam_data_source.cpp -- Stephen Fegan -- 2016-01-11

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

#include <stdexcept>
#include <memory>
#include <cctype>

#include <util/log.hpp>
#include <util/file.hpp>
#include <iact_data/actl_event_decoder.hpp>
#include <iact_data/zfits_data_source.hpp>

using namespace calin::iact_data::zfits_data_source;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::util::log;
using calin::util::file::is_file;
using calin::util::file::is_readable;
using calin::util::file::expand_filename;

#include <ProtobufIFits.h>
#include <R1.pb.h>

using namespace calin::iact_data::actl_event_decoder;

/*

                      RRRRRRRRRRRRRRRRR          1111111
                      R::::::::::::::::R        1::::::1
                      R::::::RRRRRR:::::R      1:::::::1
                      RR:::::R     R:::::R     111:::::1
                        R::::R     R:::::R        1::::1
                        R::::R     R:::::R        1::::1
                        R::::RRRRRR:::::R         1::::1
                        R:::::::::::::RR          1::::l
                        R::::RRRRRR:::::R         1::::l
                        R::::R     R:::::R        1::::l
                        R::::R     R:::::R        1::::l
                        R::::R     R:::::R        1::::l
                      RR:::::R     R:::::R     111::::::111
                      R::::::R     R:::::R     1::::::::::1
                      R::::::R     R:::::R     1::::::::::1
                      RRRRRRRR     RRRRRRR     111111111111

*/

// =============================================================================
// ZFITSSingleFileDataSource_R1 - single ZFits file with decoder
// Uses ZFITSSingleFileACTL_R1_CameraEventDataSource to read events and
// decoder to translate them.
// =============================================================================

ZFITSSingleFileDataSource_R1::
ZFITSSingleFileDataSource_R1(calin::iact_data::zfits_actl_data_source::
      ZFITSSingleFileACTL_R1_CameraEventDataSource* actl_zfits,
    bool dont_decode_run_configuration,
    ACTL_R1_CameraEventDecoder* decoder, bool adopt_decoder, bool adopt_actl_zfits):
  TelescopeRandomAccessDataSourceWithRunConfig(),
  decoder_(decoder), adopt_decoder_(adopt_decoder),
  actl_zfits_(actl_zfits), adopt_actl_zfits_(adopt_actl_zfits)
{
  if(not dont_decode_run_configuration)
  {
    const R1::CameraEvent* actl_sample_event = nullptr;
    const R1::CameraConfiguration* actl_run_header = nullptr;
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

ZFITSSingleFileDataSource_R1::
ZFITSSingleFileDataSource_R1(const std::string& filename,
    ACTL_R1_CameraEventDecoder* decoder, bool adopt_decoder,
    const config_type& config):
  ZFITSSingleFileDataSource_R1(new calin::iact_data::zfits_actl_data_source::
    ZFITSSingleFileACTL_R1_CameraEventDataSource(filename, config), false,
    decoder, adopt_decoder, true)
{
  // nothing to see here
}

ZFITSSingleFileDataSource_R1::~ZFITSSingleFileDataSource_R1()
{
  delete run_config_;
  if(adopt_actl_zfits_)delete actl_zfits_;
  if(adopt_decoder_)delete decoder_;
}

calin::ix::iact_data::telescope_event::TelescopeEvent*
ZFITSSingleFileDataSource_R1::get_next(
  uint64_t& seq_index_out, google::protobuf::Arena** arena)
{
  const R1::CameraEvent* cta_event =
    actl_zfits_->borrow_next_event(seq_index_out);
  if(!cta_event){
    if(arena)*arena = nullptr;
    return nullptr;
  }
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

uint64_t ZFITSSingleFileDataSource_R1::size()
{
  return actl_zfits_->size();
}

void ZFITSSingleFileDataSource_R1::set_next_index(uint64_t next_index)
{
  actl_zfits_->set_next_index(next_index);
}

calin::ix::iact_data::telescope_run_configuration::
  TelescopeRunConfiguration* ZFITSSingleFileDataSource_R1::get_run_configuration()
{
  if(!run_config_)return nullptr;
  auto* run_config = new TelescopeRunConfiguration();
  run_config->CopyFrom(*run_config_);
  return run_config;
}

// =============================================================================
// ZFITSDataSource_R1 - chained ZFits files with decoder
// =============================================================================

ZFITSDataSource_R1::ZFITSDataSource_R1(const std::string& filename,
    ACTL_R1_CameraEventDecoder* decoder, bool adopt_decoder,
    const config_type& config):
  calin::io::data_source::BasicChainedRandomAccessDataSource<
    calin::iact_data::telescope_data_source::
      TelescopeRandomAccessDataSourceWithRunConfig>(
    new ZFITSDataSourceOpener_R1(filename, decoder, config), true),
  decoder_(decoder), adopt_decoder_(adopt_decoder),
  run_config_(source_->get_run_configuration())
{
  if(run_config_) {
    run_config_->clear_fragment_filename();
    for(const auto& ifilename : source_names())
      run_config_->add_fragment_filename(ifilename);
  }
}

ZFITSDataSource_R1::~ZFITSDataSource_R1()
{
  delete run_config_;
  if(adopt_decoder_)delete decoder_;
}

calin::ix::iact_data::telescope_run_configuration::
TelescopeRunConfiguration* ZFITSDataSource_R1::get_run_configuration()
{
  if(!run_config_)return nullptr;
  auto* run_config = new TelescopeRunConfiguration();
  run_config->CopyFrom(*run_config_);
  return run_config;
}

// =============================================================================
// ZFITSDataSourceOpener_R1 - opener for ZFITSSingleFileDataSource_R1
// Uses ZFITSACTL_R1_CameraEventDataSourceOpener to open underlying ZFITSSingleFileACTL_R1_CameraEventDataSource
// objects
// =============================================================================

ZFITSDataSourceOpener_R1::
ZFITSDataSourceOpener_R1(std::string filename, ACTL_R1_CameraEventDecoder* decoder,
    const ZFITSDataSource_R1::config_type& config):
  calin::io::data_source::DataSourceOpener<
    calin::iact_data::telescope_data_source::
      TelescopeRandomAccessDataSourceWithRunConfig>(),
  zfits_actl_opener_(new calin::iact_data::zfits_actl_data_source::
    ZFITSACTL_R1_CameraEventDataSourceOpener(filename, config)),
  decoder_(decoder), config_(config)
{
  // nothing to see here
}

ZFITSDataSourceOpener_R1::~ZFITSDataSourceOpener_R1()
{
  delete zfits_actl_opener_;
}

unsigned ZFITSDataSourceOpener_R1::num_sources()
{
  return zfits_actl_opener_->num_sources();
}

std::string ZFITSDataSourceOpener_R1::source_name(unsigned isource)
{
  return zfits_actl_opener_->source_name(isource);
}

ZFITSSingleFileDataSource_R1* ZFITSDataSourceOpener_R1::open(unsigned isource)
{
  bool suppress_run_config = zfits_actl_opener_->has_opened_file();
  auto* zfits_actl = zfits_actl_opener_->open(isource);
  if(zfits_actl == nullptr)return nullptr;
  return new ZFITSSingleFileDataSource_R1(zfits_actl, suppress_run_config,
     decoder_, false, true);
}