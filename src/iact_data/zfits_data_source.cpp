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
#include <iact_data/zfits_data_source.hpp>

using namespace calin::iact_data::zfits_data_source;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::util::log;
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
// DecodedACTLDataSource - extract ACTL events from an ACTLDataSource and
// decode them using a supplied decoder
// =============================================================================

DecodedACTLDataSource::DecodedACTLDataSource(
    calin::iact_data::zfits_actl_data_source::ACTLDataSource* actl_src,
    CTACameraEventDecoder* decoder, bool adopt_actl_src, bool adopt_decoder):
  calin::iact_data::telescope_data_source::TelescopeDataSource(),
  decoder_(decoder), adopt_decoder_(adopt_decoder),
  actl_src_(actl_src), adopt_actl_src_(adopt_actl_src)
{
  // nothing to see here
}

DecodedACTLDataSource::~DecodedACTLDataSource()
{
  if(adopt_decoder_)delete decoder_;
  if(adopt_actl_src_)delete actl_src_;
}

calin::ix::iact_data::telescope_event::TelescopeEvent*
DecodedACTLDataSource::get_next(
  uint64_t& seq_index_out, google::protobuf::Arena** arena)
{
  const DataModel::CameraEvent* cta_event = actl_src_->get_next(seq_index_out);
  if(!cta_event) {
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
    delete cta_event;
    delete delete_arena;
    throw std::runtime_error("Could not allocate telescope event");
  }
  if(!decoder_->decode(event, cta_event))
  {
    delete cta_event;
    delete delete_arena;
    delete delete_event;
    throw std::runtime_error("Could not decode ACTL event");
  }
  event->set_source_event_index(seq_index_out);
  delete cta_event;
  return event;
}

// =============================================================================
// ZFITSSingleFileDataSource - single ZFits file with decoder
// Uses ZFITSSingleFileACTLDataSource to read events and decoder to translate
// them.
// =============================================================================

ZFITSSingleFileDataSource::
ZFITSSingleFileDataSource(calin::iact_data::zfits_actl_data_source::
      ZFITSSingleFileACTLDataSource* actl_zfits,
    bool dont_decode_run_configuration,
    CTACameraEventDecoder* decoder, bool adopt_decoder, bool adopt_actl_zfits):
  TelescopeRandomAccessDataSourceWithRunConfig(),
  decoder_(decoder), adopt_decoder_(adopt_decoder),
  actl_zfits_(actl_zfits), adopt_actl_zfits_(adopt_actl_zfits)
{
  if(not dont_decode_run_configuration)
  {
    const DataModel::CameraEvent* actl_sample_event = nullptr;
    const DataModel::CameraRunHeader* actl_run_header = nullptr;
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
    ZFITSSingleFileACTLDataSource(filename, config), false,
    decoder, adopt_decoder, true)
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
  const DataModel::CameraEvent* cta_event =
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
  bool suppress_run_config = zfits_actl_opener_->has_opened_file();
  auto* zfits_actl = zfits_actl_opener_->open(isource);
  if(zfits_actl == nullptr)return nullptr;
  return new ZFITSSingleFileDataSource(zfits_actl, suppress_run_config,
     decoder_, false, true);
}

// ----------------------------------------------------------------------------
//
// Utility decoding functions
//
// ----------------------------------------------------------------------------

void calin::iact_data::zfits_data_source::
decode_cdts_data(calin::ix::iact_data::telescope_event::CDTSData* calin_cdts_data,
  const DataModel::AnyArray& cta_array)
{
  // Reference : https://forge.in2p3.fr/projects/cta/repository/entry/ACTL/ExternalDevicesCommunication/trunk/TiCkSdecode/ticks_decode.c
  struct CDTSMessageData {
    uint32_t event_counter;
    uint32_t pps_counter;
    uint32_t clock_counter;
    uint64_t ucts_timestamp;
    uint64_t camera_timestamp;
    uint8_t trigger_type;
    uint8_t white_rabbit_status;
    uint8_t arbitrary_information;
  } __attribute__((__packed__));

  const auto& cta_cdts_data = cta_array.data();
#if TEST_ANYARRAY_TYPES
  if(cta_cdts_data.type() != DataModel::AnyArray::U32)
    throw std::runtime_error("CDTS counters type not U32");
#endif
  if(cta_cdts_data.size() != sizeof(CDTSMessageData))
    throw std::runtime_error("CDTS data array not expected size");
  const auto* cdts_data =
    reinterpret_cast<const CDTSMessageData*>(&cta_cdts_data.front());

  calin_cdts_data->set_event_counter(cdts_data->event_counter);
  calin_cdts_data->set_pps_counter(cdts_data->pps_counter);
  calin_cdts_data->set_clock_counter(cdts_data->clock_counter);
  calin_cdts_data->set_ucts_timestamp(cdts_data->ucts_timestamp);
  calin_cdts_data->set_camera_timestamp(cdts_data->camera_timestamp);
  calin_cdts_data->set_trigger_type(cdts_data->trigger_type);
  calin_cdts_data->set_white_rabbit_status(cdts_data->white_rabbit_status);
  calin_cdts_data->set_arbitrary_information(cdts_data->arbitrary_information);
}

void calin::iact_data::zfits_data_source::
decode_tib_data(calin::ix::iact_data::telescope_event::TIBData* calin_tib_data,
  const DataModel::AnyArray& cta_array)
{
  // No of bits            Data
  // 95 - 64 (32 bits)     Event Counter
  // 63 - 48 (16 bits)     PPS Counter
  // 47 - 24 (24 bits)     10 MHz Counter
  // 23 - 17 (7 bits)      Zeros
  // 16 - 8 (9 bits)       Stereo Pattern
  // 7 - 0 (8 bits)        Trigger Type
  //
  // Inside the trigger type byte, the meaning of each bit is:
  //
  // Bit   Meaning
  // 0     Mono
  // 1     Stereo
  // 2     Calibration
  // 3     Single Photo-electron
  // 4     Auxiliary trigger from UCTS
  // 5     Pedestal
  // 6     Slow Control
  // 7     Busy

  struct TIBMessageData {
    uint32_t event_counter;
    uint16_t pps_counter;
    uint16_t clock_counter_lo16;
    uint8_t  clock_counter_hi8;
    uint16_t stereo_pattern;
    uint8_t  trigger_type;
  } __attribute__((__packed__));

  const auto& cta_tib_data = cta_array.data();
#if TEST_ANYARRAY_TYPES
  if(cta_array.type() != DataModel::AnyArray::U8)
    throw std::runtime_error("TIB type not U8");
#endif
  if(cta_tib_data.size() != sizeof(TIBMessageData))
    throw std::runtime_error("TIB data array not expected size");
  const auto* tib_data =
    reinterpret_cast<const TIBMessageData*>(&cta_tib_data.front());

  calin_tib_data->set_event_counter(tib_data->event_counter);
  calin_tib_data->set_pps_counter(tib_data->pps_counter);
  calin_tib_data->set_clock_counter(tib_data->clock_counter_lo16
    + (tib_data->clock_counter_hi8<<16) );
  calin_tib_data->set_stereo_pattern(tib_data->stereo_pattern&0x0001FFFF);
  calin_tib_data->set_mono_trigger(tib_data->trigger_type & 0x01);
  calin_tib_data->set_stereo_trigger(tib_data->trigger_type & 0x02);
  calin_tib_data->set_external_calibration_trigger(tib_data->trigger_type & 0x04);
  calin_tib_data->set_internal_calibration_trigger(tib_data->trigger_type & 0x08);
  calin_tib_data->set_ucts_aux_trigger(tib_data->trigger_type & 0x10);
  calin_tib_data->set_pedestal_trigger(tib_data->trigger_type & 0x20);
  calin_tib_data->set_slow_control_trigger(tib_data->trigger_type & 0x40);
  calin_tib_data->set_busy_trigger(tib_data->trigger_type & 0x80);
  calin_tib_data->set_spare_bits(tib_data->stereo_pattern>>9);
}
