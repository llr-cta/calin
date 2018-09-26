/*

   calin/iact_data/actl_event_decoder.cpp -- Stephen Fegan -- 2018-09-21

   A decoder of ACTL event types

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
#include <memory>
#include <cctype>

#include <util/log.hpp>
#include <util/file.hpp>
#include <iact_data/actl_event_decoder.hpp>

using namespace calin::iact_data::actl_event_decoder;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::util::log;
using calin::util::file::is_file;
using calin::util::file::is_readable;
using calin::util::file::expand_filename;

#include <ProtobufIFits.h>
#include <CoreMessages.pb.h>

// ----------------------------------------------------------------------------
//
// Utility decoding functions
//
// ----------------------------------------------------------------------------

void calin::iact_data::actl_event_decoder::
decode_cdts_data(calin::ix::iact_data::telescope_event::CDTSData* calin_cdts_data,
  const DataModel::AnyArray& cta_array)
{
  // Reference : https://forge.in2p3.fr/projects/cta/repository/entry/ACTL/ExternalDevicesCommunication/trunk/TiCkSdecode/ticks_decode.c
  struct CDTSMessageData_V0 {
    uint32_t event_counter;
    uint32_t pps_counter;
    uint32_t clock_counter;
    uint64_t ucts_timestamp;
    uint64_t camera_timestamp;
    uint8_t trigger_type;
    uint8_t white_rabbit_status;
    uint8_t arbitrary_information;
  } __attribute__((__packed__));

  struct CDTSMessageData_V1 {
    uint32_t event_counter;
    uint32_t busy_counter;
    uint32_t pps_counter;
    uint32_t clock_counter;
    uint64_t ucts_timestamp;
    uint8_t trigger_type;
    uint8_t white_rabbit_reset_busy_status;
    uint16_t arbitrary_information;
    uint8_t num_in_bunch;
  } __attribute__((__packed__));

  const auto& cta_cdts_data = cta_array.data();
#if TEST_ANYARRAY_TYPES
  if(cta_cdts_data.type() != DataModel::AnyArray::U32)
    throw std::runtime_error("CDTS counters type not U32");
#endif
  if(cta_cdts_data.size() == sizeof(CDTSMessageData_V1)) {
    const auto* cdts_data =
      reinterpret_cast<const CDTSMessageData_V1*>(&cta_cdts_data.front());

    calin_cdts_data->set_event_counter(cdts_data->event_counter);
    calin_cdts_data->set_busy_counter(cdts_data->busy_counter);
    calin_cdts_data->set_pps_counter(cdts_data->pps_counter);
    calin_cdts_data->set_clock_counter(cdts_data->clock_counter);
    calin_cdts_data->set_ucts_timestamp(cdts_data->ucts_timestamp);
    calin_cdts_data->set_camera_timestamp(0);
    calin_cdts_data->set_trigger_type(cdts_data->trigger_type);
    calin_cdts_data->set_white_rabbit_status(cdts_data->white_rabbit_reset_busy_status);
    calin_cdts_data->set_arbitrary_information(cdts_data->arbitrary_information);
    calin_cdts_data->set_num_in_bunch(cdts_data->num_in_bunch);
    calin_cdts_data->set_version(1);
  } else if(cta_cdts_data.size() == sizeof(CDTSMessageData_V0)) {
    const auto* cdts_data =
      reinterpret_cast<const CDTSMessageData_V0*>(&cta_cdts_data.front());

    calin_cdts_data->set_event_counter(cdts_data->event_counter);
    calin_cdts_data->set_busy_counter(0);
    calin_cdts_data->set_pps_counter(cdts_data->pps_counter);
    calin_cdts_data->set_clock_counter(cdts_data->clock_counter);
    calin_cdts_data->set_ucts_timestamp(cdts_data->ucts_timestamp);
    calin_cdts_data->set_camera_timestamp(cdts_data->camera_timestamp);
    calin_cdts_data->set_trigger_type(cdts_data->trigger_type);
    calin_cdts_data->set_white_rabbit_status(cdts_data->white_rabbit_status);
    calin_cdts_data->set_arbitrary_information(cdts_data->arbitrary_information);
    calin_cdts_data->set_num_in_bunch(0);
    calin_cdts_data->set_version(0);
  } else {
    throw std::runtime_error("CDTS data array not expected size");
  }

}

void calin::iact_data::actl_event_decoder::
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

/*

              LLLLLLLLLLL                       000000000
              L:::::::::L                     00:::::::::00
              L:::::::::L                   00:::::::::::::00
              LL:::::::LL                  0:::::::000:::::::0
                L:::::L                    0::::::0   0::::::0
                L:::::L                    0:::::0     0:::::0
                L:::::L                    0:::::0     0:::::0
                L:::::L                    0:::::0 000 0:::::0
                L:::::L                    0:::::0 000 0:::::0
                L:::::L                    0:::::0     0:::::0
                L:::::L                    0:::::0     0:::::0
                L:::::L         LLLLLL     0::::::0   0::::::0
              LL:::::::LLLLLLLLL:::::L     0:::::::000:::::::0
              L::::::::::::::::::::::L      00:::::::::::::00
              L::::::::::::::::::::::L        00:::::::::00
              LLLLLLLLLLLLLLLLLLLLLLLL          000000000

*/

#include <L0.pb.h>

ACTL_L0_CameraEventDecoder::~ACTL_L0_CameraEventDecoder()
{
  // nothing to see here
}

// =============================================================================
// DecodedACTL_L0_CameraEventDataSource - extract ACTL events from an
// ACTL_L0_CameraEventDataSource and decode them using a supplied decoder
// =============================================================================

DecodedACTL_L0_CameraEventDataSource::DecodedACTL_L0_CameraEventDataSource(
    calin::iact_data::zfits_actl_data_source::ACTL_L0_CameraEventDataSource* actl_src,
    ACTL_L0_CameraEventDecoder* decoder, bool adopt_actl_src, bool adopt_decoder):
  calin::iact_data::telescope_data_source::TelescopeDataSource(),
  decoder_(decoder), adopt_decoder_(adopt_decoder),
  actl_src_(actl_src), adopt_actl_src_(adopt_actl_src)
{
  // nothing to see here
}

DecodedACTL_L0_CameraEventDataSource::~DecodedACTL_L0_CameraEventDataSource()
{
  if(adopt_decoder_)delete decoder_;
  if(adopt_actl_src_)delete actl_src_;
}

calin::ix::iact_data::telescope_event::TelescopeEvent*
DecodedACTL_L0_CameraEventDataSource::get_next(
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

DecodedConstACTL_L0_CameraEventDataSource::DecodedConstACTL_L0_CameraEventDataSource(
    calin::iact_data::zfits_actl_data_source::ConstACTL_L0_CameraEventDataSource* actl_src,
    calin::iact_data::zfits_actl_data_source::ConstACTL_L0_CameraEventDataSink* actl_sink,
    ACTL_L0_CameraEventDecoder* decoder,
    bool adopt_actl_src, bool adopt_actl_sink, bool adopt_decoder):
  calin::iact_data::telescope_data_source::TelescopeDataSource(),
  decoder_(decoder), adopt_decoder_(adopt_decoder),
  actl_src_(actl_src), adopt_actl_src_(adopt_actl_src),
  actl_sink_(actl_sink), adopt_actl_sink_(adopt_actl_sink)
{
  // nothing to see here
}

DecodedConstACTL_L0_CameraEventDataSource::DecodedConstACTL_L0_CameraEventDataSource(
    calin::iact_data::zfits_actl_data_source::ConstACTL_L0_CameraEventDataSource* actl_src,
    ACTL_L0_CameraEventDecoder* decoder,
    bool adopt_actl_src, bool adopt_decoder):
  calin::iact_data::telescope_data_source::TelescopeDataSource(),
  decoder_(decoder), adopt_decoder_(adopt_decoder),
  actl_src_(actl_src), adopt_actl_src_(adopt_actl_src),
  actl_sink_(nullptr), adopt_actl_sink_(false)
{
  // nothing to see here
}

DecodedConstACTL_L0_CameraEventDataSource::~DecodedConstACTL_L0_CameraEventDataSource()
{
  if(adopt_decoder_)delete decoder_;
  if(adopt_actl_src_)delete actl_src_;
  if(adopt_actl_sink_)delete actl_sink_;
}

calin::ix::iact_data::telescope_event::TelescopeEvent*
DecodedConstACTL_L0_CameraEventDataSource::get_next(
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
    if(actl_sink_)actl_sink_->put_next(cta_event, seq_index_out, nullptr, true);
    else delete cta_event;
    delete delete_arena;
    throw std::runtime_error("Could not allocate telescope event");
  }
  if(!decoder_->decode(event, cta_event))
  {
    if(actl_sink_)actl_sink_->put_next(cta_event, seq_index_out, nullptr, true);
    else delete cta_event;
    delete delete_arena;
    delete delete_event;
    throw std::runtime_error("Could not decode ACTL event");
  }
  event->set_source_event_index(seq_index_out);
  if(actl_sink_)actl_sink_->put_next(cta_event, seq_index_out, nullptr, true);
  else delete cta_event;
  return event;
}

DecodedConstACTL_L0_CameraEventDataSourceFactory::DecodedConstACTL_L0_CameraEventDataSourceFactory(
    calin::io::data_source::BidirectionalBufferedDataSourcePump<
      const DataModel::CameraEvent>* pump, ACTL_L0_CameraEventDecoder* decoder,
    bool adopt_pump, bool adopt_decoder):
  calin::iact_data::telescope_data_source::TelescopeDataSourceFactory(),
  decoder_(decoder), adopt_decoder_(adopt_decoder),
  pump_(pump), adopt_pump_(adopt_pump)
{
  // nothing to see here
}

DecodedConstACTL_L0_CameraEventDataSourceFactory::~DecodedConstACTL_L0_CameraEventDataSourceFactory()
{
  if(adopt_pump_)delete pump_;
  if(adopt_decoder_)delete decoder_;
}

DecodedConstACTL_L0_CameraEventDataSource* DecodedConstACTL_L0_CameraEventDataSourceFactory::new_data_source()
{
  return new DecodedConstACTL_L0_CameraEventDataSource(
    pump_->new_data_source(), pump_->new_data_sink(), decoder_,
    /* adopt_actl_src= */ true, /* adopt_actl_sink = */ true);
}

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

#include <R1.pb.h>

ACTL_R1_CameraEventDecoder::~ACTL_R1_CameraEventDecoder()
{
  // nothing to see here
}

// =============================================================================
// DecodedACTL_R1_CameraEventDataSource - extract ACTL events from an
// ACTL_R1_CameraEventDataSource and decode them using a supplied decoder
// =============================================================================

DecodedACTL_R1_CameraEventDataSource::DecodedACTL_R1_CameraEventDataSource(
    calin::iact_data::zfits_actl_data_source::ACTL_R1_CameraEventDataSource* actl_src,
    ACTL_R1_CameraEventDecoder* decoder, bool adopt_actl_src, bool adopt_decoder):
  calin::iact_data::telescope_data_source::TelescopeDataSource(),
  decoder_(decoder), adopt_decoder_(adopt_decoder),
  actl_src_(actl_src), adopt_actl_src_(adopt_actl_src)
{
  // nothing to see here
}

DecodedACTL_R1_CameraEventDataSource::~DecodedACTL_R1_CameraEventDataSource()
{
  if(adopt_decoder_)delete decoder_;
  if(adopt_actl_src_)delete actl_src_;
}

calin::ix::iact_data::telescope_event::TelescopeEvent*
DecodedACTL_R1_CameraEventDataSource::get_next(
  uint64_t& seq_index_out, google::protobuf::Arena** arena)
{
  const R1::CameraEvent* cta_event = actl_src_->get_next(seq_index_out);
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

DecodedConstACTL_R1_CameraEventDataSource::DecodedConstACTL_R1_CameraEventDataSource(
    calin::iact_data::zfits_actl_data_source::ConstACTL_R1_CameraEventDataSource* actl_src,
    calin::iact_data::zfits_actl_data_source::ConstACTL_R1_CameraEventDataSink* actl_sink,
    ACTL_R1_CameraEventDecoder* decoder,
    bool adopt_actl_src, bool adopt_actl_sink, bool adopt_decoder):
  calin::iact_data::telescope_data_source::TelescopeDataSource(),
  decoder_(decoder), adopt_decoder_(adopt_decoder),
  actl_src_(actl_src), adopt_actl_src_(adopt_actl_src),
  actl_sink_(actl_sink), adopt_actl_sink_(adopt_actl_sink)
{
  // nothing to see here
}

DecodedConstACTL_R1_CameraEventDataSource::DecodedConstACTL_R1_CameraEventDataSource(
    calin::iact_data::zfits_actl_data_source::ConstACTL_R1_CameraEventDataSource* actl_src,
    ACTL_R1_CameraEventDecoder* decoder,
    bool adopt_actl_src, bool adopt_decoder):
  calin::iact_data::telescope_data_source::TelescopeDataSource(),
  decoder_(decoder), adopt_decoder_(adopt_decoder),
  actl_src_(actl_src), adopt_actl_src_(adopt_actl_src),
  actl_sink_(nullptr), adopt_actl_sink_(false)
{
  // nothing to see here
}

DecodedConstACTL_R1_CameraEventDataSource::~DecodedConstACTL_R1_CameraEventDataSource()
{
  if(adopt_decoder_)delete decoder_;
  if(adopt_actl_src_)delete actl_src_;
  if(adopt_actl_sink_)delete actl_sink_;
}

calin::ix::iact_data::telescope_event::TelescopeEvent*
DecodedConstACTL_R1_CameraEventDataSource::get_next(
  uint64_t& seq_index_out, google::protobuf::Arena** arena)
{
  const R1::CameraEvent* cta_event = actl_src_->get_next(seq_index_out);
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
    if(actl_sink_)actl_sink_->put_next(cta_event, seq_index_out, nullptr, true);
    else delete cta_event;
    delete delete_arena;
    throw std::runtime_error("Could not allocate telescope event");
  }
  if(!decoder_->decode(event, cta_event))
  {
    if(actl_sink_)actl_sink_->put_next(cta_event, seq_index_out, nullptr, true);
    else delete cta_event;
    delete delete_arena;
    delete delete_event;
    throw std::runtime_error("Could not decode ACTL event");
  }
  event->set_source_event_index(seq_index_out);
  if(actl_sink_)actl_sink_->put_next(cta_event, seq_index_out, nullptr, true);
  else delete cta_event;
  return event;
}

DecodedConstACTL_R1_CameraEventDataSourceFactory::DecodedConstACTL_R1_CameraEventDataSourceFactory(
    calin::io::data_source::BidirectionalBufferedDataSourcePump<
      const R1::CameraEvent>* pump, ACTL_R1_CameraEventDecoder* decoder,
    bool adopt_pump, bool adopt_decoder):
  calin::iact_data::telescope_data_source::TelescopeDataSourceFactory(),
  decoder_(decoder), adopt_decoder_(adopt_decoder),
  pump_(pump), adopt_pump_(adopt_pump)
{
  // nothing to see here
}

DecodedConstACTL_R1_CameraEventDataSourceFactory::~DecodedConstACTL_R1_CameraEventDataSourceFactory()
{
  if(adopt_pump_)delete pump_;
  if(adopt_decoder_)delete decoder_;
}

DecodedConstACTL_R1_CameraEventDataSource* DecodedConstACTL_R1_CameraEventDataSourceFactory::new_data_source()
{
  return new DecodedConstACTL_R1_CameraEventDataSource(
    pump_->new_data_source(), pump_->new_data_sink(), decoder_,
    /* adopt_actl_src= */ true, /* adopt_actl_sink = */ true);
}
