/*

   calin/iact_data/acada_event_decoder.cpp -- Stephen Fegan -- 2024-09-05

   Base for all decoders of ACADA event types

   Copyright 2024, Stephen Fegan <sfegan@llr.in2p3.fr>
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
#include <memory>
#include <cctype>

#include <util/log.hpp>
#include <util/file.hpp>
#include <iact_data/acada_event_decoder.hpp>

using namespace calin::iact_data::acada_data_source;
using namespace calin::iact_data::acada_event_decoder;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::util::log;
using calin::util::file::is_file;
using calin::util::file::is_readable;
using calin::util::file::expand_filename;

template<typename MessageSet>
ACADACameraEventDecoder<MessageSet>::
~ACADACameraEventDecoder()
{
  // nothing to see here
}

#if 0
// =============================================================================
// DecodedACADACameraEventDataSource - extract ACTL events from an
// ACADACameraEventDataSource and decode them using a supplied decoder
// =============================================================================

template<typename EventMessage, typename HeaderMessage>
DecodedACADACameraEventDataSource<EventMessage,HeaderMessage>::
DecodedACADACameraEventDataSource(
    ACADACameraEventDataSource* acada_src, 
    ACADACameraEventDecoder* decoder,
    bool adopt_acada_src, bool adopt_decoder):
  calin::iact_data::telescope_data_source::TelescopeDataSource(),
  decoder_(decoder), adopt_decoder_(adopt_decoder),
  acada_src_(acada_src), adopt_acada_src_(adopt_acada_src)
{
  // nothing to see here
}

template<typename EventMessage, typename HeaderMessage>
DecodedACADACameraEventDataSource<EventMessage,HeaderMessage>::
~DecodedACADACameraEventDataSource()
{
  if(adopt_decoder_)delete decoder_;
  if(adopt_acada_src_)delete acada_src_;
}

template<typename EventMessage, typename HeaderMessage>
calin::ix::iact_data::telescope_event::TelescopeEvent* 
DecodedACADACameraEventDataSource<EventMessage,HeaderMessage>::get_next(
  uint64_t& seq_index_out, google::protobuf::Arena** arena)
{
  const EventMessage* cta_event = this->borrow_next_acada_event(seq_index_out);
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
    this->release_borrowed_acada_event(cta_event);
    delete delete_arena;
    throw std::runtime_error("Could not allocate telescope event");
  }
  if(!decoder_->decode(event, cta_event))
  {
    this->release_borrowed_acada_event(cta_event);
    delete delete_arena;
    delete delete_event;
    throw std::runtime_error("Could not decode ACTL event");
  }
  event->set_source_event_index(seq_index_out);
  this->release_borrowed_acada_event(cta_event);
  return event;
}

template<typename EventMessage, typename HeaderMessage>
const EventMessage* DecodedACADACameraEventDataSource<EventMessage,HeaderMessage>::
borrow_next_acada_event(uint64_t& seq_index_out)
{
  return acada_src_->borrow_next_event(seq_index_out);
}

template<typename EventMessage, typename HeaderMessage>
void DecodedACADACameraEventDataSource<EventMessage,HeaderMessage>::
release_borrowed_acada_event(const EventMessage* event)
{
  acada_src_->release_borrowed_event(event);
}

// =============================================================================
// DecodedACADACameraEventDataSourceWithRunConfig - extract ACTA events
// from an ACADACameraEventDataSourceWithRunHeader and decode them
// using a supplied decoder
// =============================================================================

template<typename EventMessage, typename HeaderMessage>
DecodedACADACameraEventDataSourceWithRunConfig<EventMessage,HeaderMessage>::
DecodedACADACameraEventDataSourceWithRunConfig(  
    ACADACameraEventDataSourceWithRunHeader* acada_src,
    ACADACameraEventDecoder* decoder,
    bool adopt_acada_src, bool adopt_decoder):
  DecodedACADACameraEventDataSource<EventMessage,HeaderMessage>(acada_src, 
    decoder, adopt_acada_src, adopt_decoder)
{
  // nothing to see here
}

template<typename EventMessage, typename HeaderMessage>
DecodedACADACameraEventDataSourceWithRunConfig<EventMessage,HeaderMessage>::
~DecodedACADACameraEventDataSourceWithRunConfig()
{
  if(saved_event_)acada_src_->release_borrowed_event(saved_event_);
  if(run_config_)delete run_config_;
}

template<typename EventMessage, typename HeaderMessage>
calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration*
DecodedACADACameraEventDataSourceWithRunConfig<EventMessage,HeaderMessage>::
get_run_configuration()
{
  ensure_run_config();
  return new TelescopeRunConfiguration(*run_config_);
}

template<typename EventMessage, typename HeaderMessage>
const EventMessage* DecodedACADACameraEventDataSourceWithRunConfig<EventMessage,HeaderMessage>::
borrow_next_acada_event(uint64_t& seq_index_out)
{
  if(saved_event_) {
    const EventMessage* event = saved_event_;
    seq_index_out = saved_seq_index_;
    saved_event_ = nullptr;
    saved_seq_index_ = 0;
    return event;
  }
  
  return DecodedACADACameraEventDataSource<EventMessage,HeaderMessage>::
    borrow_next_acada_event(seq_index_out);
}

template<typename EventMessage, typename HeaderMessage>
void DecodedACADACameraEventDataSourceWithRunConfig<EventMessage,HeaderMessage>::
ensure_run_config()
{
  if(!run_config_) {
    saved_event_ = this->borrow_next_acada_event(saved_seq_index_);
    HeaderMessage* run_header = acada_src_->get_run_header();
    run_config_ = new TelescopeRunConfiguration;
    this->decoder_->decode_run_config(run_config_, run_header, saved_event_);
    delete run_header;
  }
}
#endif

namespace calin { namespace iact_data { namespace acada_event_decoder {

template class ACADACameraEventDecoder<ACADA_MessageSet_L0>;
// template class DecodedACADACameraEventDataSource<ACADA_EventMessage_L0, ACADA_HeaderMessage_L0>;
// template class DecodedACADACameraEventDataSourceWithRunConfig<ACADA_EventMessage_L0, ACADA_HeaderMessage_L0>;

template class ACADACameraEventDecoder<ACADA_MessageSet_R1v0>;
// template class DecodedACADACameraEventDataSource<ACADA_EventMessage_R1v0, ACADA_HeaderMessage_R1v0>;
// template class DecodedACADACameraEventDataSourceWithRunConfig<ACADA_EventMessage_R1v0, ACADA_HeaderMessage_R1v0>;

template class ACADACameraEventDecoder<ACADA_MessageSet_R1v1>;

} } } // namespace calin::iact_data::acada_event_decoder

// =============================================================================
//
// Utility decoding functions
//
// =============================================================================

void calin::iact_data::acada_event_decoder::
decode_cdts_data(calin::ix::iact_data::telescope_event::CDTSData* calin_cdts_data,
  const AnyArray& cta_array)
{
  // Reference : https://forge.in2p3.fr/projects/cta/repository/entry/ACTL/ExternalDevicesCommunication/trunk/TiCkSdecode/ticks_decode.c
  struct CDTSMessageData_V0 { // 31 bytes
    /*  4 */ uint32_t event_counter;
    /*  8 */ uint32_t pps_counter;
    /* 12 */ uint32_t clock_counter;
    /* 20 */ uint64_t ucts_timestamp;
    /* 28 */ uint64_t camera_timestamp;
    /* 29 */ uint8_t trigger_type;
    /* 30 */ uint8_t white_rabbit_status;
    /* 31 */ uint8_t arbitrary_information;
  } __attribute__((__packed__));

  struct CDTSMessageData_V1 { // 29 bytes
    /*  4 */ uint32_t event_counter;
    /*  8 */ uint32_t busy_counter;
    /* 12 */ uint32_t pps_counter;
    /* 16 */ uint32_t clock_counter;
    /* 24 */ uint64_t ucts_timestamp;
    /* 25 */ uint8_t trigger_type;
    /* 26 */ uint8_t white_rabbit_reset_busy_status;
    /* 28 */ uint16_t arbitrary_information;
    /* 29 */ uint8_t num_in_bunch;
  } __attribute__((__packed__));

  struct CDTSMessageData_V2 { // 28 bytes
    /*  4 */ uint32_t event_counter;
    /*  8 */ uint32_t busy_counter;
    /* 12 */ uint32_t pps_counter;
    /* 16 */ uint32_t clock_counter;
    /* 24 */ uint64_t ucts_timestamp;
    /* 25 */ uint8_t trigger_type; // For TIB, this is the first 7 bits of SPI
    /* 26 */ uint8_t white_rabbit_reset_busy_status; // A few status flags here, see below
    /* 27 */ uint8_t stereo_pattern; // For TIB, this is the next 8 bits of the SPI string
    /* 28 */ uint8_t num_in_bunch; // This information only needed for debugging
  } __attribute__((__packed__));

  // V3 from : https://forge.in2p3.fr/projects/cta/repository/revisions/37071/entry/ACTL/ExternalDevicesCommunication/trunk/TiCkS/TiCkS_decode/ticks_decode.h

  struct CDTSMessageData_V3 { // 36 bytes
    /*  8 */ uint64_t ucts_timestamp;
    /* 12 */ uint32_t ucts_address;
    /* 16 */ uint32_t event_counter;
    /* 20 */ uint32_t busy_counter;
    /* 24 */ uint32_t pps_counter;
    /* 28 */ uint32_t clock_counter;
    /* 29 */ uint8_t trigger_type; // For TIB, this is the first 7 bits of SPI
    /* 30 */ uint8_t white_rabbit_reset_busy_status; // A few status flags here, see below
    /* 31 */ uint8_t stereo_pattern; // For TIB, this is the next 8 bits of the SPI string
    /* 32 */ uint8_t num_in_bunch; // This information only needed for debugging
    /* 36 */ uint32_t cdts_version; // Should be x.y.z where x is 2bytes and y, z are a byte each
  } __attribute__((__packed__));

  /* whiteRabbitResetBusyStatus bits defined as:

     Bit Meaning :
     0 White Rabbit Status (PLL Lock)
     1 Reset counter acknowledge
     6 Busy bit from SPI (should be the same as the busy flag, bit 7)
     7 Busy flag (defined by UCTS depending on which channel the T-S came in on)
  */

  /* Note: "triggerType" and "StereoPattern" from TIB SPI defined as follows:
     From document TIB User manual v4.6: https://portal.cta-observatory.org/WG/mst/CAM/NeCTArCam/NectarCamDoc/Shared%20Documents/07.%20Design/7.4%20Trigger%20and%20Clock/7.4.2%20TIB/TIB%20User%20Manual.pdf
     Table 9: TIB to UCTS SPI trigger type message structure

     Bit Meaning :
     0 Mono
     1 Stereo
     2 Calib
     3 Sphe
     4 Softrig Trigger Pattern
     5 Pedestal
     6 Slow
     7 Local
     8 NB1
     9 NB2
     10 NB3
     11 NB4 Stereo Pattern
     12 NB5
     13 NB6
     14 NB7  <- Replaced by Muon Candidate Flag MCF https://indico.cta-observatory.org/event/1293/attachments/23531/33935/NectarCAM-MuonCandidateFlag_210423_v3.pdf
     15 Busy Busy  (put the SPI busy into the whiteRabbitResetBusyStatus)
   */

  const auto& cta_cdts_data = cta_array.data();
#if TEST_ANYARRAY_TYPES
  if(cta_cdts_data.type() != AnyArray::U32)
    throw std::runtime_error("CDTS counters type not U32");
#endif
  if(cta_cdts_data.size() == sizeof(CDTSMessageData_V3)) {
    const auto* cdts_data =
      reinterpret_cast<const CDTSMessageData_V3*>(&cta_cdts_data.front());

    calin_cdts_data->set_event_counter(cdts_data->event_counter);
    calin_cdts_data->set_busy_counter(cdts_data->busy_counter);
    calin_cdts_data->set_pps_counter(cdts_data->pps_counter);
    calin_cdts_data->set_clock_counter(cdts_data->clock_counter);
    calin_cdts_data->set_ucts_timestamp(cdts_data->ucts_timestamp);
    calin_cdts_data->set_camera_timestamp(0);
    calin_cdts_data->set_trigger_type(cdts_data->trigger_type);
    calin_cdts_data->set_white_rabbit_status(cdts_data->white_rabbit_reset_busy_status);
    calin_cdts_data->set_stereo_pattern(cdts_data->stereo_pattern);
    calin_cdts_data->set_arbitrary_information(0);
    calin_cdts_data->set_num_in_bunch(cdts_data->num_in_bunch);
    calin_cdts_data->set_ucts_address(cdts_data->ucts_address);
    calin_cdts_data->set_cdts_version(cdts_data->cdts_version);
    calin_cdts_data->set_version(3);

    calin_cdts_data->set_mono_trigger(cdts_data->trigger_type & 0x01);
    calin_cdts_data->set_stereo_trigger(cdts_data->trigger_type & 0x02);
    calin_cdts_data->set_external_calibration_trigger(cdts_data->trigger_type & 0x04);
    calin_cdts_data->set_internal_calibration_trigger(cdts_data->trigger_type & 0x08);
    calin_cdts_data->set_ucts_aux_trigger(cdts_data->trigger_type & 0x10);
    calin_cdts_data->set_pedestal_trigger(cdts_data->trigger_type & 0x20);
    calin_cdts_data->set_slow_control_trigger(cdts_data->trigger_type & 0x40);
    calin_cdts_data->set_local_trigger(cdts_data->stereo_pattern & 0x01);
    calin_cdts_data->set_muon_candidate(cdts_data->stereo_pattern & 0x80);
    // calin_cdts_data->set_busy_trigger(cdts_data->stereo_pattern & 0x80);
  } else if(cta_cdts_data.size() == sizeof(CDTSMessageData_V2)) {
    const auto* cdts_data =
      reinterpret_cast<const CDTSMessageData_V2*>(&cta_cdts_data.front());

    calin_cdts_data->set_event_counter(cdts_data->event_counter);
    calin_cdts_data->set_busy_counter(cdts_data->busy_counter);
    calin_cdts_data->set_pps_counter(cdts_data->pps_counter);
    calin_cdts_data->set_clock_counter(cdts_data->clock_counter);
    calin_cdts_data->set_ucts_timestamp(cdts_data->ucts_timestamp);
    calin_cdts_data->set_camera_timestamp(0);
    calin_cdts_data->set_trigger_type(cdts_data->trigger_type);
    calin_cdts_data->set_white_rabbit_status(cdts_data->white_rabbit_reset_busy_status);
    calin_cdts_data->set_stereo_pattern(cdts_data->stereo_pattern);
    calin_cdts_data->set_arbitrary_information(0);
    calin_cdts_data->set_num_in_bunch(cdts_data->num_in_bunch);
    calin_cdts_data->set_ucts_address(0);
    calin_cdts_data->set_cdts_version(0);
    calin_cdts_data->set_version(2);

    calin_cdts_data->set_mono_trigger(cdts_data->trigger_type & 0x01);
    calin_cdts_data->set_stereo_trigger(cdts_data->trigger_type & 0x02);
    calin_cdts_data->set_external_calibration_trigger(cdts_data->trigger_type & 0x04);
    calin_cdts_data->set_internal_calibration_trigger(cdts_data->trigger_type & 0x08);
    calin_cdts_data->set_ucts_aux_trigger(cdts_data->trigger_type & 0x10);
    calin_cdts_data->set_pedestal_trigger(cdts_data->trigger_type & 0x20);
    calin_cdts_data->set_slow_control_trigger(cdts_data->trigger_type & 0x40);
    calin_cdts_data->set_local_trigger(cdts_data->stereo_pattern & 0x01);
    calin_cdts_data->set_muon_candidate(cdts_data->stereo_pattern & 0x80);
  } else if(cta_cdts_data.size() == sizeof(CDTSMessageData_V1)) {
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
    calin_cdts_data->set_stereo_pattern(0);
    calin_cdts_data->set_arbitrary_information(cdts_data->arbitrary_information);
    calin_cdts_data->set_num_in_bunch(cdts_data->num_in_bunch);
    calin_cdts_data->set_ucts_address(0);
    calin_cdts_data->set_cdts_version(0);
    calin_cdts_data->set_version(1);

    calin_cdts_data->set_mono_trigger(cdts_data->trigger_type & 0x01);
    calin_cdts_data->set_stereo_trigger(cdts_data->trigger_type & 0x02);
    calin_cdts_data->set_external_calibration_trigger(cdts_data->trigger_type & 0x04);
    calin_cdts_data->set_internal_calibration_trigger(cdts_data->trigger_type & 0x08);
    calin_cdts_data->set_ucts_aux_trigger(cdts_data->trigger_type & 0x10);
    calin_cdts_data->set_pedestal_trigger(cdts_data->trigger_type & 0x20);
    calin_cdts_data->set_slow_control_trigger(cdts_data->trigger_type & 0x40);
    calin_cdts_data->set_local_trigger(cdts_data->trigger_type & 0x80);
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
    calin_cdts_data->set_stereo_pattern(0);
    calin_cdts_data->set_arbitrary_information(cdts_data->arbitrary_information);
    calin_cdts_data->set_num_in_bunch(0);
    calin_cdts_data->set_ucts_address(0);
    calin_cdts_data->set_cdts_version(0);
    calin_cdts_data->set_version(0);

    calin_cdts_data->set_mono_trigger(cdts_data->trigger_type & 0x01);
    calin_cdts_data->set_stereo_trigger(cdts_data->trigger_type & 0x02);
    calin_cdts_data->set_external_calibration_trigger(cdts_data->trigger_type & 0x04);
    calin_cdts_data->set_internal_calibration_trigger(cdts_data->trigger_type & 0x08);
    calin_cdts_data->set_ucts_aux_trigger(cdts_data->trigger_type & 0x10);
    calin_cdts_data->set_pedestal_trigger(cdts_data->trigger_type & 0x20);
    calin_cdts_data->set_slow_control_trigger(cdts_data->trigger_type & 0x40);
    calin_cdts_data->set_local_trigger(cdts_data->trigger_type & 0x80);
  } else {
    throw std::runtime_error("CDTS data array not expected size");
  }
}

void calin::iact_data::acada_event_decoder::
decode_tib_data(calin::ix::iact_data::telescope_event::TIBData* calin_tib_data,
  const AnyArray& cta_array)
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
  if(cta_array.type() != AnyArray::U8)
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
  calin_tib_data->set_stereo_pattern(tib_data->stereo_pattern&0x1FF);
  calin_tib_data->set_trigger_type(tib_data->trigger_type);
  calin_tib_data->set_spare_bits(tib_data->stereo_pattern>>9);

  calin_tib_data->set_mono_trigger(tib_data->trigger_type & 0x01);
  calin_tib_data->set_stereo_trigger(tib_data->trigger_type & 0x02);
  calin_tib_data->set_external_calibration_trigger(tib_data->trigger_type & 0x04);
  calin_tib_data->set_internal_calibration_trigger(tib_data->trigger_type & 0x08);
  calin_tib_data->set_ucts_aux_trigger(tib_data->trigger_type & 0x10);
  calin_tib_data->set_pedestal_trigger(tib_data->trigger_type & 0x20);
  calin_tib_data->set_slow_control_trigger(tib_data->trigger_type & 0x40);
  calin_tib_data->set_busy_trigger(tib_data->trigger_type & 0x80);
}

calin::ix::iact_data::telescope_event::TriggerType
calin::iact_data::acada_event_decoder::determine_trigger_type(
  const calin::ix::iact_data::telescope_event::TIBData* calin_tib_data,
  const calin::ix::iact_data::telescope_event::CDTSData* calin_cdts_data)
{
  if(calin_tib_data) {
    switch(calin_tib_data->trigger_type()) {
    case 0x01:
      return calin::ix::iact_data::telescope_event::TRIGGER_PHYSICS;
    case 0x02:
      return calin::ix::iact_data::telescope_event::TRIGGER_FORCED_BY_ARRAY;
    case 0x04:
    case 0x05: // allow either pure external or external/physics
      return calin::ix::iact_data::telescope_event::TRIGGER_EXTERNAL_FLASHER;
    case 0x08:
      return calin::ix::iact_data::telescope_event::TRIGGER_INTERNAL_FLASHER;
    case 0x10:
      return calin::ix::iact_data::telescope_event::TRIGGER_UCTS_AUX;
    case 0x20:
      return calin::ix::iact_data::telescope_event::TRIGGER_PEDESTAL;
    case 0x40:
      return calin::ix::iact_data::telescope_event::TRIGGER_SOFTWARE;
    case 0x80:
      return calin::ix::iact_data::telescope_event::TRIGGER_UNKNOWN;
    default:
      return calin::ix::iact_data::telescope_event::TRIGGER_MULTIPLE;
    }
  } else if(calin_cdts_data) {
    switch(calin_cdts_data->trigger_type()) {
    case 0x01:
      return calin::ix::iact_data::telescope_event::TRIGGER_PHYSICS;
    case 0x02:
      return calin::ix::iact_data::telescope_event::TRIGGER_FORCED_BY_ARRAY;
    case 0x04:
    case 0x05: // allow either pure external or external/physics
      return calin::ix::iact_data::telescope_event::TRIGGER_EXTERNAL_FLASHER;
    case 0x08:
      return calin::ix::iact_data::telescope_event::TRIGGER_INTERNAL_FLASHER;
    case 0x10:
      return calin::ix::iact_data::telescope_event::TRIGGER_UCTS_AUX;
    case 0x20:
      return calin::ix::iact_data::telescope_event::TRIGGER_PEDESTAL;
    case 0x40:
      return calin::ix::iact_data::telescope_event::TRIGGER_SOFTWARE;
    case 0x80:
      return calin::ix::iact_data::telescope_event::TRIGGER_UNKNOWN;
    default:
      return calin::ix::iact_data::telescope_event::TRIGGER_MULTIPLE;
    }
  }
  return calin::ix::iact_data::telescope_event::TRIGGER_UNKNOWN;
}

void calin::iact_data::acada_event_decoder::
decode_swat_data(calin::ix::iact_data::telescope_event::SWATData* calin_swat_data,
  const AnyArray& cta_array)
{
  struct SWATMessageData {
    uint64_t trigger_id;
    uint64_t bunch_id;
    uint32_t trigger_type;
    uint32_t trigger_time_s;
    uint32_t trigger_time_qns;
    uint32_t readout_requested;
    uint32_t data_available;
    uint32_t hardware_stereo_trigger_mask;
  } __attribute__((__packed__));

  const auto& cta_swat_data = cta_array.data();
  if(cta_swat_data.size() != sizeof(SWATMessageData))
    throw std::runtime_error("SWAT data array not expected size : " 
      + std::to_string(cta_swat_data.size())  + " != "
      + std::to_string(sizeof(SWATMessageData)));
  const auto* swat_data =
    reinterpret_cast<const SWATMessageData*>(&cta_swat_data.front());

  calin_swat_data->set_trigger_id(swat_data->trigger_id);
  calin_swat_data->set_bunch_id(swat_data->bunch_id);
  calin_swat_data->set_trigger_type_code(swat_data->trigger_type);
  calin_swat_data->set_trigger_time_s(swat_data->trigger_time_s);
  calin_swat_data->set_trigger_time_qns(swat_data->trigger_time_qns);
  calin_swat_data->set_readout_requested(swat_data->readout_requested);
  calin_swat_data->set_data_available(swat_data->data_available);
  calin_swat_data->set_hardware_stereo_trigger_mask(swat_data->hardware_stereo_trigger_mask);

  calin_swat_data->set_mono_trigger(swat_data->trigger_type & 0x0001);
  calin_swat_data->set_stereo_trigger(swat_data->trigger_type & 0x0002);
  calin_swat_data->set_calibration_trigger(swat_data->trigger_type & 0x0004);
  calin_swat_data->set_photo_electron_trigger(swat_data->trigger_type & 0x0008);
  calin_swat_data->set_software_trigger(swat_data->trigger_type & 0x0010);
  calin_swat_data->set_pedestal_trigger(swat_data->trigger_type & 0x0020);
  calin_swat_data->set_slow_control_trigger(swat_data->trigger_type & 0x0040);
  calin_swat_data->set_neighbour_1_trigger(swat_data->trigger_type & 0x0100);
  calin_swat_data->set_neighbour_2_trigger(swat_data->trigger_type & 0x0200);
  calin_swat_data->set_neighbour_3_trigger(swat_data->trigger_type & 0x0400);
  calin_swat_data->set_neighbour_4_trigger(swat_data->trigger_type & 0x0800);
  calin_swat_data->set_neighbour_5_trigger(swat_data->trigger_type & 0x1000);
  calin_swat_data->set_neighbour_6_trigger(swat_data->trigger_type & 0x2000);
  calin_swat_data->set_neighbour_7_trigger(swat_data->trigger_type & 0x4000);
                                                
  calin_swat_data->set_trigger_type(calin::ix::iact_data::telescope_event::SwatTriggerType(swat_data->trigger_type));
}
