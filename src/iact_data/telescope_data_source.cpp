/*

   calin/iact_data/telescope_data_source.cpp -- Stephen Fegan -- 2016-01-08

   A supplier of single telescope data, for example:
   ix::iact_data::telescope_event::TelescopeEvent

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

//#define CALIN_TELESCOPE_DATA_SOURCE_NO_EXTERN
#include <iact_data/telescope_data_source.hpp>

using namespace calin::iact_data::telescope_data_source;

template class calin::io::data_source::DataSource<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;
template class calin::io::data_source::RandomAccessDataSource<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;
template class calin::io::data_source::ProtobufFileDataSource<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;
template class calin::io::data_source::BufferedDataSource<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;
template class calin::io::data_source::MultiThreadDataSourceBuffer<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;

template class calin::io::data_source::DataSink<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;
template class calin::io::data_source::ProtobufFileDataSink<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;

TelescopeRandomAccessDataSourceWithRunConfig::
~TelescopeRandomAccessDataSourceWithRunConfig()
{
  // nothing to see here
}

void calin::iact_data::telescope_data_source::
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

void calin::iact_data::telescope_data_source::
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
