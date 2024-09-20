/*

   calin/iact_data/nectarcam_acada_event_decoder_r1v0.cpp -- Stephen Fegan -- 2018-09-21

   A decoder of NectarCAM ACADA data in prototype R1 format

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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
#include <util/string.hpp>
#include <iact_data/acada_data_source.hpp>
#include <iact_data/unified_acada_event_decoder.hpp>
#include <provenance/system_info.hpp>
#include <iact_data/nectarcam_layout.hpp>
#include <iact_data/lstcam_layout.hpp>

using namespace calin::iact_data::acada_data_source;
using namespace calin::iact_data::unified_acada_event_decoder;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::util::log;

/*

        RRRRRRRRRRRRRRRRR     1111111                        1111111   
        R::::::::::::::::R   1::::::1                       1::::::1   
        R::::::RRRRRR:::::R 1:::::::1                      1:::::::1   
        RR:::::R     R:::::R111:::::1                      111:::::1   
          R::::R     R:::::R   1::::1vvvvvvv           vvvvvvv1::::1   
          R::::R     R:::::R   1::::1 v:::::v         v:::::v 1::::1   
          R::::RRRRRR:::::R    1::::1  v:::::v       v:::::v  1::::1   
          R:::::::::::::RR     1::::l   v:::::v     v:::::v   1::::l   
          R::::RRRRRR:::::R    1::::l    v:::::v   v:::::v    1::::l   
          R::::R     R:::::R   1::::l     v:::::v v:::::v     1::::l   
          R::::R     R:::::R   1::::l      v:::::v:::::v      1::::l   
          R::::R     R:::::R   1::::l       v:::::::::v       1::::l   
        RR:::::R     R:::::R111::::::111     v:::::::v     111::::::111
        R::::::R     R:::::R1::::::::::1      v:::::v      1::::::::::1
        R::::::R     R:::::R1::::::::::1       v:::v       1::::::::::1
        RRRRRRRR     RRRRRRR111111111111        vvv        111111111111
                                                               

*/

Unified_ACADACameraEventDecoder_R1v1::
Unified_ACADACameraEventDecoder_R1v1(const std::string& filename, const config_type& config):
  calin::iact_data::acada_event_decoder::ACADACameraEventDecoder<
    calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1>(),
  config_(config), filename_(filename)
{
  // nothing to see here
}

Unified_ACADACameraEventDecoder_R1v1::
~Unified_ACADACameraEventDecoder_R1v1()
{
  // nothing to see here
}

bool Unified_ACADACameraEventDecoder_R1v1::
decode(calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
  const calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1& cta_messages)
{
  // The decoder assumes that calin_events are reused - prefer to overwrite entries than clear
  // all of them which is potentially more efficient. Note we only overwrite entries that are
  // used here.

  // calin_event->Clear(); 

  const event_type* cta_event = cta_messages.event;

  if(cta_event == nullptr) {
    throw std::runtime_error("Unified_ACADACameraEventDecoder_R1v1::decode: "
      "No event found");
  }

    // optional fixed32  event_time_s                     =  5; // HighResTimestamp in DM
    // optional fixed32  event_time_qns                   =  6; // HighResTimestamp in DM
    // optional fixed32  num_channels                     =  7; // uint8 in DM
    // optional fixed32  num_samples                      =  8; // uint16 in DM
    // optional fixed32  num_pixels                       =  9; // uint16 in DM
    // optional fixed32  num_modules                      = 10; // missing in DM ? 
    // optional AnyArray waveform                         = 11; // uint16[num_channels*num_pixels*num_samples] in DM
    // optional AnyArray pixel_status                     = 12; // uint8[num_pixels] in DM
    // optional AnyArray first_cell_id                    = 13; // uint16[N] with N=1 for SST, 2120 for LST, not present otherwise
    // optional AnyArray module_hires_local_clock_counter = 14; // uint64[num_modules] in DM. 
    // optional AnyArray pedestal_intensity               = 15; // float32[num_pixels] in DM
    // optional fixed64  calibration_monitoring_id        = 16; // uint64 in DM
    // optional R1v1_debug.DebugEvent debug = 17;
 
  calin_event->set_telescope_id(cta_event->tel_id());
  calin_event->set_local_event_number(cta_event->event_id());
  switch(cta_event->event_type()) {
  case 0: calin_event->set_trigger_type(TRIGGER_EXTERNAL_FLASHER); break;
  case 1: calin_event->set_trigger_type(TRIGGER_INTERNAL_FLASHER); break;
  case 2: calin_event->set_trigger_type(TRIGGER_PEDESTAL); break;
  case 16: calin_event->set_trigger_type(TRIGGER_MUON); break;
  case 17: calin_event->set_trigger_type(TRIGGER_FORCED_BY_ARRAY); break;
  case 24: calin_event->set_trigger_type(TRIGGER_SOFTWARE); break;
  case 32:
  case 33: calin_event->set_trigger_type(TRIGGER_PHYSICS); break;
  default: calin_event->set_trigger_type(TRIGGER_UNKNOWN); break;
  }
  calin_event->set_array_trigger_received(cta_event->event_type() == 17);
  calin_event->set_array_event_number(cta_event->event_id());

  // **************************************************************************
  // Waveforms and pixel status
  // **************************************************************************

  if(cta_event->has_waveform() and cta_event->waveform().has_data() and 
    cta_event->has_pixel_status() and cta_event->pixel_status().has_data())
  {
    unsigned ngain = cta_event->num_channels();
    unsigned nsample = cta_event->num_samples();
    if(cta_event->num_pixels() != nchan_configured_) {
      throw std::runtime_error("Unified_ACADACameraEventDecoder_R1v1::decode: "
        "Number of pixels in event (" + std::to_string(cta_event->num_pixels()) + 
        ") differs from header (" + std::to_string(nchan_configured_) + ")");
    }
    if(cta_event->pixel_status().data().size() != nchan_configured_*sizeof(uint8_t)) {
      throw std::runtime_error("Unified_ACADACameraEventDecoder_R1v1::decode: "
        "Unexpected pixel-status array size (" + std::to_string(cta_event->pixel_status().data().size()) + 
        " != " + std::to_string(nchan_configured_*sizeof(uint8_t)) + ")");
    }
    if(cta_event->waveform().data().size() != ngain*nsample*nchan_configured_*sizeof(uint16_t)) {
      throw std::runtime_error("Unified_ACADACameraEventDecoder_R1v1::decode: "
        "Unexpected waveform array size (" + std::to_string(cta_event->waveform().data().size()) + 
        " != " + std::to_string(ngain*nsample*nchan_configured_*sizeof(uint16_t)) + ")");
    }
    const uint8_t* pix_status =
      reinterpret_cast<const uint8_t*>(cta_event->pixel_status().data().data());
    const int16_t* waveforms =
      reinterpret_cast<const int16_t*>(cta_event->waveform().data().data());

    if(ngain == 2) {
      copy_single_gain_waveforms(calin_event, waveforms, pix_status, nsample,
        calin_event->mutable_high_gain_image()->mutable_camera_waveforms(),
        0x04, "high");
      copy_single_gain_waveforms(calin_event, waveforms+nchan_configured_*nsample, pix_status, nsample,
        calin_event->mutable_low_gain_image()->mutable_camera_waveforms(),
        0x08, "low");
      calin_event->clear_image();
    } else if(ngain == 1) {
      calin_event->clear_high_gain_image();
      calin_event->clear_low_gain_image();
      copy_single_gain_waveforms(calin_event, waveforms, pix_status, nsample,
        calin_event->mutable_image()->mutable_camera_waveforms(),
        0x0C, "mixed");
    } else {
      throw std::runtime_error("Unified_ACADACameraEventDecoder_R1v1::decode: "
        "Unsupported number of gains : " + std::to_string(ngain));
    }
  } else {
    calin_event->clear_high_gain_image();
    calin_event->clear_low_gain_image();
    calin_event->clear_image();
  }

  // **************************************************************************
  // Camera time
  // **************************************************************************

  calin_event->mutable_absolute_event_time()->set_time_ns(
    uint64_t(cta_event->event_time_s())*uint64_t(1000000000)
      + uint64_t(cta_event->event_time_qns()>>2));
  calin_event->mutable_elapsed_event_time()->set_time_ns(
    calin_event->absolute_event_time().time_ns() - run_start_time_);

  calin_event->mutable_camera_clock_index()->Resize(ncamera_clock_,-1);
  calin_event->clear_camera_clock();

  calin_event->set_camera_clock_index(0,calin_event->camera_clock_size());
  auto* calin_clock = calin_event->add_camera_clock();
  calin_clock->set_clock_id(0);
  calin_clock->set_time_value(calin_event->absolute_event_time().time_ns());
  calin_clock->set_time_sequence_id(0);
  calin_clock->set_time_value_may_be_suspect(false);

  // **************************************************************************
  // CDTS data
  // **************************************************************************

  if(cta_event->has_debug() and cta_event->debug().has_cdts_data() 
      and cta_event->debug().cdts_data().has_data()) {
    calin::iact_data::acada_event_decoder::decode_cdts_data(
      calin_event->mutable_cdts_data(), cta_event->debug().cdts_data());

    const auto& cdts = calin_event->cdts_data();

    if(cdts.event_counter() != cta_event->event_id()) {
      calin_event->clear_cdts_data();
    }
  } else {
    calin_event->clear_cdts_data();
  }

  if(calin_event->has_cdts_data())
  {
    const auto& cdts = calin_event->cdts_data();

    bool clock_may_be_suspect =
      (calin_event->cdts_data().white_rabbit_status() & 0x01) == 0;

    calin_event->set_camera_clock_index(1,calin_event->camera_clock_size());
    auto* calin_clock = calin_event->add_camera_clock();
    calin_clock->set_clock_id(1);
    calin_clock->set_time_value(cdts.ucts_timestamp());
    calin_clock->set_time_sequence_id(0);
    calin_clock->set_time_value_may_be_suspect(clock_may_be_suspect);

    calin_event->set_camera_clock_index(2,calin_event->camera_clock_size());
    calin_clock = calin_event->add_camera_clock();
    calin_clock->set_clock_id(2);
    calin_clock->set_time_value(cdts.clock_counter());
    calin_clock->set_time_sequence_id(cdts.pps_counter());
    calin_clock->set_time_value_may_be_suspect(clock_may_be_suspect);

    calin_event->set_camera_clock_index(3,calin_event->camera_clock_size());
    calin_clock = calin_event->add_camera_clock();
    calin_clock->set_clock_id(3);
    calin_clock->set_time_value(cdts.pps_counter());
    calin_clock->set_time_sequence_id(0);
    calin_clock->set_time_value_may_be_suspect(clock_may_be_suspect);

    calin_event->set_camera_clock_index(4,calin_event->camera_clock_size());
    calin_clock = calin_event->add_camera_clock();
    calin_clock->set_clock_id(4);
    calin_clock->set_time_value(cdts.pps_counter()*10000000ULL + cdts.clock_counter());
    calin_clock->set_time_sequence_id(0);
    calin_clock->set_time_value_may_be_suspect(clock_may_be_suspect);
  } else {
    calin_event->set_camera_clock_index(1,-1);
    calin_event->set_camera_clock_index(2,-1);
    calin_event->set_camera_clock_index(3,-1);
    calin_event->set_camera_clock_index(4,-1);
  }

  // **************************************************************************
  // TIB data
  // **************************************************************************

  // SJF : TIB has bit 0x01 in effect, contrary to what is in CamerasToACTL - 2020-06-28
  if(cta_event->has_debug() and cta_event->debug().has_tib_data() 
      and cta_event->debug().tib_data().has_data())
  {
    calin::iact_data::acada_event_decoder::decode_tib_data(
      calin_event->mutable_tib_data(), cta_event->debug().tib_data());

    const auto& tib = calin_event->tib_data();

    if(tib.event_counter() != cta_event->event_id()) {
      calin_event->clear_tib_data();
    }
  } else {
    calin_event->clear_tib_data();
  }

  if(calin_event->has_tib_data()) {
    const auto& tib = calin_event->tib_data();

    calin_event->set_camera_clock_index(5,calin_event->camera_clock_size());
    auto* calin_clock = calin_event->add_camera_clock();
    calin_clock->set_clock_id(5);
    calin_clock->set_time_value(tib.clock_counter());
    calin_clock->set_time_sequence_id(tib.pps_counter());

    calin_event->set_camera_clock_index(6,calin_event->camera_clock_size());
    calin_clock = calin_event->add_camera_clock();
    calin_clock->set_clock_id(6);
    calin_clock->set_time_value(tib.pps_counter());
    calin_clock->set_time_sequence_id(0);

    calin_event->set_camera_clock_index(7,calin_event->camera_clock_size());
    calin_clock = calin_event->add_camera_clock();
    calin_clock->set_clock_id(7);
    calin_clock->set_time_value(tib.pps_counter()*10000000ULL + tib.clock_counter());
    calin_clock->set_time_sequence_id(0);
  } else {
    calin_event->set_camera_clock_index(5,-1);
    calin_event->set_camera_clock_index(6,-1);
    calin_event->set_camera_clock_index(7,-1);
  }

  return true;
}

bool Unified_ACADACameraEventDecoder_R1v1::
decode_run_config(
  calin::ix::iact_data::telescope_run_configuration:: TelescopeRunConfiguration* calin_run_config,
  const calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1& cta_messages)
{
  calin_run_config->Clear();

  const data_stream_type* cta_data_stream = cta_messages.data_stream;
  const header_type* cta_run_header = cta_messages.header;
  // const event_type* cta_event = cta_messages.event;

  if(cta_run_header == nullptr) {
    throw std::runtime_error("Unified_ACADACameraEventDecoder_R1v1::decode_run_config: "
      "No camera configuration found");
  }
  if(cta_data_stream == nullptr) {
    throw std::runtime_error("Unified_ACADACameraEventDecoder_R1v1::decode_run_config: "
      "No telescope data stream header found");
  }

  calin_run_config->set_data_transcoder(
    "calin::iact_data::unified_event_decoder::Unified_ACADACameraEventDecoder_R1v1");
  calin_run_config->set_filename(filename_);
  calin_run_config->add_fragment_filename(filename_);
  calin_run_config->set_run_number(cta_run_header->local_run_id());
  run_start_time_ = uint64_t(cta_run_header->config_time_s() * 1e9);
  calin_run_config->mutable_run_start_time()->set_time_ns(run_start_time_);
  calin_run_config->set_telescope_id(cta_data_stream->tel_id());  
  calin_run_config->set_scheduling_block_id(cta_data_stream->sb_id());
  calin_run_config->set_observation_id(cta_data_stream->obs_id());

  // **************************************************************************
  // Camera layout
  // **************************************************************************

  switch(config_.camera_type()) {
  case calin::ix::iact_data::cta_data_source::AUTO_DETECT:
    // Hard-code for NectarCAM for the moment - eventually use telescope id
    nectarcam_layout::nectarcam_layout(calin_run_config->mutable_camera_layout());
    break;
  case calin::ix::iact_data::cta_data_source::NECTARCAM:
    nectarcam_layout::nectarcam_layout(calin_run_config->mutable_camera_layout());
    break;
  case calin::ix::iact_data::cta_data_source::LSTCAM:
    lstcam_layout::lstcam_layout(calin_run_config->mutable_camera_layout());
    break;
  default:
    throw std::runtime_error("Unified_ACADACameraEventDecoder_R1v1::decode_run_config: "
      "Unsupported camera type :" + calin::ix::iact_data::cta_data_source::CameraType_Name(config_.camera_type()));
  }

  calin_run_config->mutable_camera_layout()->clear_camera_clock_name();
  calin_run_config->mutable_camera_layout()->clear_camera_clock_frequency();
  calin_run_config->mutable_camera_layout()->clear_module_clock_name();
  calin_run_config->mutable_camera_layout()->clear_module_clock_frequency();

  #define ADD_CAMERA_CLOCK(name, freq) \
    calin_run_config->mutable_camera_layout()->add_camera_clock_name(name); \
    calin_run_config->mutable_camera_layout()->add_camera_clock_frequency(freq)

  #define ADD_MODULE_CLOCK(name, freq) \
    calin_run_config->mutable_camera_layout()->add_module_clock_name(name); \
    calin_run_config->mutable_camera_layout()->add_module_clock_frequency(freq)

  ADD_CAMERA_CLOCK("EVB timestamp",                         1.0e9);

  ADD_CAMERA_CLOCK("UCTS timestamp",                        1.0e9);
  ADD_CAMERA_CLOCK("UCTS 10MHz counter",                    1.0e7);
  ADD_CAMERA_CLOCK("UCTS pps counter",                      1.0);
  ADD_CAMERA_CLOCK("UCTS combined 10MHz and pps counter",   1.0e7);

  ADD_CAMERA_CLOCK("TIB 10MHz counter",                     1.0e7);
  ADD_CAMERA_CLOCK("TIB pps counter",                       1.0);
  ADD_CAMERA_CLOCK("TIB combined 10MHz and pps counter",    1.0e7);

  ncamera_clock_ = 8;
  nmodule_clock_ = 0;

  // **************************************************************************
  // Build configured-module list
  // **************************************************************************

  nmod_configured_ = cta_run_header->num_modules();
  if(nmod_configured_ > 0) {
    if(not cta_run_header->has_module_id_map() or 
        not cta_run_header->module_id_map().has_data() or
        cta_run_header->module_id_map().data().size() != nmod_configured_*sizeof(uint16_t)) {
      throw std::runtime_error("bool Unified_ACADACameraEventDecoder_R1v1::decode_run_config: "
        "configured-module map inconsistent");
    } else {
      calin_run_config->mutable_configured_module_index()->Resize(
        calin_run_config->camera_layout().module_size(), -1);
      const uint16_t* mod_id = reinterpret_cast<const uint16_t*>(
        &cta_run_header->module_id_map().data().front());
      for(unsigned imod=0;imod<nmod_configured_;imod++) {
        if(mod_id[imod] >= calin_run_config->configured_module_index_size()) {
          throw std::runtime_error("bool Unified_ACADACameraEventDecoder_R1v1::decode_run_config: "
            "module id out of range :" + std::to_string(mod_id[imod]));
        }
        calin_run_config->add_configured_module_id(mod_id[imod]);
        calin_run_config->set_configured_module_index(mod_id[imod],imod);
      }
    }
  }

  // **************************************************************************
  // Build configured-channel list
  // **************************************************************************

  nchan_configured_ = cta_run_header->num_pixels();
  if(nchan_configured_ > 0) {
    if(not cta_run_header->has_pixel_id_map() or
        not cta_run_header->pixel_id_map().has_data() or
        cta_run_header->pixel_id_map().data().size() != nchan_configured_*sizeof(uint16_t)) {
      throw std::runtime_error("bool Unified_ACADACameraEventDecoder_R1v1::decode_run_config: "
        "configured-pixel map not consistent");
    } else {
      calin_run_config->mutable_configured_channel_index()->Resize(
        calin_run_config->camera_layout().channel_size(), -1);
      const uint16_t* chan_id = reinterpret_cast<const uint16_t*>(
        &cta_run_header->pixel_id_map().data().front());
      for(unsigned ichan=0;ichan<nchan_configured_;ichan++) {
        if(chan_id[ichan] >= calin_run_config->configured_channel_index_size()) {
          throw std::runtime_error("bool Unified_ACADACameraEventDecoder_R1v1::decode_run_config: "
            "channel id out of range :" + std::to_string(chan_id[ichan]));
        }
        calin_run_config->add_configured_channel_id(chan_id[ichan]);
        calin_run_config->set_configured_channel_index(chan_id[ichan],ichan);
      }
    }
  }

  // **************************************************************************
  // Waveform details
  // **************************************************************************

  calin_run_config->set_num_samples(cta_run_header->num_samples_nominal());
  calin_run_config->set_num_samples_long(cta_run_header->num_samples_long());
  calin_run_config->set_nominal_sampling_frequency(1000.0);
  calin_run_config->set_waveform_scale(cta_data_stream->waveform_scale());
  calin_run_config->set_waveform_offset(cta_data_stream->waveform_offset());

  // **************************************************************************
  // Various configuration parameters
  // **************************************************************************

  calin_run_config->set_configuration_id(cta_run_header->camera_config_id());
  calin_run_config->set_data_model_version(cta_run_header->data_model_version());
  calin_run_config->set_calibration_service_id(cta_run_header->calibration_service_id());
  calin_run_config->set_calibration_algorithm_id(cta_run_header->calibration_algorithm_id());
  (*calin_run_config->mutable_configuration_elements())["data_model_version"] = 
    cta_run_header->data_model_version();
  if(cta_run_header->has_debug()) {
    auto& cta_run_header_debug = cta_run_header->debug();
    if(cta_run_header_debug.cs_serial().size())
      (*calin_run_config->mutable_configuration_elements())["camera_server_serial_number"] = 
        calin::util::string::chomp(cta_run_header_debug.cs_serial());
    if(cta_run_header_debug.evb_version().size())
      (*calin_run_config->mutable_configuration_elements())["event_builder_version"] = 
        cta_run_header_debug.evb_version();
    if(cta_run_header_debug.cdhs_version().size())
      (*calin_run_config->mutable_configuration_elements())["common_data_handling_system_version"] = 
        cta_run_header_debug.cdhs_version();
  }

  // **************************************************************************
  // SERIALIZE RAW DATA
  // **************************************************************************

  if(config_.include_serialized_raw_data())
  {
    calin_run_config->set_serialized_raw_header_type(
      SerializedRawHeaderType::SERIALIZED_RAW_HEADER_ACADA_PROTOBUF_R1V1);
    if(cta_run_header)
      cta_run_header->SerializeToString(calin_run_config->mutable_serialized_raw_header());
    if(cta_data_stream)
      cta_data_stream->SerializeToString(calin_run_config->mutable_serialized_data_stream());
  } else {
    calin_run_config->set_serialized_raw_header_type(
      SerializedRawHeaderType::SERIALIZED_RAW_HEADER_NONE);
  }

  return true;
}

Unified_ACADACameraEventDecoder_R1v1* Unified_ACADACameraEventDecoder_R1v1::clone() const
{
  return new Unified_ACADACameraEventDecoder_R1v1(*this);
}

void Unified_ACADACameraEventDecoder_R1v1::
copy_single_gain_waveforms(
  const calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event, 
  const int16_t* cta_waveforms, const uint8_t* cta_pixel_mask, unsigned nsample,
  calin::ix::iact_data::telescope_event::Waveforms* calin_waveforms,
  uint8_t has_gain_mask, const std::string& which_gain) const
{
  bool all_channels_present = true;

  calin_waveforms->mutable_channel_index()->Reserve(nchan_configured_);
  calin_waveforms->mutable_channel_id()->Reserve(nchan_configured_);

  std::string* calin_wf_raw_data_string = calin_waveforms->mutable_raw_samples_array();
  unsigned simd_vec_size = calin::provenance::system_info::the_host_info()->simd_vec_size();
  calin_wf_raw_data_string->resize(nsample*nchan_configured_*sizeof(int16_t) + simd_vec_size);
  char* cp = &calin_wf_raw_data_string->front();
  std::fill(cp+nsample*nchan_configured_*sizeof(int16_t), cp+calin_wf_raw_data_string->size(), int8_t(0));
  int16_t* calin_wf_raw_data = reinterpret_cast<int16_t*>(cp);

  for(unsigned ichan=0;ichan<nchan_configured_;ichan++)
  {
    if(cta_pixel_mask[ichan] & has_gain_mask) {
      std::copy(cta_waveforms, cta_waveforms+nsample, calin_wf_raw_data);
      calin_wf_raw_data += nsample;
      calin_waveforms->add_channel_index(calin_waveforms->channel_id_size());
      if((has_gain_mask == 0x04) or (has_gain_mask == 0x0C and cta_pixel_mask[ichan] == 0x04)) {
        calin_waveforms->add_channel_signal_type(calin::ix::iact_data::telescope_event::SIGNAL_HIGH_GAIN);
      } else if ((has_gain_mask == 0x08) or (has_gain_mask == 0x0C and cta_pixel_mask[ichan] == 0x08)) {
        calin_waveforms->add_channel_signal_type(calin::ix::iact_data::telescope_event::SIGNAL_LOW_GAIN);
      } else {
        throw std::runtime_error("Unified_ACADACameraEventDecoder_R1v1::copy_single_gain_waveforms: Unhandled pixel mask: " +
          std::to_string(unsigned(has_gain_mask)) + " / " + std::to_string(unsigned(cta_pixel_mask[ichan])));
      }
      calin_waveforms->add_channel_id(ichan);
      if(config_.separate_channel_waveforms()) {
        auto* calin_samp = calin_waveforms->add_waveform()->mutable_samples();
        calin_samp->Reserve(nsample);
        for(unsigned isample=0;isample<nsample;isample++)
          calin_samp->Add(*cta_waveforms++);
      } else {
        cta_waveforms += nsample;
      }
    } else {
      std::fill(calin_wf_raw_data, calin_wf_raw_data+nsample, 0);
      calin_wf_raw_data += nsample;
      all_channels_present = false;
      calin_waveforms->add_channel_index(-1);
      calin_waveforms->add_channel_signal_type(
        calin::ix::iact_data::telescope_event::SIGNAL_NONE);
      cta_waveforms += nsample;
    }
  }

  calin_waveforms->set_all_channels_present(all_channels_present);
}
