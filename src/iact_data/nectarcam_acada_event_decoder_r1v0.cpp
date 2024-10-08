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
#include <iact_data/acada_data_source.hpp>
#include <iact_data/nectarcam_acada_event_decoder.hpp>
#include <iact_data/nectarcam_layout.hpp>
#include <iact_data/nectarcam_configuration.hpp>
#include <provenance/system_info.hpp>

using namespace calin::iact_data::acada_data_source;
using namespace calin::iact_data::nectarcam_acada_event_decoder;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::util::log;

#define TEST_ANYARRAY_TYPES 0

/*

    RRRRRRRRRRRRRRRRR     1111111                              000000000     
    R::::::::::::::::R   1::::::1                            00:::::::::00   
    R::::::RRRRRR:::::R 1:::::::1                          00:::::::::::::00 
    RR:::::R     R:::::R111:::::1                         0:::::::000:::::::0
      R::::R     R:::::R   1::::1vvvvvvv           vvvvvvv0::::::0   0::::::0
      R::::R     R:::::R   1::::1 v:::::v         v:::::v 0:::::0     0:::::0
      R::::RRRRRR:::::R    1::::1  v:::::v       v:::::v  0:::::0     0:::::0
      R:::::::::::::RR     1::::l   v:::::v     v:::::v   0:::::0 000 0:::::0
      R::::RRRRRR:::::R    1::::l    v:::::v   v:::::v    0:::::0 000 0:::::0
      R::::R     R:::::R   1::::l     v:::::v v:::::v     0:::::0     0:::::0
      R::::R     R:::::R   1::::l      v:::::v:::::v      0:::::0     0:::::0
      R::::R     R:::::R   1::::l       v:::::::::v       0::::::0   0::::::0
    RR:::::R     R:::::R111::::::111     v:::::::v        0:::::::000:::::::0
    R::::::R     R:::::R1::::::::::1      v:::::v          00:::::::::::::00 
    R::::::R     R:::::R1::::::::::1       v:::v             00:::::::::00   
    RRRRRRRR     RRRRRRR111111111111        vvv                000000000     

*/

NectarCam_ACADACameraEventDecoder_R1v0::NectarCam_ACADACameraEventDecoder_R1v0(
    const std::string& filename, unsigned run_number, const config_type& config):
  calin::iact_data::acada_event_decoder::ACADACameraEventDecoder<
    calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0>(), 
  config_(config), filename_(filename), run_number_(run_number)
{
  if(config_.demand_configured_module_id_size() != 0)
    LOG(WARNING) << "Decoder option \"demand_configured_module_id_size\" not "
      "supported in R1 data at this time.";

  if(config_.exchange_gain_channels() ==
      ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig::EXCHANGE_GAIN_MODE_FORCED)
    LOG(WARNING) << "Decoder option \"exchange_gain_channels\" not supported "
      "in R1 data at this time.";
}

NectarCam_ACADACameraEventDecoder_R1v0::~NectarCam_ACADACameraEventDecoder_R1v0()
{
  // nothing to see here
}

bool NectarCam_ACADACameraEventDecoder_R1v0::decode(
  calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
  const calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0& cta_messages)
{
  const event_type* cta_event = cta_messages.event;

  if(!cta_event->has_nectarcam())
    throw(std::runtime_error("NectarCam_ACADACameraEventDecoder_R1v0::decode: "
      "ACADA event does not have NectarCAM extension"));

  calin_event->set_telescope_id(telescope_id_);
  calin_event->set_local_event_number(cta_event->tel_event_id());
  calin_event->set_trigger_type(TRIGGER_UNKNOWN);
  calin_event->set_array_trigger_received(false);
  calin_event->set_array_event_number(cta_event->event_id());
  //calin_event->local_clock_time
  calin_event->set_image_treatment(TREATMENT_SCIENCE);
  calin_event->set_configuration_id(cta_event->configuration_id());
  calin_event->set_pedestal_dataset_id(cta_event->ped_id());

  // ---------------------------------------------------------------------------
  // DONT TRUST TIME FROM CAMERA SERVER UNTIL WE KNOW WHAT IT IS
  // ---------------------------------------------------------------------------
  // calin_event->mutable_absolute_event_time()->set_time_ns(
  //   uint64_t(cta_event->trigger_time_s())*uint64_t(1000000000)
  //   + uint64_t(cta_event->trigger_time_qns()>>2));
  // calin_event->mutable_elapsed_event_time()->set_time_ns(
  //   calin_event->absolute_event_time().time_ns() - run_start_time_);
  calin_event->clear_absolute_event_time();
  calin_event->clear_elapsed_event_time();

  bool all_modules_present = true;
  if(cta_event->nectarcam().has_module_status())
  {
    const auto& cta_status = cta_event->nectarcam().module_status();
#if TEST_ANYARRAY_TYPES
    if(cta_status.type() != DataModel::AnyArray::U8)
      throw std::runtime_error("NectarCam_ACADACameraEventDecoder_R1v0::decode: "
        "Camera status type not U8");
#endif
    if(nmod_ != cta_status.data().size())
      throw std::runtime_error("NectarCam_ACADACameraEventDecoder_R1v0::decode: "
        "Module status array size does not match number of nodukes.");
    const auto* mod_status =
      reinterpret_cast<const uint8_t*>(&cta_status.data().front());
    for(unsigned imod=0, mod_index=0;imod<nmod_;imod++)
    {
      if(*(mod_status++)&0x01)
      {
        calin_event->add_module_index(mod_index);
        calin_event->add_module_id(imod);
        mod_index++;
      }
      else
      {
        calin_event->add_module_index(-1);
        all_modules_present = false;
      }
    }
  }
  else
  {
    throw(std::runtime_error("NectarCam_ACADACameraEventDecoder_R1v0::decode: "
      "ACADA event does not have NectarCAM module_status"));
  }
  calin_event->set_all_modules_present(all_modules_present);

  // ==========================================================================
  //
  // TRANSFER IMAGE DATA
  //
  // ==========================================================================

  if(cta_event->has_waveform() and cta_event->has_pixel_status())
  {
    const unsigned npix = nmod_*7;
    unsigned single_gain_dataset_size = npix*nsample_*sizeof(int16_t);
    if(cta_event->pixel_status().data().size() != npix)
      throw(std::runtime_error("NectarCam_ACADACameraEventDecoder_R1v0::decode: "
        "Pixel status array incorrect size: "
        + std::to_string(cta_event->pixel_status().data().size())
        + ", expected: " + std::to_string(npix)));

    const uint8_t* pix_status =
      reinterpret_cast<const uint8_t*>(cta_event->pixel_status().data().data());
    const int16_t* waveforms =
      reinterpret_cast<const int16_t*>(cta_event->waveform().data().data());

    if(cta_event->waveform().data().size() == 2*single_gain_dataset_size) {
      copy_single_gain_waveforms(calin_event, waveforms, pix_status,
        calin_event->mutable_high_gain_image()->mutable_camera_waveforms(),
        0x04, "high");
      copy_single_gain_waveforms(calin_event, waveforms+npix*nsample_, pix_status,
        calin_event->mutable_low_gain_image()->mutable_camera_waveforms(),
        0x08, "low");
    } else if(cta_event->waveform().data().size() == single_gain_dataset_size) {
      copy_single_gain_waveforms(calin_event, waveforms, pix_status,
        calin_event->mutable_image()->mutable_camera_waveforms(),
        0x0C, "mixed");
    } else {
      throw(std::runtime_error("NectarCam_ACADACameraEventDecoder_R1v0::decode: "
        "Samples array incorrect size: "
        + std::to_string(cta_event->waveform().data().size())
        + ", expected: " + std::to_string(2*single_gain_dataset_size)
        + " (dual gain mode) or " + std::to_string(single_gain_dataset_size)
        + " (mixed gain)"));
    }
  }

  // ==========================================================================
  //
  // DECODE NECTARCAM COUNTERS
  //
  // ==========================================================================

  uint64_t mod_clock_sum = 0;
  uint64_t mod_clock_seq_sum = 0;
  int64_t mod_clock_num = 0;
  bool mod_clock_is_suspect = false;

  if(cta_event->nectarcam().has_counters())
  {
    struct NectarCounters {
      /*  4 */ uint32_t global_event_counter;
      /*  6 */ uint16_t bunch_counter;
      /*  8 */ uint16_t event_counter;
      /* 12 */ uint32_t ts1;
      /* 13 */ int8_t   ts2_event;
      /* 14 */ int8_t   ts2_bunch;
      /* 16 */ uint16_t ts2_empty;
    }__attribute__((packed));

    struct NectarCountersWithTriggerPattern {
      /*  4 */ uint32_t global_event_counter;
      /*  6 */ uint16_t bunch_counter;
      /*  8 */ uint16_t event_counter;
      /* 12 */ uint32_t ts1;
      /* 13 */ int8_t   ts2_event;
      /* 14 */ int8_t   ts2_bunch;
      /* 16 */ uint16_t ts2_empty;
      /* 20 */ uint32_t trigger_pattern;
    }__attribute__((packed));

    const auto& cta_counters = cta_event->nectarcam().counters();
#if TEST_ANYARRAY_TYPES
    if(cta_counters.type() != DataModel::AnyArray::U16)
      throw std::runtime_error("Camera counters type not U16");
#endif

    calin::ix::iact_data::telescope_event::ChannelTriggerMap* trigger_map = nullptr;
    unsigned data_size = sizeof(NectarCounters);
    unsigned data_version = 0;
    if(cta_counters.data().size() != nmod_ * data_size) {
      data_size = sizeof(NectarCountersWithTriggerPattern);
      data_version = 1;
      if(cta_counters.data().size() != nmod_ * data_size) {
        throw std::runtime_error("Camera counters data array correct size. "
          "Expected " + std::to_string(nmod_*sizeof(NectarCounters)) + " or " +
          std::to_string(nmod_*sizeof(NectarCountersWithTriggerPattern)) +
          ", got " + std::to_string(cta_counters.data().size()));
      }

      trigger_map = calin_event->mutable_trigger_map();
      trigger_map->mutable_hit_channel_id()->Reserve(20);
      trigger_map->mutable_trigger_image()->Resize(nmod_ * 7, 0x80000000);
    }

    const auto* mod_data = cta_counters.data().data();
    for(unsigned imod=0;imod<nmod_;imod++, mod_data+=data_size)
    {
      const auto* mod_counter = reinterpret_cast<const NectarCounters*>(mod_data);

      if(imod < static_cast<unsigned>(calin_event->module_index_size()) and
        calin_event->module_index(imod) == -1)continue;

      auto* module_counters = calin_event->add_module_counter();
      module_counters->set_module_id(imod);
#define add_mod_counter(id,val) \
      { \
        module_counters->add_counter_id(id); \
        module_counters->add_counter_value(val); \
      }
      add_mod_counter(0, mod_counter->global_event_counter);
      add_mod_counter(1, mod_counter->bunch_counter);
      add_mod_counter(2, mod_counter->event_counter);
      add_mod_counter(3, mod_counter->ts1);
      add_mod_counter(4, mod_counter->ts2_bunch);
      add_mod_counter(5, mod_counter->ts2_event);
      add_mod_counter(6, mod_counter->ts2_empty);

      auto* module_data = calin_event->add_module_data()->mutable_nectarcam();
      module_data->set_module_id(imod);
      module_data->set_global_event_counter(mod_counter->global_event_counter);
      module_data->set_bunch_counter(mod_counter->bunch_counter);
      module_data->set_event_counter(mod_counter->event_counter);
      module_data->set_ts1(mod_counter->ts1);
      module_data->set_ts2_bunch(mod_counter->ts2_bunch);
      module_data->set_ts2_event(mod_counter->ts2_event);
      module_data->set_ts2_empty(mod_counter->ts2_empty);
      module_data->set_version(data_version);

      if(data_version > 0) {
        const auto* mod_trigger = reinterpret_cast<const NectarCountersWithTriggerPattern*>(mod_data);
        module_data->set_trigger_pattern(mod_trigger->trigger_pattern & 0x7F7F7F7F);

        for(unsigned imodchan=0; imodchan<7; imodchan++) {
          // Shuffle trigger pattern for each channel into its own bit stream
          // As an example, if imodchan=4
          // and pattern starts as               XXX4XXXX XXX3XXXX XXX2XXXX XXX1XXXX
          // first apply shift and mask giving   00000004 00000003 00000002 00000001
          uint32_t chan_trigger = (mod_trigger->trigger_pattern>>imodchan) & 0x01010101;
          chan_trigger |= chan_trigger >> 7;  // 00000004 00000043 00000032 00000021
          chan_trigger |= chan_trigger >> 14; // 00000004 00000043 00000432 00004321
          chan_trigger &= 0x0F;               // 00000000 00000000 00000000 00004321
          trigger_map->set_trigger_image(imod*7+imodchan, chan_trigger);
          if(chan_trigger) {
            trigger_map->add_hit_channel_id(imod*7+imodchan);
          }
        }
      } else {
        module_data->set_trigger_pattern(0x80808080);
      }

#define ts2_decode(x) int32_t(x)

      int32_t ts2_bunch = ts2_decode(mod_counter->ts2_bunch);
      int32_t ts2_event = ts2_decode(mod_counter->ts2_event);
      int32_t ts = mod_counter->ts1*8 + ts2_event - ts2_bunch;

      module_data->set_bunch_event_time(ts);

      bool clock_is_suspect = false;
      int32_t time_seq_id = mod_counter->bunch_counter;

      if(mod_counter->event_counter == 1 and mod_counter->ts1 > 124987500)
      {
        // Here we attempt to handle events where the TS1 value seems to
        // still be at the a hight count value but the PPS event_counter
        // says it should be the first event in the new PPS bunch. We can :
        // - flag its value as potentially suspect
        // - and/or try to fix the mismatch by decreasing its sequence id

        clock_is_suspect = true;
        time_seq_id -= 1;
      }

      auto* module_clocks = calin_event->add_module_clock();
      module_clocks->set_module_id(imod);

      // Clock that combines TS1 and TS2
      auto* clock = module_clocks->add_clock();
      clock->set_clock_id(0);
      clock->set_time_value(ts);
      clock->set_time_sequence_id(time_seq_id);
      clock->set_time_value_may_be_suspect(clock_is_suspect);

      // Clock using TS1 only
      clock = module_clocks->add_clock();
      clock->set_clock_id(1);
      clock->set_time_value(mod_counter->ts1);
      clock->set_time_sequence_id(time_seq_id);
      clock->set_time_value_may_be_suspect(clock_is_suspect);

      // Clock using PPS counter only
      clock = module_clocks->add_clock();
      clock->set_clock_id(2);
      clock->set_time_value(time_seq_id);
      clock->set_time_sequence_id(0);
      clock->set_time_value_may_be_suspect(clock_is_suspect);

      mod_clock_sum += ts;
      mod_clock_seq_sum += time_seq_id;
      mod_clock_num += 1;
      mod_clock_is_suspect |= clock_is_suspect;
    }
  }

  // ==========================================================================
  //
  // DECODE NECTARCAM CDTS DATA MESSAGE
  //
  // ==========================================================================

  // SJF : UCTS has bit 0x02 in effect, contrary to what is in CamerasToACTL - 2020-06-28
  if(cta_event->nectarcam().has_cdts_data()
    and cta_event->nectarcam().cdts_data().has_data())
  {
    calin::iact_data::acada_event_decoder::decode_cdts_data(
      calin_event->mutable_cdts_data(), cta_event->nectarcam().cdts_data());

    const auto& cdts = calin_event->cdts_data();

    if(cdts.event_counter() != cta_event->tel_event_id()) {
      calin_event->clear_cdts_data();
    }
  }

  if(calin_event->has_cdts_data())
  {
    const auto& cdts = calin_event->cdts_data();

    bool clock_may_be_suspect =
      (calin_event->cdts_data().white_rabbit_status() & 0x01) == 0;

    calin_event->add_camera_clock_index(calin_event->camera_clock_size());
    auto* calin_clock = calin_event->add_camera_clock();
    calin_clock->set_clock_id(0);
    calin_clock->set_time_value(cdts.ucts_timestamp());
    calin_clock->set_time_sequence_id(0);
    calin_clock->set_time_value_may_be_suspect(clock_may_be_suspect);

    calin_event->add_camera_clock_index(calin_event->camera_clock_size());
    calin_clock = calin_event->add_camera_clock();
    calin_clock->set_clock_id(1);
    calin_clock->set_time_value(cdts.clock_counter());
    calin_clock->set_time_sequence_id(cdts.pps_counter());
    calin_clock->set_time_value_may_be_suspect(clock_may_be_suspect);

    calin_event->add_camera_clock_index(calin_event->camera_clock_size());
    calin_clock = calin_event->add_camera_clock();
    calin_clock->set_clock_id(2);
    calin_clock->set_time_value(cdts.pps_counter());
    calin_clock->set_time_sequence_id(0);
    calin_clock->set_time_value_may_be_suspect(clock_may_be_suspect);

    calin_event->add_camera_clock_index(calin_event->camera_clock_size());
    calin_clock = calin_event->add_camera_clock();
    calin_clock->set_clock_id(3);
    calin_clock->set_time_value(cdts.pps_counter()*10000000ULL + cdts.clock_counter());
    calin_clock->set_time_sequence_id(0);
    calin_clock->set_time_value_may_be_suspect(clock_may_be_suspect);

    calin_event->set_is_muon_candidate(calin_event->is_muon_candidate() 
      || cdts.muon_candidate());
  } else {
    calin_event->add_camera_clock_index(-1);
    calin_event->add_camera_clock_index(-1);
    calin_event->add_camera_clock_index(-1);
    calin_event->add_camera_clock_index(-1);
  }

  // ==========================================================================
  //
  // DECODE NECTARCAM TIB DATA MESSAGE
  //
  // ==========================================================================

  // SJF : TIB has bit 0x01 in effect, contrary to what is in CamerasToACTL - 2020-06-28
  if(cta_event->nectarcam().has_tib_data()
    and cta_event->nectarcam().tib_data().has_data())
  {
    calin::iact_data::acada_event_decoder::decode_tib_data(
      calin_event->mutable_tib_data(), cta_event->nectarcam().tib_data());

    const auto& tib = calin_event->tib_data();

    if(tib.event_counter() != cta_event->tel_event_id()) {
      calin_event->clear_tib_data();
    }
  }

  if(calin_event->has_tib_data()) {
    const auto& tib = calin_event->tib_data();

    calin_event->add_camera_clock_index(calin_event->camera_clock_size());
    auto* calin_clock = calin_event->add_camera_clock();
    calin_clock->set_clock_id(4);
    calin_clock->set_time_value(tib.clock_counter());
    calin_clock->set_time_sequence_id(tib.pps_counter());

    calin_event->add_camera_clock_index(calin_event->camera_clock_size());
    calin_clock = calin_event->add_camera_clock();
    calin_clock->set_clock_id(5);
    calin_clock->set_time_value(tib.pps_counter());
    calin_clock->set_time_sequence_id(0);

    calin_event->add_camera_clock_index(calin_event->camera_clock_size());
    calin_clock = calin_event->add_camera_clock();
    calin_clock->set_clock_id(6);
    calin_clock->set_time_value(tib.pps_counter()*10000000ULL + tib.clock_counter());
    calin_clock->set_time_sequence_id(0);
  } else {
    calin_event->add_camera_clock_index(-1);
    calin_event->add_camera_clock_index(-1);
    calin_event->add_camera_clock_index(-1);
  }

  // ==========================================================================
  //
  // ADD MODULE CLOCK SUM AFTER CDTS AND TIB, IF IT IS VALID
  //
  // ==========================================================================

  if(mod_clock_num==calin_event->module_index_size()) {
    calin_event->add_camera_clock_index(calin_event->camera_clock_size());
    auto* calin_clock = calin_event->add_camera_clock();
    calin_clock->set_clock_id(7);
    calin_clock->set_time_value(mod_clock_sum);
    calin_clock->set_time_sequence_id(mod_clock_seq_sum);
    calin_clock->set_time_value_may_be_suspect(mod_clock_is_suspect);

    calin_event->add_camera_clock_index(calin_event->camera_clock_size());
    calin_clock = calin_event->add_camera_clock();
    calin_clock->set_clock_id(8);
    calin_clock->set_time_value(mod_clock_seq_sum*1000000000ULL + mod_clock_sum);
    calin_clock->set_time_sequence_id(0);
    calin_clock->set_time_value_may_be_suspect(mod_clock_is_suspect);
  } else {
    calin_event->add_camera_clock_index(-1);
    calin_event->add_camera_clock_index(-1);
  }

  // ==========================================================================
  //
  // FIGURE OUT EVENT TIME
  //
  // ==========================================================================

  if(calin_event->has_cdts_data()) {
    calin_event->mutable_absolute_event_time()->set_time_ns(
      calin_event->cdts_data().ucts_timestamp());
  } else {
    // Now what cat? Now what?
  }

  if(calin_event->has_absolute_event_time() and run_start_time_!=0) {
    calin_event->mutable_elapsed_event_time()->set_time_ns(
      calin_event->absolute_event_time().time_ns() - run_start_time_);
  }

  // ==========================================================================
  //
  // FIGURE OUT EVENT TYPE
  //
  // ==========================================================================

  if(calin_event->has_tib_data() and calin_event->has_cdts_data()) {
    calin_event->set_trigger_type(
      calin::iact_data::acada_event_decoder::determine_trigger_type(
        &calin_event->tib_data(), &calin_event->cdts_data()));
  } else if(calin_event->has_tib_data()) {
    calin_event->set_trigger_type(
      calin::iact_data::acada_event_decoder::determine_trigger_type(
        &calin_event->tib_data(), nullptr));
  } else if(calin_event->has_cdts_data()) {
    calin_event->set_trigger_type(
      calin::iact_data::acada_event_decoder::determine_trigger_type(
        nullptr, &calin_event->cdts_data()));
  } else {
    // Now what cat? Now what?
  }

  // ==========================================================================
  //
  // SERIALIZE RAW DATA
  //
  // ==========================================================================

  if(config_.include_serialized_raw_data())
  {
    calin_event->set_serialized_raw_event_type(
      SerializedRawEventType::SERIALIZED_RAW_EVENT_ACADA_PROTOBUF_R1V0);
    cta_event->SerializeToString(calin_event->mutable_serialized_raw_event());
  } else {
    calin_event->set_serialized_raw_event_type(
      SerializedRawEventType::SERIALIZED_RAW_EVENT_NONE);
  }

  return true;
}

bool NectarCam_ACADACameraEventDecoder_R1v0::decode_run_config(
  calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* calin_run_config,
  const calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0& cta_messages)
{
  const header_type* cta_run_header = cta_messages.header;
  const event_type* cta_event = cta_messages.event;

  calin_run_config->set_data_transcoder(
    "calin::iact_data::nectarcam_acada_event_decoder::NectarCam_ACADACameraEventDecoder_R1v0");
  calin_run_config->set_filename(filename_);
  calin_run_config->add_fragment_filename(filename_);
  calin_run_config->set_run_number(run_number_);

  switch(config_.camera_type())
  {
  case calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig::AUTOMATIC:
  case calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig::NECTARCAM:
    nectarcam_layout::nectarcam_layout(
      calin_run_config->mutable_camera_layout());
    break;
  case calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig::NECTARCAM_TESTBENCH_19CHANNEL:
  default:
    nectarcam_layout::nectarcam_19module_layout(
      calin_run_config->mutable_camera_layout());
    break;
  }

  // ---------------------------------------------------------------------------
  //
  // Try to read the NectarCam module configuration XML file
  //
  // ---------------------------------------------------------------------------

  std::vector<std::string> nmc_file_tried;
  std::string nmc_file;

  if(not config_.demand_nmc_xml_file().empty()) {
    if(calin::util::file::is_readable(config_.demand_nmc_xml_file())) {
      nmc_file = config_.demand_nmc_xml_file();
    } else {
      nmc_file_tried.emplace_back(config_.demand_nmc_xml_file());
    }
  } else {
    std::string nmc_dirname = calin::util::file::dirname(filename_);
    if(nmc_dirname == ".") {
      nmc_dirname = "";
    } else {
      nmc_dirname += '/';
    }
    std::string nmc_basename = calin::util::file::basename(filename_);
    while(not nmc_basename.empty()) {
      std::string test_file = nmc_dirname + nmc_basename + config_.nmc_xml_suffix();
      if(calin::util::file::is_readable(test_file)) {
        nmc_file = test_file;
        break;
      } else {
        nmc_file_tried.emplace_back(test_file);
      }
      nmc_basename = calin::util::file::strip_extension(nmc_basename);
    }
  }

  if(not nmc_file.empty()) {
    calin::ix::iact_data::nectarcam_configuration::NectarCamCameraConfiguration* nccc =
      calin::iact_data::nectarcam_configuration::decode_nmc_xml_file(nmc_file);
    if(nccc) {
      calin_run_config->mutable_nectarcam()->CopyFrom(*nccc);
      delete nccc;
    } else {
      LOG(WARNING) << "Could not parse NectarCAM module configuration XML file "
        << nmc_file;
    }
  } else {
    auto logger = LOG(WARNING);
    logger << "Could not find NectarCAM module configuration XML file, tried:\n";
    for(auto try_fn : nmc_file_tried) {
      logger << "- " << try_fn << '\n';
    }
    logger << "Set the \"demand_nmc_xml_file\" decoder option if you wish to "
      "specify a different file.";
  }

  // ---------------------------------------------------------------------------
  //
  // Extract CTA run header data
  //
  // ---------------------------------------------------------------------------

  if(cta_run_header)
  {
    telescope_id_ = cta_run_header->telescope_id();
    calin_run_config->set_telescope_id(telescope_id_);
    if(cta_run_header->cs_serial().size())
      (*calin_run_config->mutable_configuration_elements())["camera_server_serial_number"] = 
        cta_run_header->cs_serial();
    calin_run_config->set_configuration_id(cta_run_header->configuration_id());
    calin_run_config->set_data_model_version(cta_run_header->data_model_version());
    if(cta_run_header->data_model_version().size())
      (*calin_run_config->mutable_configuration_elements())["data_model_version"] = 
        cta_run_header->data_model_version();
    run_start_time_ = uint64_t(cta_run_header->date()) * uint64_t(1000000000);
    calin_run_config->mutable_run_start_time()->set_time_ns(run_start_time_);
    if(cta_run_header->has_nectarcam()) {
      const auto* cta_nectarcam_header = &cta_run_header->nectarcam();
      auto* ncc = calin_run_config->mutable_nectarcam();
      if(ncc->daq_mode() == "") {
        switch(cta_nectarcam_header->acquisition_mode()) {
          // See MST-CAM-TN-0008-LPNHE
          case 0: ncc->set_daq_mode("DAQCHARGE_C"); break;
          case 1: ncc->set_daq_mode("DAQCHARGE_T"); break;
          case 2: ncc->set_daq_mode("DAQSAMPLE"); break;
          case 3: ncc->set_daq_mode("DAQSAMPLE_C"); break;
          case 4: ncc->set_daq_mode("DAQSAMPLE_D"); break;
          case 5: ncc->set_daq_mode("DAQCHARGESAMPLE"); break;
          case 6: ncc->set_daq_mode("DAQCHARGESAMPLE_D"); break;
          default:
            ncc->set_daq_mode("Unrecognised DAQ mode: "
              + std::to_string(cta_nectarcam_header->acquisition_mode()));
            break;
        }
      }
      ncc->set_daq_processing_algorithms(cta_nectarcam_header->algorithms());
      ncc->set_idaq_version(cta_nectarcam_header->idaq_version());
      ncc->set_cdhs_version(cta_nectarcam_header->cdhs_version());
    }
  }

  // ---------------------------------------------------------------------------
  //
  // Get the list of configured modules
  //
  // ---------------------------------------------------------------------------

  // First try to determine number of modules - require this to come from data
  nmod_ = 0;
  if(nmod_==0 and cta_run_header and cta_run_header->has_nectarcam()) {
    nmod_ = cta_run_header->nectarcam().num_modules();
  }
  if(nmod_==0 and cta_event and cta_event->has_nectarcam()
      and cta_event->nectarcam().has_module_status()) {
    nmod_ = cta_event->nectarcam().module_status().data().size();
  }
  if(nmod_ == 0) {
    throw std::runtime_error("NectarCam_ACADACameraEventDecoder_R1v0::decode_run_config: "
      "Could not determine number of modules");
  }

  // Next get list of IDs
  unsigned nmod_camera = calin_run_config->camera_layout().module_size();
  std::vector<unsigned> config_mod_id;
  config_mod_id.reserve(nmod_camera);
  if(config_.demand_configured_module_id_size() != 0)
  {
    if(unsigned(config_.demand_configured_module_id_size()) != nmod_)
      throw std::runtime_error("NectarCam_ACADACameraEventDecoder_R1v0::decode_run_config: "
        "Demand module list size must equal number of modules in data.");
    for(unsigned imod=0;imod<nmod_;imod++) {
      unsigned mod_id = config_.demand_configured_module_id(imod);
      if(mod_id >= nmod_camera)
        throw std::runtime_error("NectarCam_ACADACameraEventDecoder_R1v0::decode_run_config: "
          "Demand module id out of range: " + std::to_string(mod_id) + " >= " +
          std::to_string(nmod_camera));
      config_mod_id.push_back(mod_id);
    }
  }
  else if(cta_run_header
    and cta_run_header->has_nectarcam()
    and cta_run_header->nectarcam().has_expected_modules_id()
    and cta_run_header->nectarcam().expected_modules_id().data().size() == nmod_*sizeof(uint16_t))
  {
    const uint16_t* mod_id =
      reinterpret_cast<const uint16_t*>(&
        cta_run_header->nectarcam().expected_modules_id().data().front());
    for(unsigned imod=0;imod<nmod_;imod++)config_mod_id.push_back(mod_id[imod]);
  }
  else if(calin_run_config->has_nectarcam() and
    unsigned(calin_run_config->nectarcam().module_size()) == nmod_)
  {
    for(unsigned imod=0; imod<nmod_; imod++) {
      unsigned mod_id = calin_run_config->nectarcam().module(imod).module_id();
      if(mod_id >= nmod_camera)
        throw std::runtime_error("NectarCam_ACADACameraEventDecoder_R1v0::decode_run_config: "
          "NMC module id out of range: " + std::to_string(mod_id) + " >= " +
          std::to_string(nmod_camera));
      config_mod_id.push_back(mod_id);
    }
  }
  else
  {
    for(unsigned imod=0;imod<nmod_;imod++)config_mod_id.push_back(imod);
  }

  for(unsigned mod_id=0; mod_id<nmod_camera; mod_id++)
  {
    calin_run_config->add_configured_module_index(-1);
    for(unsigned ipix=0; ipix<7; ipix++)
      calin_run_config->add_configured_channel_index(-1);
  }
  for(unsigned imod=0; imod<config_mod_id.size();imod++)
  {
    unsigned mod_id = config_mod_id[imod];
    calin_run_config->add_configured_module_id(mod_id);
    calin_run_config->set_configured_module_index(mod_id, imod);
    for(unsigned ipix=0; ipix<7; ipix++) {
      calin_run_config->add_configured_channel_id(mod_id*7+ipix);
      calin_run_config->set_configured_channel_index(mod_id*7+ipix, imod*7+ipix);
    }
  }

  calin_run_config->mutable_camera_layout()->set_camera_clock_frequency(7, double(config_mod_id.size())*1e9);
  calin_run_config->mutable_camera_layout()->set_camera_clock_frequency(8, double(config_mod_id.size())*1e9);

  calin::iact_data::telescope_data_source::report_run_configuration_problems(calin_run_config);

  // ---------------------------------------------------------------------------
  //
  // Determine nsample
  //
  // ---------------------------------------------------------------------------

  nsample_ = config_.demand_nsample();
  if(nsample_ == 0 and cta_run_header)
    nsample_ = cta_run_header->num_samples();
  if(nsample_ == 0 and calin_run_config->has_nectarcam() and
      calin_run_config->nectarcam().module_size()>0) {
    nsample_ = calin_run_config->nectarcam().module(0).num_samples();
    for(int imod=1; imod<calin_run_config->nectarcam().module_size(); imod++)
      if(calin_run_config->nectarcam().module(imod).num_samples() != nsample_)
        nsample_ = 0;
  }
#if 0
  if(nsample == 0 and cta_event and cta_event->has_logain() and
      cta_event->logain().has_waveforms())
    nsample = cta_event->logain().waveforms().num_samples();
  if(nsample == 0 and cta_event and cta_event->has_higain() and
      cta_event->higain().has_waveforms())
    nsample = cta_event->higain().waveforms().num_samples();
#endif
  if(nsample_ == 0) {
    throw std::runtime_error("NectarCam_ACADACameraEventDecoder_R1v0::decode_run_config: "
      "Could not determine number of samples");
  }
  calin_run_config->set_num_samples(nsample_);

  // ==========================================================================
  //
  // RUN SAMPLING FREQUENCY
  //
  // ==========================================================================

  double nominal_sampling_frequency = config_.demand_sampling_frequency();
  if(nominal_sampling_frequency == 0.0 and calin_run_config->has_nectarcam() and
      calin_run_config->nectarcam().module_size()>0) {
    nominal_sampling_frequency =
      calin_run_config->nectarcam().module(0).nominal_sampling_frequency();
    for(int imod=1; imod<calin_run_config->nectarcam().module_size(); imod++)
      if(calin_run_config->nectarcam().module(imod).nominal_sampling_frequency()
        != nominal_sampling_frequency)nominal_sampling_frequency = 0;
  }
  calin_run_config->set_nominal_sampling_frequency(nominal_sampling_frequency);

  // ==========================================================================
  //
  // RUN START TIME
  //
  // ==========================================================================

  if(cta_event and cta_event->has_nectarcam()
    and cta_event->nectarcam().has_cdts_data()
    and cta_event->nectarcam().cdts_data().has_data())
  {
    calin::ix::iact_data::telescope_event::CDTSData calin_cdts_data;
    calin::iact_data::acada_event_decoder::decode_cdts_data(
      &calin_cdts_data, cta_event->nectarcam().cdts_data());

    if(calin_cdts_data.event_counter() == cta_event->tel_event_id()) {
      run_start_time_ = calin_cdts_data.ucts_timestamp();
      calin_run_config->mutable_run_start_time()->set_time_ns(run_start_time_);
    }
  }

  // ==========================================================================
  //
  // SERIALIZE RAW DATA
  //
  // ==========================================================================

  if(cta_run_header and config_.include_serialized_raw_data())
  {
    calin_run_config->set_serialized_raw_header_type(
      SerializedRawHeaderType::SERIALIZED_RAW_HEADER_ACADA_PROTOBUF_R1V0);
    cta_run_header->SerializeToString(calin_run_config->mutable_serialized_raw_header());
  } else {
    calin_run_config->set_serialized_raw_header_type(
      SerializedRawHeaderType::SERIALIZED_RAW_HEADER_NONE);
  }
  return true;
}

void NectarCam_ACADACameraEventDecoder_R1v0::
copy_single_gain_integrals(const event_type* cta_event,
  const calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
  const int16_t* cta_charges, const uint8_t* cta_pixel_mask,
  calin::ix::iact_data::telescope_event::DigitizedSkyImage* calin_image,
  const std::string& which_gain) const
{
}

void NectarCam_ACADACameraEventDecoder_R1v0::
copy_single_gain_waveforms(
  const calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
  const int16_t* cta_waveforms, const uint8_t* cta_pixel_mask,
  calin::ix::iact_data::telescope_event::Waveforms* calin_waveforms,
  uint8_t has_gain_mask, const std::string& which_gain) const
{
  unsigned npix = nmod_*7;
  bool all_channels_present = true;

  calin_waveforms->mutable_channel_index()->Reserve(npix);
  calin_waveforms->mutable_channel_id()->Reserve(npix);

  std::string* calin_wf_raw_data_string = calin_waveforms->mutable_raw_samples_array();
  unsigned simd_vec_size = calin::provenance::system_info::the_host_info()->simd_vec_size();
  calin_wf_raw_data_string->resize(nsample_*npix*sizeof(int16_t) + simd_vec_size);
  char* cp = &calin_wf_raw_data_string->front();
  std::fill(cp+nsample_*npix*sizeof(int16_t), cp+calin_wf_raw_data_string->size(), int8_t(0));
  int16_t* calin_wf_raw_data = reinterpret_cast<int16_t*>(cp);

  for(unsigned ipix=0;ipix<npix;ipix++)
  {
    if(cta_pixel_mask[ipix] & has_gain_mask) {
      std::copy(cta_waveforms, cta_waveforms+nsample_, calin_wf_raw_data);
      calin_wf_raw_data += nsample_;
      calin_waveforms->add_channel_index(calin_waveforms->channel_id_size());
      if((has_gain_mask == 0x04) or (has_gain_mask == 0x0C and cta_pixel_mask[ipix] == 0x04)) {
        calin_waveforms->add_channel_signal_type(calin::ix::iact_data::telescope_event::SIGNAL_HIGH_GAIN);
      } else if ((has_gain_mask == 0x08) or (has_gain_mask == 0x0C and cta_pixel_mask[ipix] == 0x08)) {
        calin_waveforms->add_channel_signal_type(calin::ix::iact_data::telescope_event::SIGNAL_LOW_GAIN);
      } else {
        throw std::runtime_error("NectarCam_ACADACameraEventDecoder_R1v0::copy_single_gain_waveforms: Unhandled pixel mask: " +
          std::to_string(unsigned(has_gain_mask)) + " / " + std::to_string(unsigned(cta_pixel_mask[ipix])));
      }
      calin_waveforms->add_channel_id(ipix);
      if(config_.separate_channel_waveforms()) {
        auto* calin_samp = calin_waveforms->add_waveform()->mutable_samples();
        calin_samp->Reserve(nsample_);
        for(unsigned isample=0;isample<nsample_;isample++)
          calin_samp->Add(*cta_waveforms++);
      } else {
        cta_waveforms += nsample_;
      }
    } else {
      std::fill(calin_wf_raw_data, calin_wf_raw_data+nsample_, 0);
      calin_wf_raw_data += nsample_;
      all_channels_present = false;
      calin_waveforms->add_channel_index(-1);
      calin_waveforms->add_channel_signal_type(
        calin::ix::iact_data::telescope_event::SIGNAL_NONE);
      cta_waveforms += nsample_;
    }
  }

  calin_waveforms->set_all_channels_present(all_channels_present);
}

NectarCam_ACADACameraEventDecoder_R1v0* NectarCam_ACADACameraEventDecoder_R1v0::clone() const {
  return new NectarCam_ACADACameraEventDecoder_R1v0(*this);
}
