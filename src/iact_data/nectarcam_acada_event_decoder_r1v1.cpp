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

NectarCam_ACADACameraEventDecoder_R1v1::
NectarCam_ACADACameraEventDecoder_R1v1(const std::string& filename, 
    const config_type& config):
  calin::iact_data::unified_acada_event_decoder::Unified_ACADACameraEventDecoder_R1v1(
    filename, unified_decoder_config(config)),
  config_(config)
{
  // nothing to see here
}

NectarCam_ACADACameraEventDecoder_R1v1::
~NectarCam_ACADACameraEventDecoder_R1v1()
{
  // nothing to see here
}

bool NectarCam_ACADACameraEventDecoder_R1v1::decode(
  calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
  const calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1& cta_messages)
{
  if(not Unified_ACADACameraEventDecoder_R1v1::decode(calin_event, cta_messages)) {
    return false;
  }

  const event_type* cta_event = cta_messages.event;

  bool all_modules_present = true;
  bool mod_sum_clock_set = false;
  if(cta_event->has_debug()) {

    if(cta_event->debug().has_module_status() and cta_event->debug().module_status().has_data()) {
      // **************************************************************************
      // Decode NectarCAM module status (presence)
      // **************************************************************************

      const auto& cta_status = cta_event->debug().module_status();

      calin_event->mutable_module_index()->Resize(nmod_configured_,-1);
      calin_event->mutable_module_id()->Resize(nmod_configured_,-1);

      if(nmod_configured_ != cta_status.data().size())
        throw std::runtime_error("NectarCam_ACADACameraEventDecoder_R1v1::decode: "
          "Module status array size does not match number of modules (" + 
          std::to_string(cta_status.data().size()) + " != " + 
          std::to_string(nmod_configured_) + ")");
    
      const auto* mod_status = reinterpret_cast<const uint8_t*>(&cta_status.data().front());
      unsigned mod_index=0;
      for(unsigned imod=0;imod<nmod_configured_;imod++)
      {
        if(*(mod_status++)&0x01)
        {
          calin_event->set_module_index(imod,mod_index);
          calin_event->set_module_id(mod_index,imod);
          mod_index++;
        }
        else
        {
          calin_event->set_module_index(imod,-1);
          all_modules_present = false;
        }
      }
      calin_event->mutable_module_id()->Resize(mod_index,-1);
    } else {
      throw(std::runtime_error("NectarCam_ACADACameraEventDecoder_R1v1::decode: "
        "ACADA event does not have NectarCAM module_status"));
    }
    calin_event->set_all_modules_present(all_modules_present);
    
    if(cta_event->debug().has_counters() and cta_event->debug().counters().has_data()) {
      // **************************************************************************
      // Decode NectarCAM module counter data
      // **************************************************************************

      uint64_t mod_clock_sum = 0;
      uint64_t mod_clock_seq_sum = 0;
      int64_t mod_clock_num = 0;
      bool mod_clock_is_suspect = false;

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

      const auto& cta_counters = cta_event->debug().counters();

      const unsigned raw_data_version = 1;
      const unsigned raw_data_size = sizeof(NectarCountersWithTriggerPattern);

      if(cta_counters.data().size() != nmod_configured_ * raw_data_size) {
        throw std::runtime_error("Camera counters data array correct size. "
          "Expected " +
            std::to_string(nmod_configured_*raw_data_size) +
            ", got " + std::to_string(cta_counters.data().size()));
      }

      auto* trigger_map = calin_event->mutable_trigger_map();
      trigger_map->mutable_hit_channel_id()->Reserve(20);
      trigger_map->mutable_trigger_image()->Resize(nmod_configured_ * 7, 0x80000000);

      const auto* mod_data_string = cta_counters.data().data();
      for(unsigned imod=0;imod<nmod_configured_;imod++, mod_data_string+=raw_data_size)
      {
        const auto* mod_data_struct = reinterpret_cast<const NectarCountersWithTriggerPattern*>(mod_data_string);

        if(imod >= static_cast<unsigned>(calin_event->module_index_size()) or
          calin_event->module_index(imod) == -1)continue;

        auto* module_counters = calin_event->add_module_counter();
        module_counters->set_module_id(imod);
#define add_mod_counter(id,val) \
        { \
          module_counters->add_counter_id(id); \
          module_counters->add_counter_value(val); \
        }
        add_mod_counter(0, mod_data_struct->global_event_counter);
        add_mod_counter(1, mod_data_struct->bunch_counter);
        add_mod_counter(2, mod_data_struct->event_counter);
        add_mod_counter(3, mod_data_struct->ts1);
        add_mod_counter(4, mod_data_struct->ts2_bunch);
        add_mod_counter(5, mod_data_struct->ts2_event);
        add_mod_counter(6, mod_data_struct->ts2_empty);

        auto* module_data = calin_event->add_module_data()->mutable_nectarcam();
        module_data->set_module_id(imod);
        module_data->set_global_event_counter(mod_data_struct->global_event_counter);
        module_data->set_bunch_counter(mod_data_struct->bunch_counter);
        module_data->set_event_counter(mod_data_struct->event_counter);
        module_data->set_ts1(mod_data_struct->ts1);
        module_data->set_ts2_bunch(mod_data_struct->ts2_bunch);
        module_data->set_ts2_event(mod_data_struct->ts2_event);
        module_data->set_ts2_empty(mod_data_struct->ts2_empty);
        module_data->set_version(raw_data_version);

        module_data->set_trigger_pattern(mod_data_struct->trigger_pattern & 0x7F7F7F7F);

        for(unsigned imodchan=0; imodchan<7; imodchan++) {
          // Shuffle trigger pattern for each channel into its own bit stream
          // As an example, if imodchan=4
          // and pattern starts as               XXX4XXXX XXX3XXXX XXX2XXXX XXX1XXXX
          // first apply shift and mask giving   00000004 00000003 00000002 00000001
          uint32_t chan_trigger = (mod_data_struct->trigger_pattern>>imodchan) & 0x01010101;
          chan_trigger |= chan_trigger >> 7;  // 00000004 00000043 00000032 00000021
          chan_trigger |= chan_trigger >> 14; // 00000004 00000043 00000432 00004321
          chan_trigger &= 0x0F;               // 00000000 00000000 00000000 00004321
          trigger_map->set_trigger_image(imod*7+imodchan, chan_trigger);
          if(chan_trigger) {
            trigger_map->add_hit_channel_id(imod*7+imodchan);
          }
        }

#define ts2_decode(x) int32_t(x)

        int32_t ts2_bunch = ts2_decode(mod_data_struct->ts2_bunch);
        int32_t ts2_event = ts2_decode(mod_data_struct->ts2_event);
        int32_t ts = mod_data_struct->ts1*8 + ts2_event - ts2_bunch;

        module_data->set_bunch_event_time(ts);

        bool clock_is_suspect = false;
        int32_t time_seq_id = mod_data_struct->bunch_counter;

        if(mod_data_struct->event_counter == 1 and mod_data_struct->ts1 > 124987500)
        {
          // Here we attempt to handle events where the TS1 value seems to
          // still be at the a high count value but the PPS event_counter
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
        clock->set_time_value(mod_data_struct->ts1);
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

      // **************************************************************************
      // Set the summed module clock
      // **************************************************************************

      if(mod_clock_num==nmod_configured_) {
        calin_event->set_camera_clock_index(9, calin_event->camera_clock_size());
        auto* calin_clock = calin_event->add_camera_clock();
        calin_clock->set_clock_id(9);
        calin_clock->set_time_value(mod_clock_sum);
        calin_clock->set_time_sequence_id(mod_clock_seq_sum);
        calin_clock->set_time_value_may_be_suspect(mod_clock_is_suspect);

        calin_event->set_camera_clock_index(10, calin_event->camera_clock_size());
        calin_clock = calin_event->add_camera_clock();
        calin_clock->set_clock_id(10);
        calin_clock->set_time_value(mod_clock_seq_sum*1000000000ULL + mod_clock_sum);
        calin_clock->set_time_sequence_id(0);
        calin_clock->set_time_value_may_be_suspect(mod_clock_is_suspect);

        mod_sum_clock_set = true;
      }
    }
  }

  if(not mod_sum_clock_set) {
    calin_event->set_camera_clock_index(9,-1);
    calin_event->set_camera_clock_index(10,-1);
  }

  return true;
}

bool NectarCam_ACADACameraEventDecoder_R1v1::decode_run_config(
  calin::ix::iact_data::telescope_run_configuration:: TelescopeRunConfiguration* calin_run_config,
  const calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1& cta_messages)
{
  if(not Unified_ACADACameraEventDecoder_R1v1::decode_run_config(calin_run_config, cta_messages)) {
    return false;
  }

  calin_run_config->set_data_transcoder(
    "calin::iact_data::unified_event_decoder::NectarCam_ACADACameraEventDecoder_R1v1");

  // const data_stream_type* cta_data_stream = cta_messages.data_stream;
  // const header_type* cta_run_header = cta_messages.header;
  // const event_type* cta_event = cta_messages.event;

  #define ADD_CAMERA_CLOCK(name, freq) \
    calin_run_config->mutable_camera_layout()->add_camera_clock_name(name); \
    calin_run_config->mutable_camera_layout()->add_camera_clock_frequency(freq)

  #define ADD_MODULE_CLOCK(name, freq) \
    calin_run_config->mutable_camera_layout()->add_module_clock_name(name); \
    calin_run_config->mutable_camera_layout()->add_module_clock_frequency(freq)

  ADD_CAMERA_CLOCK("FEB local 2ns TDC counter sum",                    // 9
    double(calin_run_config->configured_module_id_size())*1.0e9); 
  ADD_CAMERA_CLOCK("FEB local pps and 2ns TDC counter sum",            // 10
    double(calin_run_config->configured_module_id_size())*1.0e9); 

  ADD_MODULE_CLOCK("local 2ns TDC time",                    1.0e9);   // 0
  ADD_MODULE_CLOCK("local 125MHz oscillator counter",       1.25e8);  // 1
  ADD_MODULE_CLOCK("pps counter",                           1.0);     // 2
 
  ncamera_clock_ = calin_run_config->camera_layout().camera_clock_name_size();
  nmodule_clock_ = calin_run_config->camera_layout().module_clock_name_size();

  calin_run_config->mutable_camera_layout()->add_module_counter_name("global_event_counter");
  calin_run_config->mutable_camera_layout()->add_module_counter_name("bunch_counter");
  calin_run_config->mutable_camera_layout()->add_module_counter_name("event_counter");
  calin_run_config->mutable_camera_layout()->add_module_counter_name("ts1"); 
  calin_run_config->mutable_camera_layout()->add_module_counter_name("ts2_bunch");
  calin_run_config->mutable_camera_layout()->add_module_counter_name("ts2_event");
  calin_run_config->mutable_camera_layout()->add_module_counter_name("ts2_empty");

  ncamera_clock_ = calin_run_config->camera_layout().camera_clock_name_size();
  nmodule_clock_ = calin_run_config->camera_layout().module_clock_name_size();

  // **************************************************************************
  // Try to read the NectarCam module configuration XML file
  // **************************************************************************

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

  return true;
}

NectarCam_ACADACameraEventDecoder_R1v1* NectarCam_ACADACameraEventDecoder_R1v1::clone() const
{
  return new NectarCam_ACADACameraEventDecoder_R1v1(*this);
}

calin::ix::iact_data::cta_data_source::UnifiedCameraEventDecoderConfig 
NectarCam_ACADACameraEventDecoder_R1v1::unified_decoder_config(config_type config) {
  auto unified_config = 
    calin::iact_data::unified_acada_event_decoder::Unified_ACADACameraEventDecoder_R1v1::default_config();  
  unified_config.set_camera_type(calin::ix::iact_data::cta_data_source::NECTARCAM);
  unified_config.set_separate_channel_waveforms(config.separate_channel_waveforms());
  unified_config.set_include_serialized_raw_data(config.separate_channel_waveforms());
  return unified_config;
}
