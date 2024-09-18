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
#include <iact_data/unified_acada_event_decoder.hpp>
#include <provenance/system_info.hpp>
#include <iact_data/nectarcam_layout.hpp>

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
decode(calin::ix::iact_data::telescope_event::TelescopeEvent* event,
  const calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1& cta_messages)
{
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
  const event_type* cta_event = cta_messages.event;

  if(cta_run_header == nullptr) {
    throw std::runtime_error("Unified_ACADACameraEventDecoder_R1v1::decode_run_config: "
      "No camera configurtaion found");
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

  // Hard-code for NectarCAM for the moment - eventually use telescope id
  nectarcam_layout::nectarcam_layout(calin_run_config->mutable_camera_layout());

  // **************************************************************************
  // Build configured-module list
  // **************************************************************************

  nmod_configured_ = cta_run_header->num_modules();
  if(nmod_configured_ > 0) {
    if(not cta_run_header->has_module_id_map() or
      cta_run_header->module_id_map().data().size() != nmod_configured_*sizeof(uint16_t)) {
      throw std::runtime_error("bool Unified_ACADACameraEventDecoder_R1v1::decode_run_config: "
        "configured-module map inconsistant");
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

  return true;
}

Unified_ACADACameraEventDecoder_R1v1* Unified_ACADACameraEventDecoder_R1v1::clone() const
{
  return new Unified_ACADACameraEventDecoder_R1v1(*this);
}

