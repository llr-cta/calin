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

ZFITSDataSource::ZFITSDataSource(const std::string& filename,
  CTACameraEventDecoder* decoder, bool adopt_decoder,
  const config_type& config):
  TelescopeRandomAccessDataSourceWithRunConfig(),
  decoder_(decoder), adopt_decoder_(adopt_decoder)
{
  actl_zfits_ = new zfits_actl_data_source::
    ZFITSACTLDataSource(filename, config);
  actl_zfits_->set_next_index(1);
  uint64_t unused_seq_index = 0;
  auto* sample_event = actl_zfits_->get_next(unused_seq_index);
  run_config_ = new TelescopeRunConfiguration;
  decoder_->decode_run_config(run_config_, actl_zfits_->get_run_header(),
    sample_event);
  delete sample_event;
  actl_zfits_->set_next_index(0);
}

ZFITSDataSource::~ZFITSDataSource()
{
  delete run_config_;
  delete actl_zfits_;
  if(adopt_decoder_)delete decoder_;
}

calin::ix::iact_data::telescope_event::TelescopeEvent*
ZFITSDataSource::get_next(uint64_t& seq_index_out,
  google::protobuf::Arena** arena)
{
  std::unique_ptr<DataModel::CameraEvent> cta_event {
    actl_zfits_->get_next(seq_index_out, arena) };
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
    throw std::runtime_error("Could not allocate telescpe event");
  }
  if(!decoder_->decode(event, cta_event.get()))
  {
    delete delete_arena;
    delete delete_event;
    throw std::runtime_error("Could not decode ACTL event");
  }
  event->set_source_event_index(seq_index_out);
  return event;
}

calin::ix::iact_data::telescope_run_configuration::
TelescopeRunConfiguration* ZFITSDataSource::get_run_configuration()
{
  if(!run_config_)return nullptr;
  auto* run_config = new TelescopeRunConfiguration();
  run_config->CopyFrom(*run_config_);
  return run_config;
}

uint64_t ZFITSDataSource::size()
{
  return actl_zfits_->size();
}

void ZFITSDataSource::set_next_index(uint64_t next_index)
{
  return actl_zfits_->set_next_index(next_index);
}
