/*

   calin/iact_data/nectarcam_data_source.cpp -- Stephen Fegan -- 2024-09-06

   A supplier of single telescope data from ACADA ZFits data files

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
#include <iact_data/zfits_data_source.hpp>

using namespace calin::iact_data::acada_data_source;
using namespace calin::iact_data::zfits_data_source;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::util::log;
using calin::util::file::is_file;
using calin::util::file::is_readable;
using calin::util::file::expand_filename;

namespace { // anonymous
  template<typename Message> inline void delete_message(Message* message)
  {
    delete message;
  }

  template<> inline void delete_message<const void>(const void*)
  {
    // nothing to see here
  }
} // anonymous namespace

// =============================================================================
// ZFITSDataSource - chained ZFits files with decoder
// =============================================================================

template<typename MessageSet>
ZFITSDataSource<MessageSet>::
ZFITSDataSource(
    calin::iact_data::zfits_acada_data_source::
      ZFITSACADACameraEventDataSource<MessageSet>* acada_zfits,
    calin::iact_data::acada_event_decoder::
      ACADACameraEventDecoder<MessageSet>* decoder,
    bool adopt_acada_zfits, bool adopt_decoder):
  TelescopeRandomAccessDataSourceWithRunConfig(),
  decoder_(decoder), adopt_decoder_(adopt_decoder),
  acada_zfits_(acada_zfits),  adopt_acada_zfits_(adopt_acada_zfits)
{
  MessageSet acada_message_set;
  try {
    acada_zfits_->set_next_index(0);
    uint64_t unused_seq_index = 0;
    acada_message_set.event = acada_zfits_->borrow_next_event(unused_seq_index);
  } catch(...) {
    // ignore errors that occur reading sample event;
  }
  try {
    acada_message_set.header = acada_zfits_->get_run_header();
  } catch(...) {
    // ignore errors that occur reading run header
  }
  try {
    acada_message_set.data_stream = acada_zfits_->get_data_stream();
  } catch(...) {
    // ignore errors that occur reading run header
  }
  run_config_ = new TelescopeRunConfiguration;
  decoder_->decode_run_config(run_config_, acada_message_set);
  delete acada_message_set.header;
  delete_message(acada_message_set.data_stream);
  if(acada_message_set.event) {
    acada_zfits_->release_borrowed_event(acada_message_set.event);
  }
  run_config_->clear_fragment_filename();
  for(const auto& ffn : acada_zfits_->all_fragment_names()) {
    run_config_->add_fragment_filename(ffn);
  }
  run_config_->set_file_size(calin::util::file::total_size(
    acada_zfits_->all_fragment_names()));
  run_config_->set_num_events(acada_zfits->size());
  acada_zfits_->set_next_index(0);
}

template<typename MessageSet>
ZFITSDataSource<MessageSet>::
ZFITSDataSource(const std::string& filename,
    calin::iact_data::acada_event_decoder::
      ACADACameraEventDecoder<MessageSet>* decoder,
    bool adopt_decoder, const config_type& config):
  ZFITSDataSource(new calin::iact_data::zfits_acada_data_source::
      ZFITSACADACameraEventDataSource<MessageSet>(filename, config),
    decoder, /* adopt_acada_zfits= */ true, adopt_decoder)
{
  // nothing to see here
}

template<typename MessageSet>
ZFITSDataSource<MessageSet>::
ZFITSDataSource(const std::string& filename,
    ZFITSDataSource<MessageSet>* base_datasource, 
    const config_type& config):
  TelescopeRandomAccessDataSourceWithRunConfig(),
  decoder_(base_datasource->decoder_->clone()), adopt_decoder_(true),
  acada_zfits_(new calin::iact_data::zfits_acada_data_source::
    ZFITSACADACameraEventDataSource<MessageSet>(filename, config)), 
  adopt_acada_zfits_(true),
  run_config_(base_datasource->get_run_configuration())
{
  // nothing to see here
}

template<typename MessageSet>
ZFITSDataSource<MessageSet>::~ZFITSDataSource()
{
  delete run_config_;
  if(adopt_acada_zfits_)delete acada_zfits_;
  if(adopt_decoder_)delete decoder_;
}

template<typename MessageSet>
calin::ix::iact_data::telescope_event::TelescopeEvent* 
ZFITSDataSource<MessageSet>::get_next(
  uint64_t& seq_index_out, google::protobuf::Arena** arena)
{
  MessageSet acada_message_set;
  acada_message_set.event =
    acada_zfits_->borrow_next_event(seq_index_out);
  if(!acada_message_set.event){
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
    acada_zfits_->release_borrowed_event(acada_message_set.event);
    throw std::runtime_error("Could not allocate telescope event");
  }
  if(!decoder_->decode(event, acada_message_set))
  {
    delete delete_arena;
    delete delete_event;
    acada_zfits_->release_borrowed_event(acada_message_set.event);
    throw std::runtime_error("Could not decode ACTL event");
  }
  acada_zfits_->release_borrowed_event(acada_message_set.event);
  event->set_source_event_index(seq_index_out);
  return event;
}

template<typename MessageSet>
uint64_t ZFITSDataSource<MessageSet>::size()
{
  return acada_zfits_->size();
}

template<typename MessageSet>
void ZFITSDataSource<MessageSet>::
set_next_index(uint64_t next_index)
{
  acada_zfits_->set_next_index(next_index);
}

template<typename MessageSet>
calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* 
ZFITSDataSource<MessageSet>::get_run_configuration()
{
  if(!run_config_)return nullptr;
  auto* run_config = new TelescopeRunConfiguration();
  run_config->CopyFrom(*run_config_);
  return run_config;
}

template<typename MessageSet>
unsigned ZFITSDataSource<MessageSet>::
current_fragment_index() const
{
  return acada_zfits_->current_fragment_index();
}

template<typename MessageSet>
unsigned ZFITSDataSource<MessageSet>::
num_fragments() const
{
  return acada_zfits_->num_fragments();
}

template<typename MessageSet>
std::string ZFITSDataSource<MessageSet>::
fragment_name(unsigned index) const
{
  return acada_zfits_->fragment_name(index);
}

template<typename MessageSet>
typename ZFITSDataSource<MessageSet>::config_type ZFITSDataSource<MessageSet>::
default_config()
{
  return calin::iact_data::zfits_acada_data_source::
    ZFITSACADACameraEventDataSource<MessageSet>::default_config();
}

namespace calin { namespace iact_data { namespace zfits_data_source {

// template class ZFITSSingleFileDataSource<ACADA_EventMessage_L0, ACADA_HeaderMessage_L0>;
template class ZFITSDataSource<ACADA_MessageSet_L0>;

// template class ZFITSSingleFileDataSource<ACADA_EventMessage_R1v0, ACADA_HeaderMessage_R1v0>;
template class ZFITSDataSource<ACADA_MessageSet_R1v0>;

// template class ZFITSSingleFileDataSource<ACADA_EventMessage_R1v1, ACADA_HeaderMessage_R1v1>;
template class ZFITSDataSource<ACADA_MessageSet_R1v1>;

} } } // namespace calin::iact_data::zfits_data_source
