/*

   calin/iact_data/zfits_acada_data_source.cpp -- Stephen Fegan -- 2024-09-05

   Source of "raw" ACADA data types from ZFITS files

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
#include <provenance/chronicle.hpp>
#include <util/file.hpp>
#include <iact_data/zfits_acada_data_source.hpp>
#include <ProtobufIFits.h>
#include <L0.pb.h>

using namespace calin::iact_data::zfits_acada_data_source;
using namespace calin::iact_data::acada_data_source;
using namespace calin::util::log;
using calin::util::file::is_file;
using calin::util::file::is_readable;
using calin::util::file::expand_filename;

namespace {
  template<typename Message> std::string default_message_table_name() { return "unknown"; }
  template<> std::string default_message_table_name<ACADA_L0_HeaderMessage>() { return "RunHeader"; }
  template<> std::string default_message_table_name<ACADA_L0_EventMessage>() { return "Events"; }
} // anonymous namespace

template<typename EventMessage, typename HeaderMessage>
ZFITSSingleFileACADACameraEventDataSource<EventMessage,HeaderMessage>::
ZFITSSingleFileACADACameraEventDataSource(const std::string& filename, config_type config):
  ACADACameraEventRandomAccessDataSourceWithRunHeader<EventMessage,HeaderMessage>(),
  filename_(expand_filename(filename)), config_(config)
{
  if(config.run_header_table_name().empty())
    config.set_run_header_table_name(default_message_table_name<HeaderMessage>());
  if(config.events_table_name().empty())
    config.set_events_table_name(default_message_table_name<EventMessage>());

  if(!is_file(filename_))
    throw std::runtime_error(std::string("No such file: ")+filename_);
  if(!is_readable(filename_))
    throw std::runtime_error(std::string("File not readable: ")+filename_);

  if(!config.dont_read_run_header())
  {
    try
    {
      ACTL::IO::ProtobufIFits rh_zfits(filename_.c_str(),
        config.run_header_table_name(), HeaderMessage::descriptor());
      if(rh_zfits.eof() && !rh_zfits.bad())
        throw std::runtime_error("ZFits reader found no table:" +
          config.run_header_table_name());
      if(config.verify_file_after_open() or config.repair_broken_file())
      {
        try
        {
          rh_zfits.CheckIfFileIsConsistent(false);
        }
        catch (std::exception& e)
        {
          if(config.repair_broken_file())
          {
            LOG(WARNING) << "ZFits file " + filename_ +
              ": integrity verification failed, attempting to repair.";
            rh_zfits.CheckIfFileIsConsistent(true);
          }
        }
      }
      if(rh_zfits.getNumMessagesInTable() > 0)
      {
        HeaderMessage* run_header =
          rh_zfits.readTypedMessage<HeaderMessage>(1);
        run_header_ = new HeaderMessage(*run_header);
        rh_zfits.recycleMessage(run_header);
        //LOG(INFO) << run_header_->DebugString();
      }
    }
    catch(...)
    {
      if(!config.ignore_run_header_errors())
        LOG(WARNING)
          << "ZFITSSingleFileACADACameraEventDataSource: Could not read run header from "
          << filename_;
    }
  }

  zfits_ = new ACTL::IO::ProtobufIFits(filename_.c_str(),
    config.events_table_name(), EventMessage::descriptor());
  if(zfits_->eof() && !zfits_->bad())
    throw std::runtime_error("ZFits file " + filename_ + " has no table: " +
      config.events_table_name());

  file_record_ = calin::provenance::chronicle::register_file_open(filename_,
    calin::ix::provenance::chronicle::AT_READ, __PRETTY_FUNCTION__);

  if(config.verify_file_after_open() or config.repair_broken_file())
  {
    try
    {
      zfits_->CheckIfFileIsConsistent(false);
    }
    catch (std::exception& e)
    {
      if(config.repair_broken_file())
      {
        LOG(WARNING) << "ZFits file " + filename_ +
          ": integrity verification failed, attempting to repair.";
        zfits_->CheckIfFileIsConsistent(true);
      }
    }
  }
}

template<typename EventMessage, typename HeaderMessage>
ZFITSSingleFileACADACameraEventDataSource<EventMessage,HeaderMessage>::
~ZFITSSingleFileACADACameraEventDataSource()
{
  delete zfits_;
  calin::provenance::chronicle::register_file_close(file_record_);
  delete run_header_;
}

template<typename EventMessage, typename HeaderMessage> const EventMessage* 
ZFITSSingleFileACADACameraEventDataSource<EventMessage,HeaderMessage>::
borrow_next_event(uint64_t& seq_index_out)
{
  if(zfits_ == nullptr)
    throw std::runtime_error(std::string("File not open: ")+filename_);

  uint64_t max_seq_index = zfits_->getNumMessagesInTable();
  if(config_.max_seq_index())
    max_seq_index = std::min(max_seq_index, config_.max_seq_index());
  if(next_event_index_ >= max_seq_index)return nullptr;

  seq_index_out = next_event_index_;
  const EventMessage* event {
    zfits_->borrowTypedMessage<EventMessage>(++next_event_index_) };
  if(!event)throw std::runtime_error("ZFits reader returned NULL");

  return event;
}

template<typename EventMessage, typename HeaderMessage>
void ZFITSSingleFileACADACameraEventDataSource<EventMessage,HeaderMessage>::
release_borrowed_event(const EventMessage* event)
{
  zfits_->returnBorrowedMessage(event);
}

template<typename EventMessage, typename HeaderMessage>
EventMessage* ZFITSSingleFileACADACameraEventDataSource<EventMessage,HeaderMessage>::
get_next(uint64_t& seq_index_out, google::protobuf::Arena** arena)
{
  if(arena)*arena = nullptr;
  const EventMessage* event = borrow_next_event(seq_index_out);
  if(event == nullptr)return nullptr;

  EventMessage* event_copy = nullptr;
#if 0
  if(arena) {
    if(!*arena)*arena = new google::protobuf::Arena;
    event_copy =
      google::protobuf::Arena::CreateMessage<EventMessage>(*arena);
  }
  else event_copy = new EventMessage;
#else
  if(arena && *arena)
    throw std::runtime_error("ZFITSSingleFileACADACameraEventDataSource::get_next: "
      " pre-allocated arena not supported.");
  event_copy = new EventMessage;
#endif
  event_copy->CopyFrom(*event);
  release_borrowed_event(event);

  return event_copy;
}

template<typename EventMessage, typename HeaderMessage>
uint64_t ZFITSSingleFileACADACameraEventDataSource<EventMessage,HeaderMessage>::size()
{
  uint64_t max_seq_index = zfits_->getNumMessagesInTable();
  if(config_.max_seq_index())
    max_seq_index = std::min(max_seq_index, config_.max_seq_index());
  return max_seq_index;
}

template<typename EventMessage, typename HeaderMessage>
void ZFITSSingleFileACADACameraEventDataSource<EventMessage,HeaderMessage>::
set_next_index(uint64_t next_index)
{
  if(zfits_ == nullptr)next_event_index_ = 0;
  else {
    uint64_t max_seq_index = zfits_->getNumMessagesInTable();
    if(config_.max_seq_index())
      max_seq_index = std::min(max_seq_index, config_.max_seq_index());
    next_event_index_ = std::min(next_index, max_seq_index);
  }
}

template<typename EventMessage, typename HeaderMessage>
HeaderMessage* ZFITSSingleFileACADACameraEventDataSource<EventMessage,HeaderMessage>::
get_run_header()
{
  if(!run_header_)return nullptr;
  auto* run_header = new HeaderMessage;
  run_header->CopyFrom(*run_header_);
  return run_header;
}

template<typename EventMessage, typename HeaderMessage>
typename ZFITSSingleFileACADACameraEventDataSource<EventMessage,HeaderMessage>::config_type
ZFITSSingleFileACADACameraEventDataSource<EventMessage,HeaderMessage>::default_config()
{
  config_type config = config_type::default_instance();
  config.set_data_model(calin::ix::iact_data::zfits_data_source::ACTL_DATA_MODEL_L0);
  config.set_extension(".fits.fz");
  config.set_run_header_table_name(default_message_table_name<HeaderMessage>());
  config.set_events_table_name(default_message_table_name<EventMessage>());
  config.set_file_fragment_stride(1);
  return config;
}
namespace calin { namespace iact_data { namespace zfits_acada_data_source {

template class ZFITSSingleFileACADACameraEventDataSource<ACADA_L0_EventMessage, ACADA_L0_HeaderMessage>;
template class ZFITSSingleFileACADACameraEventDataSource<ACADA_R1v0_EventMessage, ACADA_R1v0_HeaderMessage>;

} } } // namespace calin::iact_data::zfits_acada_data_source
