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

using namespace calin::iact_data::zfits_acada_data_source;
using namespace calin::iact_data::acada_data_source;
using namespace calin::util::log;
using calin::util::file::is_file;
using calin::util::file::is_readable;
using calin::util::file::expand_filename;

namespace {
  template<typename Message> calin::ix::iact_data::zfits_data_source::ACADADataModel 
  default_data_model() { 
    return calin::ix::iact_data::zfits_data_source::ACADA_DATA_MODEL_AUTO_DETECT; }
  template<> calin::ix::iact_data::zfits_data_source::ACADADataModel 
  default_data_model<ACADA_EventMessage_L0>() {
    return calin::ix::iact_data::zfits_data_source::ACADA_DATA_MODEL_L0; }
  template<> calin::ix::iact_data::zfits_data_source::ACADADataModel 
  default_data_model<ACADA_EventMessage_R1v0>() {
    return calin::ix::iact_data::zfits_data_source::ACADA_DATA_MODEL_R1V0; }
  template<> calin::ix::iact_data::zfits_data_source::ACADADataModel 
  default_data_model<ACADA_EventMessage_R1v1>() {
    return calin::ix::iact_data::zfits_data_source::ACADA_DATA_MODEL_R1V1; }

  template<typename Message> std::string default_message_table_name() { return "unknown"; }
  template<> std::string default_message_table_name<ACADA_HeaderMessage_L0>() { return "RunHeader"; }
  template<> std::string default_message_table_name<ACADA_EventMessage_L0>() { return "Events"; }
  template<> std::string default_message_table_name<ACADA_HeaderMessage_R1v0>() { return "CameraConfig"; }
  template<> std::string default_message_table_name<ACADA_EventMessage_R1v0>() { return "Events"; }
  template<> std::string default_message_table_name<ACADA_HeaderMessage_R1v1>() { return "CameraConfiguration"; }
  template<> std::string default_message_table_name<ACADA_EventMessage_R1v1>() { return "Events"; }
  template<> std::string default_message_table_name<ACADA_DataStreamMessage_R1v1>() { return "DataStream"; }

  template<typename Message> Message* read_one_message_from_zfits_table(
    const std::string& filename, const std::string& tablename, bool suppress_file_record = false) 
  {
    ZFITSSingleFileSingleMessageDataSource<Message> zfits(filename, tablename, suppress_file_record);
    uint64_t seq_index_out;
    return zfits.get_next(seq_index_out);
  }

  template<> void* read_one_message_from_zfits_table<void>(
    const std::string& filename, const std::string& tablename, bool suppress_file_record) 
  {
    return nullptr;
  }

  template<typename Message> Message* copy_message(const Message* message)
  {
    if(message == nullptr)return nullptr;
    Message* new_message = new Message;
    new_message->CopyFrom(*message);
    return new_message;
  }

  template<> void* copy_message<void>(const void*)
  {
    return nullptr;
  }

  template<typename Message> void delete_message(Message* message)
  {
    delete message;
  }

  template<> void delete_message<void>(void*)
  {
    // nothing to see here
  }

  template<> void delete_message<const void>(const void*)
  {
    // nothing to see here
  }

} // anonymous namespace

std::vector<std::string> calin::iact_data::zfits_acada_data_source::
get_zfits_table_names(std::string filename)
{
  calin::util::file::expand_filename_in_place(filename);
  std::vector<std::string> tables;
  try {
    IFits ifits(filename, "", /* force= */ true);
    if(ifits) {
      tables.push_back(ifits.GetTable().name);
    }
    while(ifits.hasNextTable()) {
      ifits.openNextTable(/* force= */ true);
      tables.push_back(ifits.GetTable().name);
    }
  } catch(...) {
    throw std::runtime_error(std::string("Could not open ZFits file: ")+filename);
  }
  return tables;
}

std::vector<std::string> calin::iact_data::zfits_acada_data_source::
get_zfits_table_column_names(std::string filename, std::string tablename)
{
  calin::util::file::expand_filename_in_place(filename);
  std::vector<std::string> columns;
  try {
    IFits ifits(filename, tablename, /* force= */ true);
    if(ifits) {
      for(const auto& col : ifits.GetColumns()) {
        columns.push_back(col.first);
      }
    }
  } catch(...) {
    throw std::runtime_error(std::string("Could not open ZFits file table: ")+filename + " -> " + tablename);
  }
  return columns;
}

std::vector<std::string> calin::iact_data::zfits_acada_data_source::
get_zfits_table_keys(std::string filename, std::string tablename)
{
  calin::util::file::expand_filename_in_place(filename);
  std::vector<std::string> keys;
  try {
    IFits ifits(filename, tablename, /* force= */ true);
    if(ifits) {
      for(const auto& key : ifits.GetKeys()) {
        keys.push_back(key.first);
      }
    }
  } catch(...) {
    throw std::runtime_error(std::string("Could not open ZFits file table: ")+filename + " -> " + tablename);
  }
  return keys;
}

std::map<std::string,std::string> calin::iact_data::zfits_acada_data_source::
get_zfits_table_key_values(std::string filename, std::string tablename)
{
  calin::util::file::expand_filename_in_place(filename);
  std::map<std::string,std::string> key_values;
  try {
    IFits ifits(filename, tablename, /* force= */ true);
    if(ifits) {
      for(const auto& key : ifits.GetKeys()) {
        key_values[key.first] = key.second.value;
      }
    }
  } catch(...) {
    throw std::runtime_error(std::string("Could not open ZFits file table: ")+filename + " -> " + tablename);
  }
  return key_values;
}

// =============================================================================
// =============================================================================
// =============================================================================
//
// ZFITSSingleFileSingleMessageDataSource - single ZFits file source of 
// single message type (for debugging purposes)
// 
// =============================================================================
// =============================================================================
// =============================================================================

template<typename Message> ZFITSSingleFileSingleMessageDataSource<Message>::
ZFITSSingleFileSingleMessageDataSource(
    const std::string& filename, const std::string& tablename, bool suppress_file_record):
  calin::io::data_source::RandomAccessDataSource<Message>(),
  filename_(expand_filename(filename)), 
  tablename_(tablename.empty() ? default_message_table_name<Message>() : tablename)
{
  if(!is_file(filename_))
    throw std::runtime_error(std::string("No such file: ")+filename_);
  if(!is_readable(filename_))
    throw std::runtime_error(std::string("File not readable: ")+filename_);

  zfits_ = new ADH::IO::ProtobufIFits(filename_, tablename_, Message::descriptor());
  if(zfits_->eof() && !zfits_->bad())
    throw std::runtime_error("ZFits file " + filename_ + " has no table: " + tablename_);

  if(not suppress_file_record) {
    file_record_ = calin::provenance::chronicle::register_file_open(filename_,
      calin::ix::provenance::chronicle::AT_READ, __PRETTY_FUNCTION__);
  }
}

template<typename Message> ZFITSSingleFileSingleMessageDataSource<Message>::
~ZFITSSingleFileSingleMessageDataSource()
{
  delete zfits_;
  calin::provenance::chronicle::register_file_close(file_record_);
}

template<typename Message> const Message* 
ZFITSSingleFileSingleMessageDataSource<Message>::
borrow_next_message(uint64_t& seq_index_out)
{
  uint64_t max_seq_index = zfits_->getNumMessagesInTable();
  if(next_message_index_ >= max_seq_index)return nullptr;

  seq_index_out = next_message_index_;
  const Message* message {
    zfits_->borrowTypedMessage<Message>(++next_message_index_) };
  if(message == nullptr)throw std::runtime_error("ZFits reader returned NULL");

  return message;
}
  
template<typename Message> void 
ZFITSSingleFileSingleMessageDataSource<Message>::
release_borrowed_message(const message_type* message)
{
  zfits_->returnBorrowedMessage(message);
}

template<typename Message> Message* 
ZFITSSingleFileSingleMessageDataSource<Message>::get_next(uint64_t& seq_index_out,
  google::protobuf::Arena** arena)
{
  if(arena)*arena = nullptr;

  uint64_t max_seq_index = zfits_->getNumMessagesInTable();
  if(next_message_index_ >= max_seq_index)return nullptr;

  seq_index_out = next_message_index_;
  const Message* message {
    zfits_->borrowTypedMessage<Message>(++next_message_index_) };
  if(message == nullptr)throw std::runtime_error("ZFits reader returned NULL");

  Message* message_copy = nullptr;
  if(arena && *arena)
    throw std::runtime_error("ZFITSSingleFileSingleMessageDataSource::get_next: "
      " pre-allocated arena not supported."); // ADH library does not support arenas
  
  message_copy = new Message;
  message_copy->CopyFrom(*message);

  zfits_->returnBorrowedMessage(message);

  return message_copy;
}

template<typename Message> uint64_t
ZFITSSingleFileSingleMessageDataSource<Message>::size()
{
  return zfits_->getNumMessagesInTable();
}

template<typename Message> void
ZFITSSingleFileSingleMessageDataSource<Message>::
set_next_index(uint64_t next_index)
{
  uint64_t max_seq_index = zfits_->getNumMessagesInTable();
  next_message_index_ = std::min(next_index, max_seq_index);
}

// =============================================================================
// =============================================================================
// =============================================================================
//
// ZFITSSingleFileACADACameraEventDataSource - single ZFits file
// 
// =============================================================================
// =============================================================================
// =============================================================================

template<typename MessageSet>
ZFITSSingleFileACADACameraEventDataSource<MessageSet>::
ZFITSSingleFileACADACameraEventDataSource(const std::string& filename, const config_type& config):
  ACADACameraEventRandomAccessDataSourceWithRunHeader<MessageSet>(),
  filename_(expand_filename(filename)), config_(config)
{
  if(config_.data_stream_table_name().empty())
    config_.set_data_stream_table_name(default_message_table_name<data_stream_type>());
  if(config_.run_header_table_name().empty())
    config_.set_run_header_table_name(default_message_table_name<header_type>());
  if(config_.events_table_name().empty())
    config_.set_events_table_name(default_message_table_name<event_type>());

  if(!is_file(filename_))
    throw std::runtime_error(std::string("No such file: ")+filename_);
  if(!is_readable(filename_))
    throw std::runtime_error(std::string("File not readable: ")+filename_);

  if(!config_.dont_read_run_header())
  {
    try {
      run_header_ = read_one_message_from_zfits_table<header_type>(
        filename_, config_.run_header_table_name(), /* suppress_file_record = */true);
    } catch(...) {
      if(!config_.ignore_run_header_errors())
        LOG(WARNING)
          << "ZFITSSingleFileACADACameraEventDataSource: Could not read run header from "
          << filename_;
    }

    try {
      data_stream_ = read_one_message_from_zfits_table<data_stream_type>(
        filename_, config_.data_stream_table_name(), /* suppress_file_record = */true);
    } catch(...) {
      if(!config_.ignore_run_header_errors())
        LOG(WARNING)
          << "ZFITSSingleFileACADACameraEventDataSource: Could not read data stream from "
          << filename_;
    }
  }

  zfits_ = new ZFITSSingleFileSingleMessageDataSource<event_type>(filename_, config_.events_table_name());
}

template<typename MessageSet>
ZFITSSingleFileACADACameraEventDataSource<MessageSet>::
~ZFITSSingleFileACADACameraEventDataSource()
{
  delete zfits_;
  delete run_header_;
  delete_message(data_stream_);
}

template<typename MessageSet> const typename MessageSet::event_type* 
ZFITSSingleFileACADACameraEventDataSource<MessageSet>::
borrow_next_event(uint64_t& seq_index_out)
{
  if(config_.max_seq_index() and zfits_->get_next_index()>=config_.max_seq_index()) {
    return nullptr;
  } else {
    return zfits_->borrow_next_message(seq_index_out);
  }
}

template<typename MessageSet>
void ZFITSSingleFileACADACameraEventDataSource<MessageSet>::
release_borrowed_event(const event_type* event)
{
  zfits_->release_borrowed_message(event);
}

template<typename MessageSet>
typename MessageSet::event_type* ZFITSSingleFileACADACameraEventDataSource<MessageSet>::
get_next(uint64_t& seq_index_out, google::protobuf::Arena** arena)
{
  if(config_.max_seq_index() and zfits_->get_next_index()>=config_.max_seq_index()) {
    return nullptr;
  } else {
    return zfits_->get_next(seq_index_out, arena);
  }
}

template<typename MessageSet>
uint64_t ZFITSSingleFileACADACameraEventDataSource<MessageSet>::size()
{
  uint64_t max_seq_index = zfits_->size();
  if(config_.max_seq_index())
    max_seq_index = std::min(max_seq_index, config_.max_seq_index());
  return max_seq_index;
}

template<typename MessageSet>
void ZFITSSingleFileACADACameraEventDataSource<MessageSet>::
set_next_index(uint64_t next_index)
{
  if(config_.max_seq_index()) {
    next_index = std::min(next_index, config_.max_seq_index());
  }
  zfits_->set_next_index(next_index);
}

template<typename MessageSet>
typename MessageSet::header_type* ZFITSSingleFileACADACameraEventDataSource<MessageSet>::
get_run_header()
{
  return copy_message(run_header_);
}

template<typename MessageSet>
typename MessageSet::data_stream_type* ZFITSSingleFileACADACameraEventDataSource<MessageSet>::
get_data_stream()
{
  return copy_message(data_stream_);
}


template<typename MessageSet>
typename ZFITSSingleFileACADACameraEventDataSource<MessageSet>::config_type
ZFITSSingleFileACADACameraEventDataSource<MessageSet>::default_config()
{
  config_type config = config_type::default_instance();
  config.set_data_model(default_data_model<event_type>());
  config.set_extension(".fits.fz");
  config.set_data_stream_table_name(default_message_table_name<data_stream_type>());
  config.set_run_header_table_name(default_message_table_name<header_type>());
  config.set_events_table_name(default_message_table_name<event_type>());
  config.set_file_fragment_stride(1);
  return config;
}

// =============================================================================
// =============================================================================
// =============================================================================
//
// ZFITSACADACameraEventDataSource - data source chaining multiple ZFits file
//
// =============================================================================
// =============================================================================
// =============================================================================

template<typename MessageSet>
ZFITSACADACameraEventDataSource<MessageSet>::
ZFITSACADACameraEventDataSource(const std::string& filename, const config_type& config):
  calin::io::data_source::BasicChainedRandomAccessDataSource<
    ACADACameraEventRandomAccessDataSourceWithRunHeader<MessageSet> >(
      new ZFITSACADACameraEventDataSourceOpener<MessageSet>(
        filename, config), true),
  config_(config), run_header_(nullptr)
{
  if(source_) {
    run_header_ = source_->get_run_header();
    data_stream_ = source_->get_data_stream();
  }
  if((run_header_ == nullptr or data_stream_ == nullptr) and opener_->num_sources() > 1) {
    // In this case the first file fragment is missing the RunHeader or DataStream, search
    // for it in later fragments then reopen the first fragment

    while((run_header_ == nullptr or data_stream_ == nullptr) and isource_ < opener_->num_sources())
    {
      ++isource_;
      open_file();
      if(run_header_ == nullptr)
        run_header_ = source_ ? source_->get_run_header() : nullptr;
      if(data_stream_ == nullptr)
        data_stream_ = source_ ? source_->get_data_stream() : nullptr;
    }

    isource_ = 0;
    open_file();
  }
}

template<typename MessageSet>
ZFITSACADACameraEventDataSource<MessageSet>::
~ZFITSACADACameraEventDataSource()
{
  delete run_header_;
  delete_message(data_stream_);
}

template<typename MessageSet>
typename MessageSet::header_type* ZFITSACADACameraEventDataSource<MessageSet>::
get_run_header()
{
  if(!run_header_)return nullptr;
  auto* run_header = new header_type();
  run_header->CopyFrom(*run_header_);
  return run_header;
}

template<typename MessageSet>
typename MessageSet::data_stream_type* ZFITSACADACameraEventDataSource<MessageSet>::
get_data_stream()
{
  return copy_message(data_stream_);
}

template<typename MessageSet>
const typename MessageSet::event_type* ZFITSACADACameraEventDataSource<MessageSet>::
borrow_next_event(uint64_t& seq_index_out)
{
  if(config_.max_seq_index() and seq_index_>=config_.max_seq_index())
    return nullptr;
  while(isource_ < opener_->num_sources())
  {
    uint64_t unused_index = 0;
    if(const event_type* next = source_->borrow_next_event(unused_index))
    {
      seq_index_out = seq_index_;
      ++seq_index_;
      return next;
    }
    ++isource_;
    open_file();
  }
  return nullptr;
}

template<typename MessageSet>
void ZFITSACADACameraEventDataSource<MessageSet>::
release_borrowed_event(const event_type* event)
{
  if(source_)
    source_->release_borrowed_event(event);
  else
    delete event;
}

template<typename MessageSet>
typename MessageSet::event_type * ZFITSACADACameraEventDataSource<MessageSet>::
get_next(uint64_t& seq_index_out, google::protobuf::Arena** arena)
{
  if(config_.max_seq_index() and seq_index_>=config_.max_seq_index())
    return nullptr;
  else
    return BaseDataSource::get_next(seq_index_out, arena);
}

template<typename MessageSet>
uint64_t ZFITSACADACameraEventDataSource<MessageSet>::size()
{
  uint64_t max_seq_index = BaseDataSource::size();
  if(config_.max_seq_index() and max_seq_index>config_.max_seq_index())
    max_seq_index = config_.max_seq_index();
  return max_seq_index;
}

template<typename MessageSet>
void ZFITSACADACameraEventDataSource<MessageSet>::
set_next_index(uint64_t next_index)
{
  if(config_.max_seq_index() and next_index>config_.max_seq_index())
    next_index = config_.max_seq_index();
  BaseDataSource::set_next_index(next_index);
}

template<typename MessageSet>
typename ZFITSACADACameraEventDataSource<MessageSet>::config_type
ZFITSACADACameraEventDataSource<MessageSet>::default_config()
{
  return ZFITSSingleFileACADACameraEventDataSource<MessageSet>::default_config();
}

// =============================================================================
// =============================================================================
// =============================================================================
//
// ZFITSACADACameraEventDataSourceOpener - open Zits file fragments
//
// =============================================================================
// =============================================================================
// =============================================================================

template<typename MessageSet>
ZFITSACADACameraEventDataSourceOpener<MessageSet>::
ZFITSACADACameraEventDataSourceOpener(std::string filename, const config_type& config):
  calin::io::data_source::DataSourceOpener<
    ACADACameraEventRandomAccessDataSourceWithRunHeader<MessageSet> >(),
  config_(config)
{
  const unsigned istride = std::max(1U,config.file_fragment_stride());
  filename = expand_filename(filename);
  if(is_file(filename))
    filenames_.emplace_back(filename);
  else
    throw(std::runtime_error("File not found: " + filename));

  if(not config_.exact_filename_only())
  {
    const std::string extension = config_.extension();
    auto ifind = filename.rfind(extension);
    if(ifind == filename.size()-extension.size())
    {
      filename = filename.substr(0, ifind);

      unsigned istart = 0;
      if(not is_file(filename+".1"+extension))
      {
        ifind = filename.rfind('.');
        if(ifind != std::string::npos and
          std::all_of(filename.begin() + ifind + 1, filename.end(), ::isdigit))
        {
          istart = std::stoi(filename.substr(ifind + 1));
          filename = filename.substr(0, ifind);
        }
      }

      bool fragment_found = true;
      for(unsigned i=istart+istride; fragment_found and (config_.max_file_fragments()==0 or
        filenames_.size()<config_.max_file_fragments()) ; i+=istride)
      {
        fragment_found = false;
        std::string fragment_i { std::to_string(i) };
        do {
          std::string filename_i { filename+"."+fragment_i+extension };
          if(is_file(filename_i)) {
            filenames_.emplace_back(filename_i);
            fragment_found = true;
          } else {
            fragment_i = std::string("0") + fragment_i;
          }
        }while(not fragment_found and fragment_i.size() <= 6);
      }
    }
  }
}

template<typename MessageSet>
ZFITSACADACameraEventDataSourceOpener<MessageSet>::
~ZFITSACADACameraEventDataSourceOpener()
{
    // nothing to see here
}

template<typename MessageSet>
unsigned ZFITSACADACameraEventDataSourceOpener<MessageSet>::num_sources() const
{
  return filenames_.size();
}

template<typename MessageSet>
std::string ZFITSACADACameraEventDataSourceOpener<MessageSet>::
source_name(unsigned isource) const
{
  if(isource >= filenames_.size())return {};
  return filenames_[isource];
}

template<typename MessageSet>
ZFITSSingleFileACADACameraEventDataSource<MessageSet>*
ZFITSACADACameraEventDataSourceOpener<MessageSet>::
open(unsigned isource)
{
  if(isource >= filenames_.size())return nullptr;
  if(config_.log_on_file_open())
    LOG(INFO) << "Opening file: " << filenames_[isource];
  auto config = config_;
  if(has_opened_file_)config.set_dont_read_run_header(true);
  config.set_max_seq_index(0);
  has_opened_file_ = true;
  return new ZFITSSingleFileACADACameraEventDataSource<MessageSet>(filenames_[isource], config);
}

template<typename MessageSet>
typename ZFITSACADACameraEventDataSourceOpener<MessageSet>::config_type
ZFITSACADACameraEventDataSourceOpener<MessageSet>::default_config()
{
  return ZFITSSingleFileACADACameraEventDataSource<MessageSet>::default_config();
}

namespace calin { namespace iact_data { namespace zfits_acada_data_source {

template class ZFITSSingleFileSingleMessageDataSource<ACADA_EventMessage_L0>;
template class ZFITSSingleFileSingleMessageDataSource<ACADA_HeaderMessage_L0>;
template class ZFITSSingleFileACADACameraEventDataSource<ACADA_MessageSet_L0>;
template class ZFITSACADACameraEventDataSource<ACADA_MessageSet_L0>;
template class ZFITSACADACameraEventDataSourceOpener<ACADA_MessageSet_L0>;

template class ZFITSSingleFileSingleMessageDataSource<ACADA_EventMessage_R1v0>;
template class ZFITSSingleFileSingleMessageDataSource<ACADA_HeaderMessage_R1v0>;
template class ZFITSSingleFileACADACameraEventDataSource<ACADA_MessageSet_R1v0>;
template class ZFITSACADACameraEventDataSource<ACADA_MessageSet_R1v0>;
template class ZFITSACADACameraEventDataSourceOpener<ACADA_MessageSet_R1v0>;

template class ZFITSSingleFileSingleMessageDataSource<ACADA_EventMessage_R1v1>;
template class ZFITSSingleFileSingleMessageDataSource<ACADA_HeaderMessage_R1v1>;
template class ZFITSSingleFileSingleMessageDataSource<ACADA_DataStreamMessage_R1v1>;
template class ZFITSSingleFileACADACameraEventDataSource<ACADA_MessageSet_R1v1>;
template class ZFITSACADACameraEventDataSource<ACADA_MessageSet_R1v1>;
template class ZFITSACADACameraEventDataSourceOpener<ACADA_MessageSet_R1v1>;

} } } // namespace calin::iact_data::zfits_acada_data_source
