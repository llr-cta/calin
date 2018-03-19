/*

   calin/iact_data/zfits_actl_data_source.cpp -- Stephen Fegan -- 2016-05-04

   A supplier of single telescope ACTL data from ZFits data files

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

#include <stdexcept>
#include <memory>
#include <cctype>

#include <util/log.hpp>
#include <provenance/chronicle.hpp>
#include <util/file.hpp>
#include <iact_data/zfits_actl_data_source.hpp>
#include <ProtobufIFits.h>
#include <L0.pb.h>

using namespace calin::iact_data::zfits_actl_data_source;
using namespace calin::util::log;
using calin::util::file::is_file;
using calin::util::file::is_readable;
using calin::util::file::expand_filename;

ACTLRandomAccessDataSourceWithRunHeader::
~ACTLRandomAccessDataSourceWithRunHeader()
{
  // nothing to see here
}

namespace {
  static std::string default_run_header_table_name("RunHeader");
  static std::string default_events_table_name("Events");
} // anonymous namespace

ZFITSSingleFileACTLDataSource::
ZFITSSingleFileACTLDataSource(const std::string& filename, config_type config):
  ACTLRandomAccessDataSourceWithRunHeader(),
  filename_(expand_filename(filename)), config_(config)
{
  if(config.run_header_table_name().empty())
    config.set_run_header_table_name(default_run_header_table_name);
  if(config.events_table_name().empty())
    config.set_events_table_name(default_events_table_name);

  if(!is_file(filename_))
    throw std::runtime_error(std::string("No such file: ")+filename_);
  if(!is_readable(filename_))
    throw std::runtime_error(std::string("File not readable: ")+filename_);

  if(!config.dont_read_run_header())
  {
    try
    {
      ACTL::IO::ProtobufIFits rh_zfits(filename_.c_str(),
        config.run_header_table_name(),
        DataModel::CameraRunHeader::descriptor());
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
        DataModel::CameraRunHeader* run_header =
          rh_zfits.readTypedMessage<DataModel::CameraRunHeader>(1);
        run_header_ = new DataModel::CameraRunHeader(*run_header);
        rh_zfits.recycleMessage(run_header);
        //LOG(INFO) << run_header_->DebugString();
      }
    }
    catch(...)
    {
      if(!config.ignore_run_header_errors())
        LOG(WARNING)
          << "ZFITSSingleFileACTLDataSource: Could not read run header from "
          << filename_;
    }
  }

  zfits_ = new ACTL::IO::ProtobufIFits(filename_.c_str(),
    config.events_table_name(), DataModel::CameraEvent::descriptor());
  if(zfits_->eof() && !zfits_->bad())
    throw std::runtime_error("ZFits file " + filename_ + " has no table: " +
      config.events_table_name());

  calin::provenance::chronicle::register_file_open(filename_,
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

ZFITSSingleFileACTLDataSource::~ZFITSSingleFileACTLDataSource()
{
  delete zfits_;
  delete run_header_;
}

ZFITSSingleFileACTLDataSource::config_type
ZFITSSingleFileACTLDataSource::default_config()
{
  config_type config = config_type::default_instance();
  config.set_extension(".fits.fz");
  config.set_run_header_table_name(default_run_header_table_name);
  config.set_events_table_name(default_events_table_name);
  return config;
}

const DataModel::CameraEvent*
ZFITSSingleFileACTLDataSource::borrow_next_event(uint64_t& seq_index_out)
{
  if(zfits_ == nullptr)
    throw std::runtime_error(std::string("File not open: ")+filename_);

  uint64_t max_seq_index = zfits_->getNumMessagesInTable();
  if(config_.max_seq_index())
    max_seq_index = std::min(max_seq_index, config_.max_seq_index());
  if(next_event_index_ >= max_seq_index)return nullptr;

  seq_index_out = next_event_index_;
  const DataModel::CameraEvent* event {
    zfits_->borrowTypedMessage<DataModel::CameraEvent>(++next_event_index_) };
  if(!event)throw std::runtime_error("ZFits reader returned NULL");

  return event;
}

void ZFITSSingleFileACTLDataSource::
release_borrowed_event(const DataModel::CameraEvent* event)
{
  zfits_->returnBorrowedMessage(event);
}

DataModel::CameraEvent* ZFITSSingleFileACTLDataSource::
get_next(uint64_t& seq_index_out, google::protobuf::Arena** arena)
{
  if(arena)*arena = nullptr;
  const DataModel::CameraEvent* event = borrow_next_event(seq_index_out);
  if(event == nullptr)return nullptr;

  DataModel::CameraEvent* event_copy = nullptr;
#if 0
  if(arena) {
    if(!*arena)*arena = new google::protobuf::Arena;
    event_copy =
      google::protobuf::Arena::CreateMessage<DataModel::CameraEvent>(*arena);
  }
  else event_copy = new DataModel::CameraEvent;
#else
  if(arena && *arena)
    throw std::runtime_error("ZFITSSingleFileACTLDataSource::get_next: "
      " pre-allocated arena not supported.");
  event_copy = new DataModel::CameraEvent;
#endif
  event_copy->CopyFrom(*event);
  release_borrowed_event(event);

  return event_copy;
}

uint64_t ZFITSSingleFileACTLDataSource::size()
{
  uint64_t max_seq_index = zfits_->getNumMessagesInTable();
  if(config_.max_seq_index())
    max_seq_index = std::min(max_seq_index, config_.max_seq_index());
  return max_seq_index;
}

void ZFITSSingleFileACTLDataSource::set_next_index(uint64_t next_index)
{
  if(zfits_ == nullptr)next_event_index_ = 0;
  else {
    uint64_t max_seq_index = zfits_->getNumMessagesInTable();
    if(config_.max_seq_index())
      max_seq_index = std::min(max_seq_index, config_.max_seq_index());
    next_event_index_ = std::min(next_index, max_seq_index);
  }
}

DataModel::CameraRunHeader* ZFITSSingleFileACTLDataSource::get_run_header()
{
  if(!run_header_)return nullptr;
  auto* run_header = new DataModel::CameraRunHeader();
  run_header->CopyFrom(*run_header_);
  return run_header;
}

ZFITSACTLDataSourceOpener::ZFITSACTLDataSourceOpener(std::string filename,
  const ZFITSACTLDataSource::config_type& config):
  calin::io::data_source::DataSourceOpener<
    ACTLRandomAccessDataSourceWithRunHeader>(),
  config_(config)
{
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
      if(not is_file(filename+".1"+extension)
        and not is_file(filename+".01"+extension)
        and not is_file(filename+".001"+extension)
        and not is_file(filename+".0001"+extension)
        and not is_file(filename+".00001"+extension))
      {
        ifind = filename.rfind('.');
        if(ifind != std::string::npos and
          std::all_of(filename.begin() + ifind + 1, filename.end(), ::isdigit))
        {
          istart = std::stoi(filename.substr(ifind + 1));
          filename = filename.substr(0, ifind);
        }
      }

      for(unsigned i=istart+1; config_.max_file_fragments()==0 or
        filenames_.size()<config_.max_file_fragments() ; ++i)
      {
        std::string filename_i { filename+"."+std::to_string(i)+extension };
        if(is_file(filename_i)) {
          filenames_.emplace_back(filename_i);
          continue;
        }

        // This is verbose but pretty simple
        if(i<10) {
          filename_i = filename + ".";
          filename_i += "0";
          filename_i += std::to_string(i)+extension;

          if(is_file(filename_i)) {
            filenames_.emplace_back(filename_i);
            continue;
          }
        }

        if(i<100) {
          filename_i = filename + ".";
          if(i<10) filename_i += "00";
          else filename_i += "0";
          filename_i += std::to_string(i)+extension;

          if(is_file(filename_i)) {
            filenames_.emplace_back(filename_i);
            continue;
          }
        }

        if(i<1000) {
          filename_i = filename + ".";
          if(i<10) filename_i += "000";
          else if(i<100) filename_i += "00";
          else filename_i += "0";
          filename_i += std::to_string(i)+extension;

          if(is_file(filename_i)) {
            filenames_.emplace_back(filename_i);
            continue;
          }
        }

        if(i<10000) {
          filename_i = filename + ".";
          if(i<10) filename_i += "0000";
          else if(i<100) filename_i += "000";
          else if(i<1000) filename_i += "00";
          else filename_i += "0";
          filename_i += std::to_string(i)+extension;

          if(is_file(filename_i)) {
            filenames_.emplace_back(filename_i);
            continue;
          }
        }

        break;
      }
    }
  }
}

ZFITSACTLDataSourceOpener::~ZFITSACTLDataSourceOpener()
{
  // nothing to see here
}

unsigned ZFITSACTLDataSourceOpener::num_sources()
{
  return filenames_.size();
}

calin::iact_data::zfits_actl_data_source::
ZFITSSingleFileACTLDataSource*
ZFITSACTLDataSourceOpener::open(unsigned isource)
{
  if(isource >= filenames_.size())return nullptr;
  if(config_.log_on_file_open())
    LOG(INFO) << "Opening file: " << filenames_[isource];
  auto config = config_;
  if(has_opened_file_)config.set_dont_read_run_header(true);
  config.set_max_seq_index(0);
  has_opened_file_ = true;
  return new ZFITSSingleFileACTLDataSource(filenames_[isource], config);
}

ZFITSACTLDataSource::ZFITSACTLDataSource(const std::string& filename,
  const config_type& config):
  calin::io::data_source::BasicChainedRandomAccessDataSource<
    ACTLRandomAccessDataSourceWithRunHeader>(
    new ZFITSACTLDataSourceOpener(filename, config), true),
  config_(config), run_header_(source_->get_run_header())
{
  // nothing to see here
}

ZFITSACTLDataSource::~ZFITSACTLDataSource()
{
  delete run_header_;
}

DataModel::CameraRunHeader* ZFITSACTLDataSource::get_run_header()
{
  if(!run_header_)return nullptr;
  auto* run_header = new DataModel::CameraRunHeader();
  run_header->CopyFrom(*run_header_);
  return run_header;
}

const DataModel::CameraEvent* ZFITSACTLDataSource::
borrow_next_event(uint64_t& seq_index_out)
{
  if(config_.max_seq_index() and seq_index_>=config_.max_seq_index())
    return nullptr;
  while(isource_ < opener_->num_sources())
  {
    uint64_t unused_index = 0;
    if(const data_type* next = source_->borrow_next_event(unused_index))
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

void ZFITSACTLDataSource::
release_borrowed_event(const DataModel::CameraEvent* event)
{
  if(source_)
    source_->release_borrowed_event(event);
  else
    delete event;
}

DataModel::CameraEvent* ZFITSACTLDataSource::get_next(uint64_t& seq_index_out,
  google::protobuf::Arena** arena)
{
  if(config_.max_seq_index() and seq_index_>=config_.max_seq_index())
    return nullptr;
  else
    return calin::io::data_source::BasicChainedRandomAccessDataSource<
      ACTLRandomAccessDataSourceWithRunHeader>::get_next(seq_index_out, arena);
}

uint64_t ZFITSACTLDataSource::size()
{
  uint64_t max_seq_index = calin::io::data_source::BasicChainedRandomAccessDataSource<
    ACTLRandomAccessDataSourceWithRunHeader>::size();
  if(config_.max_seq_index() and max_seq_index>config_.max_seq_index())
    max_seq_index = config_.max_seq_index();
  return max_seq_index;
}

void ZFITSACTLDataSource::set_next_index(uint64_t next_index)
{
  if(config_.max_seq_index() and next_index>config_.max_seq_index())
    next_index = config_.max_seq_index();
  calin::io::data_source::BasicChainedRandomAccessDataSource<
    ACTLRandomAccessDataSourceWithRunHeader>::set_next_index(next_index);
}

ZFITSConstACTLDataSourceBorrowAdapter::
ZFITSConstACTLDataSourceBorrowAdapter(ZFITSACTLDataSource* src):
  calin::io::data_source::DataSource<const DataModel::CameraEvent>(),
  src_(src)
{
  // nothing to see here
}

ZFITSConstACTLDataSourceBorrowAdapter::~ZFITSConstACTLDataSourceBorrowAdapter()
{
  // nothing to see here
}

const DataModel::CameraEvent* ZFITSConstACTLDataSourceBorrowAdapter::
get_next(uint64_t& seq_index_out, google::protobuf::Arena** arena)
{
  assert(arena==nullptr or *arena==nullptr);
  return src_->borrow_next_event(seq_index_out);
}

ZFITSConstACTLDataSourceReleaseAdapter::
ZFITSConstACTLDataSourceReleaseAdapter(ZFITSACTLDataSource* src):
  calin::io::data_source::DataSink<const DataModel::CameraEvent>(),
  src_(src)
{
  // nothing to see here
}

ZFITSConstACTLDataSourceReleaseAdapter::~ZFITSConstACTLDataSourceReleaseAdapter()
{
  // nothing to see here
}

bool ZFITSConstACTLDataSourceReleaseAdapter::
put_next(const DataModel::CameraEvent* data, uint64_t seq_index,
  google::protobuf::Arena* arena, bool adopt_data)
{
  assert(adopt_data);
  assert(arena==nullptr);
  src_->release_borrowed_event(data);
  return true;
}
