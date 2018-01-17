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
  filename_(expand_filename(filename))
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
  if(next_event_index_ >= zfits_->getNumMessagesInTable())return nullptr;

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
  return zfits_->getNumMessagesInTable();
}

void ZFITSSingleFileACTLDataSource::set_next_index(uint64_t next_index)
{
  if(zfits_ == nullptr)next_event_index_ = 0;
  else next_event_index_ =
    std::min(next_index, uint64_t(zfits_->getNumMessagesInTable()));
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
      if(not is_file(filename+".1"+extension) and not is_file(filename+".001"+extension))
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

        filename_i = filename + ".";
        if(i<10) filename_i += "00";
        else if(i<100) filename_i += "0";
        else break; // no point in rechecking same file as above
        filename_i += std::to_string(i)+extension;

        if(not is_file(filename_i))break;
        filenames_.emplace_back(filename_i);
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
