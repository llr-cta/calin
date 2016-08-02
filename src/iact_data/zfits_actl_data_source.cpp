/*

   calin/iact_data/zfits_actl_data_source.cpp -- Stephen Fegan -- 2016-05-04

   A supplier of single telescope ACTL data from ZFits data files

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
#include <iact_data/zfits_actl_data_source.hpp>
#include <ProtobufIFits.h>
#include <L0.pb.h>

using namespace calin::iact_data::zfits_actl_data_source;
using namespace calin::io::log;
using calin::util::file::is_file;
using calin::util::file::is_readable;
using calin::util::file::expand_filename;

ACTLRandomAccessDataSourceWithRunHeader::
~ACTLRandomAccessDataSourceWithRunHeader()
{
  // nothing to see here
}

ZFITSSingleFileACTLDataSource::
ZFITSSingleFileACTLDataSource(const std::string& filename,
    const config_type& config):
  ACTLRandomAccessDataSourceWithRunHeader(),
  filename_(expand_filename(filename))
{
  if(!is_file(filename_))
    throw std::runtime_error(std::string("No such file: ")+filename_);
  if(!is_readable(filename_))
    throw std::runtime_error(std::string("File not readable: ")+filename_);

  if(!config.dont_read_run_header())
  {
    try
    {
      std::unique_ptr<ACTL::IO::ProtobufIFits> rh_zfits {
        new ACTL::IO::ProtobufIFits(filename_.c_str(),
          "RunHeader", DataModel::CameraRunHeader::descriptor()) };
      run_header_ = zfits_->readTypedMessage<DataModel::CameraRunHeader>(1);
    }
    catch(...)
    {
      if(!config.ignore_run_header_errors())
        LOG(WARNING)
          << "ZFITSSingleFileACTLDataSource: Could not read RunHeader from "
          << filename_;
    }
  }

  zfits_ = new ACTL::IO::ProtobufIFits(filename_.c_str(), "Events");
}

ZFITSSingleFileACTLDataSource::~ZFITSSingleFileACTLDataSource()
{
  delete zfits_;
  delete run_header_;
}

ZFITSSingleFileACTLDataSource::config_type
ZFITSSingleFileACTLDataSource::default_config()
{
  return ZFITSACTLDataSource::default_config();
}

DataModel::CameraEvent* ZFITSSingleFileACTLDataSource::
get_next(uint64_t& seq_index_out, google::protobuf::Arena** arena)
{
  if(zfits_ == nullptr)
    throw std::runtime_error(std::string("File not open: ")+filename_);
  if(next_event_index_ >= zfits_->getNumMessagesInTable())return nullptr;

  seq_index_out = next_event_index_;
  DataModel::CameraEvent* event {
    zfits_->readTypedMessage<DataModel::CameraEvent>(++next_event_index_) };
  if(!event)throw runtime_error("ZFits reader returned NULL");

  return event;
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

      for(unsigned i=istart+1; config_.max_file_fragments()==0 or
        filenames_.size()<config_.max_file_fragments() ; ++i)
      {
        std::string filename_i { filename+"."+std::to_string(i)+extension };
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
ACTLRandomAccessDataSourceWithRunHeader*
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
