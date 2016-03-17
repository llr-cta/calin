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

ZFITSSingleFileDataSource::
ZFITSSingleFileDataSource(const std::string& filename,
    CTACameraEventDecoder* decoder, bool adopt_decoder,
    const config_type& config):
  calin::iact_data::telescope_data_source::
    TelescopeRandomAccessDataSourceWithRunConfig(),
  filename_(expand_filename(filename)), decoder_(decoder),
  adopt_decoder_(adopt_decoder)
{
  if(!is_file(filename_))
    throw std::runtime_error(std::string("No such file: ")+filename_);
  if(!is_readable(filename_))
    throw std::runtime_error(std::string("File not readable: ")+filename_);

  const DataModel::CameraRunHeader* cta_run_header = nullptr;
  if(!config.dont_read_run_header())
  {
    try
    {
      std::unique_ptr<ACTL::IO::ProtobufIFits> rh_zfits {
        new ACTL::IO::ProtobufIFits(filename_.c_str(),
          "RunHeader", DataModel::CameraRunHeader::descriptor()) };
      cta_run_header = zfits_->readTypedMessage<DataModel::CameraRunHeader>(1);
    }
    catch(...)
    {
      if(!config.ignore_run_header_errors())
        LOG(WARNING) << "ZFITSSingleFileDataSource: Could not read RunHeader from "
          << filename_;
    }
  }

  zfits_ = new ACTL::IO::ProtobufIFits(filename_.c_str());

  const DataModel::CameraEvent* cta_event =
    zfits_->readTypedMessage<DataModel::CameraEvent>(1);
  run_config_ = decoder_->decode_run_config(cta_run_header, cta_event);
}

ZFITSSingleFileDataSource::~ZFITSSingleFileDataSource()
{
  delete zfits_;
  delete run_config_;
  if(adopt_decoder_)delete decoder_;
}

ZFITSSingleFileDataSource::config_type
ZFITSSingleFileDataSource::default_config()
{
  return ZFITSDataSource::default_config();
}

TelescopeEvent* ZFITSSingleFileDataSource::get_next()
{
  if(zfits_ == nullptr)
    throw std::runtime_error(std::string("File not open: ")+filename_);
  if(next_event_index_ >= zfits_->getNumMessagesInTable())return nullptr;

  std::unique_ptr<DataModel::CameraEvent> cta_event {
    zfits_->readTypedMessage<DataModel::CameraEvent>(++next_event_index_) };

  if(!cta_event)throw runtime_error("ZFits reader returned NULL");

  TelescopeEvent* event = decoder_->decode(cta_event.get());
  if(event)event->set_source_event_index(next_event_index_-1);
  return event;
}

uint64_t ZFITSSingleFileDataSource::size()
{
  return zfits_->getNumMessagesInTable();
}

void ZFITSSingleFileDataSource::set_next_index(uint64_t next_index)
{
  next_event_index_ = next_index;
}

TelescopeRunConfiguration* ZFITSSingleFileDataSource::get_run_configuration()
{
  if(!run_config_)return nullptr;
  auto* run_config = new TelescopeRunConfiguration();
  run_config->CopyFrom(*run_config_);
  return run_config;
}

ZFitsDataSourceOpener::ZFitsDataSourceOpener(std::string filename,
  bool exact_filename_only,
  CTACameraEventDecoder* decoder, bool adopt_decoder,
  const ZFITSDataSource::config_type& config):
  DataSourceOpener<calin::iact_data::
    telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig>(),
  decoder_(decoder), adopt_decoder_(adopt_decoder), config_(config)
{
  if(is_file(filename))
  {
    filenames_.emplace_back(filename);
  }
  else if(exact_filename_only)
  {
    throw(std::runtime_error("File not found: "+ filename));
  }

  std::string extension = config_.extension();
  if(filename.size() > extension.size() and
    filename.rfind(extension) == filename.size()-extension.size())
  {
    filename = filename.substr(0, filename.size()-extension.size());
  }
  else if(is_file(filename+extension))
  {
    filenames_.emplace_back(filename+extension);
  }

  for(unsigned i=1; true; ++i)
  {
    std::string filename_i { filename+"."+std::to_string(i)+extension };
    if(not is_file(filename_i))break;
    filenames_.emplace_back(filename_i);
  }

  if(filenames_.empty())
    throw(std::runtime_error("File not found: "+ filename+extension
      + " and " + filename+".1"+extension));
}

ZFitsDataSourceOpener::~ZFitsDataSourceOpener()
{
  if(adopt_decoder_)delete decoder_;
}

unsigned ZFitsDataSourceOpener::num_sources()
{
  return filenames_.size();
}

calin::iact_data::telescope_data_source::
TelescopeRandomAccessDataSourceWithRunConfig*
ZFitsDataSourceOpener::open(unsigned isource)
{
  if(isource >= filenames_.size())return nullptr;
  if(config_.log_on_file_open())
    LOG(INFO) << "Opening file: " << filenames_[isource];
  auto config = config_;
  if(isource != 0)config.set_dont_read_run_header(true);
  return new ZFITSSingleFileDataSource(filenames_[isource], decoder_,
    false, config);
}

//ZFITSDataSource::config_helper ZFITSDataSource::default_config_;

ZFITSDataSource::ZFITSDataSource(const std::string& filename,
  bool exact_filename_only,
  CTACameraEventDecoder* decoder, bool adopt_decoder,
  const config_type& config):
  BasicChainedRandomAccessDataSource<calin::iact_data::
      telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig>(
    new ZFitsDataSourceOpener(filename, exact_filename_only, decoder,
      adopt_decoder, config), true),
  config_(config), run_config_(source_->get_run_configuration())
{
  // nothing to see here
}

ZFITSDataSource::~ZFITSDataSource()
{
  delete run_config_;
}

calin::ix::iact_data::telescope_event::TelescopeEvent*
ZFITSDataSource::get_next()
{
  calin::ix::iact_data::telescope_event::TelescopeEvent* event =
    BasicChainedRandomAccessDataSource<calin::iact_data::
      telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig>::get_next();
  if(event and isource_)
    event->set_source_event_index(event->source_event_index() +
      chained_file_index_[isource_-1]);
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
