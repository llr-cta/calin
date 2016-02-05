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
    CTACameraEventDecoder* decoder, bool adopt_decoder):
  calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSource(),
  filename_(expand_filename(filename)), decoder_(decoder),
  adopt_decoder_(adopt_decoder)
{
  if(!is_file(filename_))
    throw std::runtime_error(std::string("No such file: ")+filename_);
  if(!is_readable(filename_))
    throw std::runtime_error(std::string("File not readable: ")+filename_);
  zfits_ = new ACTL::IO::ProtobufIFits(filename_.c_str());
}

ZFITSSingleFileDataSource::~ZFITSSingleFileDataSource()
{
  delete zfits_;
  if(adopt_decoder_)delete decoder_;
}

TelescopeEvent* ZFITSSingleFileDataSource::get_next()
{
  if(zfits_ == nullptr)
    throw std::runtime_error(std::string("File not open: ")+filename_);
  if(next_event_index_ >= zfits_->getNumMessagesInTable())return nullptr;

  std::unique_ptr<DataModel::CameraEvent> cta_event {
    zfits_->readTypedMessage<DataModel::CameraEvent>(++next_event_index_) };

  if(!cta_event)throw runtime_error("ZFits reader returned NULL");

  return decoder_->decode(cta_event.get());
}

uint64_t ZFITSSingleFileDataSource::size()
{
  return zfits_->getNumMessagesInTable();
}

void ZFITSSingleFileDataSource::set_next_index(uint64_t next_index)
{
  next_event_index_ = next_index;
}

ZFitsDataSourceOpener::ZFitsDataSourceOpener(std::string filename,
  CTACameraEventDecoder* decoder, bool adopt_decoder,
  const ZFITSDataSource::config_type& config):
  DataSourceOpener<calin::iact_data::
    telescope_data_source::TelescopeRandomAccessDataSource>(),
  decoder_(decoder), adopt_decoder_(adopt_decoder), config_(config)
{
  if(is_file(filename))
  {
    filenames_.emplace_back(filename);
    return;
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
    return;
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

calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSource*
ZFitsDataSourceOpener::open(unsigned isource)
{
  if(isource >= filenames_.size())return nullptr;
  if(config_.log_on_file_open())
    LOG(INFO) << "Opening file: " << filenames_[isource];
  return new ZFITSSingleFileDataSource(filenames_[isource], decoder_, false);
}

//ZFITSDataSource::config_helper ZFITSDataSource::default_config_;

ZFITSDataSource::ZFITSDataSource(const std::string& filename,
  CTACameraEventDecoder* decoder, bool adopt_decoder,
  const config_type& config):
  BasicChaninedRandomAccessDataSource<calin::iact_data::
      telescope_data_source::TelescopeRandomAccessDataSource>(
    new ZFitsDataSourceOpener(filename,decoder,adopt_decoder, config), true),
  config_(config)
{
  // nothing to see here
}

ZFITSDataSource::~ZFITSDataSource()
{
  // nothing to see here
}
