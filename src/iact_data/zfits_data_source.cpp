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
    CTACameraEventDecoder* decoder, bool adopt_decoder,
    const config_type& config):
  calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSource(),
  filename_(expand_filename(filename)), config_(config), decoder_(decoder),
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
  if(config_.next_event_index() >= zfits_->getNumMessagesInTable())
    return nullptr;

  std::unique_ptr<DataModel::CameraEvent> cta_event {
    zfits_->readTypedMessage<DataModel::CameraEvent>(
      config_.next_event_index()+1) };

  config_.set_next_event_index(config_.next_event_index()+1);

  if(!cta_event)throw runtime_error("ZFits reader returned NULL");
  if(auto num_events_left = config_.num_events_max())
  {
    if(num_events_left == 1)
      config_.set_next_event_index(zfits_->getNumMessagesInTable());
    config_.set_num_events_max(num_events_left-1);
  }

  return decoder_->decode(cta_event.get());
}

uint64_t ZFITSSingleFileDataSource::size()
{
  return zfits_->getNumMessagesInTable();
}

void ZFITSSingleFileDataSource::set_next_index(uint64_t next_index)
{
  config_.set_next_event_index(next_index);
}

ZFitsDataSourceOpener::ZFitsDataSourceOpener(std::string filename,
  CTACameraEventDecoder* decoder, bool adopt_decoder,
  const std::string& extension):
  DataSourceOpener<calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSource>(),
  decoder_(decoder), adopt_decoder_(adopt_decoder)
{
  if(is_file(filename))
  {
    filenames_.push_back(filename);
    return;
  }

  if(filename.size() > extension.size()
      and filename.rfind(extension) == filename.size()-extension.size())
  {
    filename = filename.substr(0, filename.size()-extension.size());
  }
  else if(is_file(filename+extension))
  {
    filenames_.push_back(filename+extension);
    return;
  }

  for(unsigned i=1; is_file(filename+"."+std::to_string(i)+extension); i++)
    filenames_.push_back(filename+"."+std::to_string(i)+extension);
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
  return new ZFITSSingleFileDataSource(filenames_[isource], decoder_, false);
}

ZFITSDataSource::ZFITSDataSource(const std::string& filename,
  CTACameraEventDecoder* decoder, bool adopt_decoder,
  const std::string& extension):
  BasicChaninedRandomAccessDataSource<calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSource>(
    new ZFitsDataSourceOpener(filename,decoder,adopt_decoder,extension), true)
{
  // nothing to see here
}

ZFITSDataSource::~ZFITSDataSource()
{
  // nothing to see here
}
