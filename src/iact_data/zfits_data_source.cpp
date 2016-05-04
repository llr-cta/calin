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
  auto* sample_event = actl_zfits_->get_next();
  run_config_ = decoder_->decode_run_config(actl_zfits_->get_run_header(),
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
ZFITSDataSource::get_next()
{
  auto index = actl_zfits_->next_index();
  std::unique_ptr<DataModel::CameraEvent> cta_event {
    actl_zfits_->get_next() };
  if(!cta_event)return nullptr;
  TelescopeEvent* event = decoder_->decode(cta_event.get());
  if(event)event->set_source_event_index(index-1);
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

uint64_t ZFITSDataSource::next_index()
{
  return actl_zfits_->next_index();
}

void ZFITSDataSource::set_next_index(uint64_t next_index)
{
  return actl_zfits_->set_next_index(next_index);
}
