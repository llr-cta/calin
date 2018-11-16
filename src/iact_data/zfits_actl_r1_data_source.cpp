/*

   calin/iact_data/zfits_actl_r1_data_source.cpp -- Stephen Fegan -- 2018-09-20

   A supplier of single telescope ACTL R1 CameraEvent data from ZFits data files

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
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

/*

                      RRRRRRRRRRRRRRRRR          1111111
                      R::::::::::::::::R        1::::::1
                      R::::::RRRRRR:::::R      1:::::::1
                      RR:::::R     R:::::R     111:::::1
                        R::::R     R:::::R        1::::1
                        R::::R     R:::::R        1::::1
                        R::::RRRRRR:::::R         1::::1
                        R:::::::::::::RR          1::::l
                        R::::RRRRRR:::::R         1::::l
                        R::::R     R:::::R        1::::l
                        R::::R     R:::::R        1::::l
                        R::::R     R:::::R        1::::l
                      RR:::::R     R:::::R     111::::::111
                      R::::::R     R:::::R     1::::::::::1
                      R::::::R     R:::::R     1::::::::::1
                      RRRRRRRR     RRRRRRR     111111111111

*/

#include <stdexcept>
#include <memory>
#include <cctype>

#include <util/log.hpp>
#include <provenance/chronicle.hpp>
#include <util/file.hpp>
#include <iact_data/zfits_actl_data_source.hpp>
#include <ProtobufIFits.h>
#include <R1.pb.h>

using namespace calin::iact_data::zfits_actl_data_source;
using namespace calin::util::log;
using calin::util::file::is_file;
using calin::util::file::is_readable;
using calin::util::file::expand_filename;

bool calin::iact_data::zfits_actl_data_source::
is_zfits_r1(std::string filename, std::string events_table_name)
{
  calin::util::file::expand_filename_in_place(filename);
  if(!is_file(filename))return false;

  if(events_table_name == "")events_table_name =
    ZFITSSingleFileACTL_R1_CameraEventDataSource::
      default_config().events_table_name();

  try {
    IFits ifits(filename, events_table_name, /* force= */ true);
    const std::string& message_name = ifits.GetStr("PBFHEAD");
    if (message_name == "R1.CameraEvent") {
      return true;
    }
    // fallthrough to return false
  } catch(...) {
    // fallthrough to return false
  }
  return false;
}

ACTL_R1_CameraEventDataSourceWithRunHeader::
~ACTL_R1_CameraEventDataSourceWithRunHeader()
{
  // nothing to see here
}

ACTL_R1_CameraEventRandomAccessDataSourceWithRunHeader::
~ACTL_R1_CameraEventRandomAccessDataSourceWithRunHeader()
{
  // nothing to see here
}

namespace {
  static std::string default_R1_run_header_table_name("CameraConfig");
  static std::string default_R1_events_table_name("Events");
} // anonymous namespace


ZFITSSingleFileACTL_R1_CameraEventDataSource::
ZFITSSingleFileACTL_R1_CameraEventDataSource(const std::string& filename, config_type config):
  ACTL_R1_CameraEventRandomAccessDataSourceWithRunHeader(),
  filename_(expand_filename(filename)), config_(config)
{
  if(config.run_header_table_name().empty())
    config.set_run_header_table_name(default_R1_run_header_table_name);
  if(config.events_table_name().empty())
    config.set_events_table_name(default_R1_events_table_name);

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
        R1::CameraConfiguration::descriptor());
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
        R1::CameraConfiguration* run_header =
          rh_zfits.readTypedMessage<R1::CameraConfiguration>(1);
        run_header_ = new R1::CameraConfiguration(*run_header);
        rh_zfits.recycleMessage(run_header);
        //LOG(INFO) << run_header_->DebugString();
      }
    }
    catch(...)
    {
      if(!config.ignore_run_header_errors())
        LOG(WARNING)
          << "ZFITSSingleFileACTL_R1_CameraEventDataSource: Could not read run header from "
          << filename_;
    }
  }

  zfits_ = new ACTL::IO::ProtobufIFits(filename_.c_str(),
    config.events_table_name(), R1::CameraEvent::descriptor());
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

ZFITSSingleFileACTL_R1_CameraEventDataSource::~ZFITSSingleFileACTL_R1_CameraEventDataSource()
{
  delete zfits_;
  delete run_header_;
}

ZFITSSingleFileACTL_R1_CameraEventDataSource::config_type
ZFITSSingleFileACTL_R1_CameraEventDataSource::default_config()
{
  config_type config = config_type::default_instance();
  config.set_data_model(calin::ix::iact_data::zfits_data_source::ACTL_DATA_MODEL_R1);
  config.set_extension(".fits.fz");
  config.set_run_header_table_name(default_R1_run_header_table_name);
  config.set_events_table_name(default_R1_events_table_name);
  config.set_file_fragment_stride(1);
  return config;
}

const R1::CameraEvent*
ZFITSSingleFileACTL_R1_CameraEventDataSource::borrow_next_event(uint64_t& seq_index_out)
{
  if(zfits_ == nullptr)
    throw std::runtime_error(std::string("File not open: ")+filename_);

  uint64_t max_seq_index = zfits_->getNumMessagesInTable();
  if(config_.max_seq_index())
    max_seq_index = std::min(max_seq_index, config_.max_seq_index());
  if(next_event_index_ >= max_seq_index)return nullptr;

  seq_index_out = next_event_index_;
  const R1::CameraEvent* event {
    zfits_->borrowTypedMessage<R1::CameraEvent>(++next_event_index_) };
  if(!event)throw std::runtime_error("ZFits reader returned NULL");

  return event;
}

void ZFITSSingleFileACTL_R1_CameraEventDataSource::
release_borrowed_event(const R1::CameraEvent* event)
{
  zfits_->returnBorrowedMessage(event);
}

R1::CameraEvent* ZFITSSingleFileACTL_R1_CameraEventDataSource::
get_next(uint64_t& seq_index_out, google::protobuf::Arena** arena)
{
  if(arena)*arena = nullptr;
  const R1::CameraEvent* event = borrow_next_event(seq_index_out);
  if(event == nullptr)return nullptr;

  R1::CameraEvent* event_copy = nullptr;
#if 0
  if(arena) {
    if(!*arena)*arena = new google::protobuf::Arena;
    event_copy =
      google::protobuf::Arena::CreateMessage<R1::CameraEvent>(*arena);
  }
  else event_copy = new R1::CameraEvent;
#else
  if(arena && *arena)
    throw std::runtime_error("ZFITSSingleFileACTL_R1_CameraEventDataSource::get_next: "
      " pre-allocated arena not supported.");
  event_copy = new R1::CameraEvent;
#endif
  event_copy->CopyFrom(*event);
  release_borrowed_event(event);

  return event_copy;
}

uint64_t ZFITSSingleFileACTL_R1_CameraEventDataSource::size()
{
  uint64_t max_seq_index = zfits_->getNumMessagesInTable();
  if(config_.max_seq_index())
    max_seq_index = std::min(max_seq_index, config_.max_seq_index());
  return max_seq_index;
}

void ZFITSSingleFileACTL_R1_CameraEventDataSource::set_next_index(uint64_t next_index)
{
  if(zfits_ == nullptr)next_event_index_ = 0;
  else {
    uint64_t max_seq_index = zfits_->getNumMessagesInTable();
    if(config_.max_seq_index())
      max_seq_index = std::min(max_seq_index, config_.max_seq_index());
    next_event_index_ = std::min(next_index, max_seq_index);
  }
}

R1::CameraConfiguration* ZFITSSingleFileACTL_R1_CameraEventDataSource::get_run_header()
{
  if(!run_header_)return nullptr;
  auto* run_header = new R1::CameraConfiguration();
  run_header->CopyFrom(*run_header_);
  return run_header;
}

ZFITSACTL_R1_CameraEventDataSourceOpener::ZFITSACTL_R1_CameraEventDataSourceOpener(std::string filename,
  const ZFITSACTL_R1_CameraEventDataSource::config_type& config):
  calin::io::data_source::DataSourceOpener<
    ACTL_R1_CameraEventRandomAccessDataSourceWithRunHeader>(),
  config_(config)
{
  const unsigned istride = std::max(1U,config.file_fragment_stride());
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

ZFITSACTL_R1_CameraEventDataSourceOpener::~ZFITSACTL_R1_CameraEventDataSourceOpener()
{
  // nothing to see here
}

unsigned ZFITSACTL_R1_CameraEventDataSourceOpener::num_sources()
{
  return filenames_.size();
}

std::string ZFITSACTL_R1_CameraEventDataSourceOpener::source_name(unsigned isource)
{
  if(isource >= filenames_.size())return {};
  return filenames_[isource];
}

calin::iact_data::zfits_actl_data_source::
ZFITSSingleFileACTL_R1_CameraEventDataSource*
ZFITSACTL_R1_CameraEventDataSourceOpener::open(unsigned isource)
{
  if(isource >= filenames_.size())return nullptr;
  if(config_.log_on_file_open())
    LOG(INFO) << "Opening file: " << filenames_[isource];
  auto config = config_;
  if(has_opened_file_)config.set_dont_read_run_header(true);
  config.set_max_seq_index(0);
  has_opened_file_ = true;
  return new ZFITSSingleFileACTL_R1_CameraEventDataSource(filenames_[isource], config);
}

ZFITSACTL_R1_CameraEventDataSource::ZFITSACTL_R1_CameraEventDataSource(const std::string& filename,
  const config_type& config):
  calin::io::data_source::BasicChainedRandomAccessDataSource<
    ACTL_R1_CameraEventRandomAccessDataSourceWithRunHeader>(
    new ZFITSACTL_R1_CameraEventDataSourceOpener(filename, config), true),
  config_(config), run_header_(source_->get_run_header())
{
  // nothing to see here
}

ZFITSACTL_R1_CameraEventDataSource::~ZFITSACTL_R1_CameraEventDataSource()
{
  delete run_header_;
}

R1::CameraConfiguration* ZFITSACTL_R1_CameraEventDataSource::get_run_header()
{
  if(!run_header_)return nullptr;
  auto* run_header = new R1::CameraConfiguration();
  run_header->CopyFrom(*run_header_);
  return run_header;
}

const R1::CameraEvent* ZFITSACTL_R1_CameraEventDataSource::
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

void ZFITSACTL_R1_CameraEventDataSource::
release_borrowed_event(const R1::CameraEvent* event)
{
  if(source_)
    source_->release_borrowed_event(event);
  else
    delete event;
}

R1::CameraEvent* ZFITSACTL_R1_CameraEventDataSource::get_next(uint64_t& seq_index_out,
  google::protobuf::Arena** arena)
{
  if(config_.max_seq_index() and seq_index_>=config_.max_seq_index())
    return nullptr;
  else
    return calin::io::data_source::BasicChainedRandomAccessDataSource<
      ACTL_R1_CameraEventRandomAccessDataSourceWithRunHeader>::get_next(seq_index_out, arena);
}

uint64_t ZFITSACTL_R1_CameraEventDataSource::size()
{
  uint64_t max_seq_index = calin::io::data_source::BasicChainedRandomAccessDataSource<
    ACTL_R1_CameraEventRandomAccessDataSourceWithRunHeader>::size();
  if(config_.max_seq_index() and max_seq_index>config_.max_seq_index())
    max_seq_index = config_.max_seq_index();
  return max_seq_index;
}

void ZFITSACTL_R1_CameraEventDataSource::set_next_index(uint64_t next_index)
{
  if(config_.max_seq_index() and next_index>config_.max_seq_index())
    next_index = config_.max_seq_index();
  calin::io::data_source::BasicChainedRandomAccessDataSource<
    ACTL_R1_CameraEventRandomAccessDataSourceWithRunHeader>::set_next_index(next_index);
}

ZFITSConstACTL_R1_CameraEventDataSourceBorrowAdapter::
ZFITSConstACTL_R1_CameraEventDataSourceBorrowAdapter(ZFITSACTL_R1_CameraEventDataSource* src):
  calin::io::data_source::DataSource<const R1::CameraEvent>(),
  src_(src)
{
  // nothing to see here
}

ZFITSConstACTL_R1_CameraEventDataSourceBorrowAdapter::~ZFITSConstACTL_R1_CameraEventDataSourceBorrowAdapter()
{
  // nothing to see here
}

const R1::CameraEvent* ZFITSConstACTL_R1_CameraEventDataSourceBorrowAdapter::
get_next(uint64_t& seq_index_out, google::protobuf::Arena** arena)
{
  assert(arena==nullptr or *arena==nullptr);
  return src_->borrow_next_event(seq_index_out);
}

ZFITSConstACTL_R1_CameraEventDataSourceReleaseAdapter::
ZFITSConstACTL_R1_CameraEventDataSourceReleaseAdapter(ZFITSACTL_R1_CameraEventDataSource* src):
  calin::io::data_source::DataSink<const R1::CameraEvent>(),
  src_(src)
{
  // nothing to see here
}

ZFITSConstACTL_R1_CameraEventDataSourceReleaseAdapter::~ZFITSConstACTL_R1_CameraEventDataSourceReleaseAdapter()
{
  // nothing to see here
}

bool ZFITSConstACTL_R1_CameraEventDataSourceReleaseAdapter::
put_next(const R1::CameraEvent* data, uint64_t seq_index,
  google::protobuf::Arena* arena, bool adopt_data)
{
  assert(adopt_data);
  assert(arena==nullptr);
  src_->release_borrowed_event(data);
  return true;
}


ZMQACTL_R1_CameraEventDataSource::
ZMQACTL_R1_CameraEventDataSource(
    const std::string& endpoint, void* zmq_ctx, int buffer_size):
  ACTL_R1_CameraEventDataSourceWithRunHeader(),
  zmq_ctx_((zmq_ctx == nullptr)?calin::io::zmq_inproc::new_zmq_ctx():zmq_ctx),
  my_zmq_ctx_((zmq_ctx == nullptr)?zmq_ctx_:nullptr),
  zmq_source_(new calin::io::data_source::ZMQProtobufDataSource<DataModel::CTAMessage>(
    zmq_ctx_, endpoint, buffer_size))
{
  // nothing to see here
}

ZMQACTL_R1_CameraEventDataSource::~ZMQACTL_R1_CameraEventDataSource()
{
  delete zmq_source_;
  if(my_zmq_ctx_)calin::io::zmq_inproc::destroy_zmq_ctx(my_zmq_ctx_);
}

R1::CameraConfiguration* ZMQACTL_R1_CameraEventDataSource::get_run_header()
{
  if(run_header_)return new R1::CameraConfiguration(*run_header_);
  if(!saved_event_)saved_event_ = get_next(saved_seq_index_, nullptr);
  if(run_header_)return new R1::CameraConfiguration(*run_header_);
  return nullptr;
}

R1::CameraEvent* ZMQACTL_R1_CameraEventDataSource::
get_next(uint64_t& seq_index_out, google::protobuf::Arena** arena)
{
  if(saved_event_) {
    R1::CameraEvent* event = saved_event_;
    seq_index_out = saved_seq_index_;
    saved_event_ = nullptr;
    return event;
  }
  if(eos_received_)return nullptr;

  DataModel::CTAMessage* message = nullptr;
  while(1) {
    message = zmq_source_->get_next(seq_index_out);
    if(message == nullptr)return nullptr;
    if(message->payload_type_size() == 0) {
      delete message;
      return nullptr;
    }
    switch(message->payload_type(0)) {
      case DataModel::EMPTY_MESSAGE:
        delete message;
        return nullptr;
      case DataModel::END_OF_STREAM:
        delete message;
        eos_received_ = true;
        return nullptr;
      case DataModel::R1_CAMERA_EVENT:
      {
        R1::CameraEvent* event = new R1::CameraEvent();
        if(event->ParseFromString(message->payload_data(0))) {
          delete message;
          return event;
        } else {
          delete message;
          throw std::runtime_error("ZMQACTL_R1_CameraEventDataSource::get_next: "
            "could not parse event payload");
        }
      }
      case DataModel::R1_CAMERA_CONFIG:
      {
        if(run_header_ == nullptr)run_header_ = new R1::CameraConfiguration();
        if(run_header_->ParseFromString(message->payload_data(0))) {
          delete message;
          continue;
        } else {
          delete message;
          throw std::runtime_error("ZMQACTL_R1_CameraEventDataSource::get_next: "
            "could not parse caamera configuration payload");
        }
      }
      default:
        delete message;
        throw std::runtime_error("ZMQACTL_R1_CameraEventDataSource::get_next: "
          "unknown message payload type: "
          + DataModel::MessageType_Name(message->payload_type(0)));
    }
  }
}
