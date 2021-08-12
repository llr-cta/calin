/*

   calin/provenance/chronicle.cpp -- Stephen Fegan -- 2017-10-30

   Provenance system for recording noteworthy run-time events (io files etc)

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

   Based on original, copyright 2006, Stephen Fegan, see notice below

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <mutex>
#include <memory>
#include <sys/stat.h>

#include <provenance/chronicle.hpp>
#include <util/timestamp.hpp>

using namespace calin::provenance::chronicle;

namespace {

std::mutex chronicle_mutex;

std::unique_ptr<calin::ix::provenance::chronicle::Chronicle> singleton_chronicle {
  new calin::ix::provenance::chronicle::Chronicle };

template<typename RepeatedFieldType> void prune_repeated_field(RepeatedFieldType* rf)
{
  int from = 0;
  int to = 0;
  while(from != rf->size()) {
    while(to != rf->size() and not rf->Get(to).has_close_timestamp())++to;
    from = std::max(to, from);
    while(from != rf->size() and rf->Get(from).has_close_timestamp())++from;
    if(from != rf->size()) {
      rf->SwapElements(to, from);
      ++to;
      ++from;
    }
  }
  if(to == 0) {
    rf->Clear();
  } else {
    while(rf->size() > to)rf->RemoveLast();
  }
}

} // anonymous namespace

const calin::ix::provenance::chronicle::Chronicle*
calin::provenance::chronicle::the_chronicle()
{
  return singleton_chronicle.get();
}

void calin::provenance::chronicle::prune_the_chronicle()
{
  std::lock_guard<std::mutex> lock { chronicle_mutex };
  prune_repeated_field(singleton_chronicle->mutable_file_io_record());
  prune_repeated_field(singleton_chronicle->mutable_network_io_record());
  prune_repeated_field(singleton_chronicle->mutable_rng_record());
  prune_repeated_field(singleton_chronicle->mutable_processing_record());
}

calin::ix::provenance::chronicle::Chronicle*
calin::provenance::chronicle::copy_the_chronicle(calin::ix::provenance::chronicle::Chronicle* x)
{
  if(x == nullptr) x = new calin::ix::provenance::chronicle::Chronicle;
  std::lock_guard<std::mutex> lock { chronicle_mutex };
  x->CopyFrom(*the_chronicle());
  return x;
}

calin::ix::provenance::chronicle::FileIORecord* calin::provenance::chronicle::
register_file_open(const std::string& file_name,
  calin::ix::provenance::chronicle::AccessType access,
  const std::string& opened_by, const std::string& comment)
{
  calin::util::timestamp::Timestamp ts = calin::util::timestamp::Timestamp::now();
  calin::ix::provenance::chronicle::FileIORecord* record = nullptr;
  {
    std::lock_guard<std::mutex> lock { chronicle_mutex };
    record = singleton_chronicle->add_file_io_record();
  }
  ts.as_proto(record->mutable_open_timestamp());
  record->set_access(access);
  record->set_file_name(file_name);
  record->set_comment(comment);
  struct stat stat_buffer;
  if(::stat(file_name.c_str(), &stat_buffer) >= 0) {
    record->set_file_size(stat_buffer.st_size);
    calin::util::timestamp::Timestamp(stat_buffer.st_mtime, 1000000000LL).as_proto(record->mutable_file_mtime());
  }
  record->set_opened_by(opened_by);
  return record;
}

void calin::provenance::chronicle::
register_file_close(calin::ix::provenance::chronicle::FileIORecord* record)
{
  calin::util::timestamp::Timestamp ts = calin::util::timestamp::Timestamp::now();
  ts.as_proto(record->mutable_close_timestamp());
}

calin::ix::provenance::chronicle::NetworkIORecord* calin::provenance::chronicle::
register_network_open(const std::string& endpoint,
  const std::string& opened_by, const std::string& comment)
{
  calin::util::timestamp::Timestamp ts = calin::util::timestamp::Timestamp::now();
  calin::ix::provenance::chronicle::NetworkIORecord* record = nullptr;
  {
    std::lock_guard<std::mutex> lock { chronicle_mutex };
    record = singleton_chronicle->add_network_io_record();
  }
  ts.as_proto(record->mutable_open_timestamp());
  record->set_endpoint_name(endpoint);
  record->set_comment(comment);
  record->set_nbytes_received(-1);
  record->set_nbytes_sent(-1);
  record->set_opened_by(opened_by);
  return record;
}

void calin::provenance::chronicle::
register_network_close(calin::ix::provenance::chronicle::NetworkIORecord* record,
  int64_t nbytes_received, int64_t nbytes_sent)
{
  calin::util::timestamp::Timestamp ts = calin::util::timestamp::Timestamp::now();
  ts.as_proto(record->mutable_close_timestamp());
  if(nbytes_received>=0)record->set_nbytes_received(nbytes_received);
  if(nbytes_sent>=0)record->set_nbytes_sent(nbytes_sent);
}

calin::ix::provenance::chronicle::RNGRecord* calin::provenance::chronicle::
register_calin_rng_open(const calin::ix::math::rng::RNGData& rng_data,
  const std::string& created_by, const std::string& comment)
{
  calin::util::timestamp::Timestamp ts = calin::util::timestamp::Timestamp::now();
  calin::ix::provenance::chronicle::RNGRecord* record = nullptr;
  {
    std::lock_guard<std::mutex> lock { chronicle_mutex };
    record = singleton_chronicle->add_rng_record();
  }
  ts.as_proto(record->mutable_open_timestamp());
  record->mutable_calin_rng()->CopyFrom(rng_data);
  record->set_comment(comment);
  record->set_created_by(created_by);
  return record;
}

calin::ix::provenance::chronicle::RNGRecord* calin::provenance::chronicle::
register_rng_core_open(const calin::ix::math::rng::RNGCoreData& rng_core_data,
  const std::string& created_by, const std::string& comment)
{
  calin::util::timestamp::Timestamp ts = calin::util::timestamp::Timestamp::now();
  calin::ix::provenance::chronicle::RNGRecord* record = nullptr;
  {
    std::lock_guard<std::mutex> lock { chronicle_mutex };
    record = singleton_chronicle->add_rng_record();
  }
  ts.as_proto(record->mutable_open_timestamp());
  record->mutable_rng_core()->CopyFrom(rng_core_data);
  record->set_comment(comment);
  record->set_created_by(created_by);
  return record;
}

calin::ix::provenance::chronicle::RNGRecord* calin::provenance::chronicle::
register_vcl_rng_core_open(const calin::ix::math::rng::VCLRNGCoreData& vcl_rng_core_data,
  const std::string& created_by, const std::string& comment)
{
  calin::util::timestamp::Timestamp ts = calin::util::timestamp::Timestamp::now();
  calin::ix::provenance::chronicle::RNGRecord* record = nullptr;
  {
    std::lock_guard<std::mutex> lock { chronicle_mutex };
    record = singleton_chronicle->add_rng_record();
  }
  ts.as_proto(record->mutable_open_timestamp());
  record->mutable_vcl_rng_core()->CopyFrom(vcl_rng_core_data);
  record->set_comment(comment);
  record->set_created_by(created_by);
  return record;
}

calin::ix::provenance::chronicle::RNGRecord* calin::provenance::chronicle::
register_external_rng_open(uint64_t seed, const std::string& rng_type,
  const std::string& created_by, void* rng_data, unsigned rng_data_size,
  const std::string& comment)
{
  calin::util::timestamp::Timestamp ts = calin::util::timestamp::Timestamp::now();
  calin::ix::provenance::chronicle::RNGRecord* record = nullptr;
  {
    std::lock_guard<std::mutex> lock { chronicle_mutex };
    record = singleton_chronicle->add_rng_record();
  }
  ts.as_proto(record->mutable_open_timestamp());
  auto* rng_data_proto = record->mutable_external_rng();
  rng_data_proto->set_seed(seed);
  rng_data_proto->set_rng_type(rng_type);
  if(rng_data!=nullptr and rng_data_size>0)
    rng_data_proto->set_rng_data(rng_data, rng_data_size);
  record->set_comment(comment);
  record->set_created_by(created_by);
  return record;
}

void calin::provenance::chronicle::
register_rng_close(calin::ix::provenance::chronicle::RNGRecord* record, uint64_t ncore_calls)
{
  calin::util::timestamp::Timestamp ts = calin::util::timestamp::Timestamp::now();
  ts.as_proto(record->mutable_close_timestamp());
  if(ncore_calls>=0)record->set_ncore_calls(ncore_calls);
}

calin::ix::provenance::chronicle::CommandLineProcessingRecord*
calin::provenance::chronicle::register_command_line_processing(
  const std::string& processed_by, const std::string& comment)
{
  calin::util::timestamp::Timestamp ts = calin::util::timestamp::Timestamp::now();
  calin::ix::provenance::chronicle::CommandLineProcessingRecord* record = nullptr;
  {
    std::lock_guard<std::mutex> lock { chronicle_mutex };
    record = singleton_chronicle->add_command_line_record();
  }
  ts.as_proto(record->mutable_timestamp());
  record->set_comment(comment);
  record->set_processed_by(processed_by);
  return record;
}

calin::ix::provenance::chronicle::ProcessingRecord*
calin::provenance::chronicle::register_subprocessing_start(
  calin::ix::provenance::chronicle::ProcessingRecord* parent_processing_record,
  const std::string& type, const std::string& description,
  const std::string& created_by, const std::string& comment)
{
  calin::util::timestamp::Timestamp ts = calin::util::timestamp::Timestamp::now();
  calin::ix::provenance::chronicle::ProcessingRecord* record;
  if(parent_processing_record == nullptr) {
    std::lock_guard<std::mutex> lock { chronicle_mutex };
    record = singleton_chronicle->add_processing_record();
  } else {
    record = parent_processing_record->add_sub_processing_record();
  }
  ts.as_proto(record->mutable_open_timestamp());
  record->set_type(type);
  record->set_description(description);
  record->set_created_by(created_by);
  record->set_comment(comment);
  return record;
}

calin::ix::provenance::chronicle::ProcessingRecord*
calin::provenance::chronicle::register_processing_start(
  const std::string& type, const std::string& description,
  const std::string& created_by, const std::string& comment)
{
  return register_subprocessing_start(nullptr, type, description, created_by, comment);
}

void calin::provenance::chronicle::register_processing_finish(
  calin::ix::provenance::chronicle::ProcessingRecord* record)
{
  calin::util::timestamp::Timestamp ts = calin::util::timestamp::Timestamp::now();
  ts.as_proto(record->mutable_close_timestamp());
}
