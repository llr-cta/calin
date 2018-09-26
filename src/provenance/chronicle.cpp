/*

   calin/provenance/chronicle.cpp -- Stephen Fegan -- 2017-10-30

   Provenance system for recording noteworthy run-time events (io files etc)

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

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

} // anonymous namespace

const calin::ix::provenance::chronicle::Chronicle*
calin::provenance::chronicle::the_chronicle()
{
  return singleton_chronicle.get();
}

calin::ix::provenance::chronicle::Chronicle*
calin::provenance::chronicle::copy_the_chronicle(calin::ix::provenance::chronicle::Chronicle* x)
{
  if(x == nullptr) x = new calin::ix::provenance::chronicle::Chronicle;
  std::lock_guard<std::mutex> lock { chronicle_mutex };
  x->CopyFrom(*the_chronicle());
  return x;
}

void calin::provenance::chronicle::
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
  ts.as_proto(record->mutable_timestamp());
  record->set_access(access);
  record->set_file_name(file_name);
  record->set_comment(comment);
  struct stat stat_buffer;
  if(::stat(file_name.c_str(), &stat_buffer) >= 0) {
    record->set_file_size(stat_buffer.st_size);
    calin::util::timestamp::Timestamp(stat_buffer.st_mtime).as_proto(record->mutable_file_mtime());
  }
  record->set_opened_by(opened_by);
}

void calin::provenance::chronicle::
register_calin_rng(const calin::ix::math::rng::RNGData& rng_data,
  const std::string& created_by, const std::string& comment)
{
  calin::util::timestamp::Timestamp ts = calin::util::timestamp::Timestamp::now();
  calin::ix::provenance::chronicle::RNGRecord* record = nullptr;
  {
    std::lock_guard<std::mutex> lock { chronicle_mutex };
    record = singleton_chronicle->add_rng_record();
  }
  ts.as_proto(record->mutable_timestamp());
  record->mutable_calin_rng()->CopyFrom(rng_data);
  record->set_comment(comment);
  record->set_created_by(created_by);
}

void calin::provenance::chronicle::
register_rng_core(const calin::ix::math::rng::RNGCoreData& rng_core_data,
  const std::string& created_by, const std::string& comment)
{
  calin::util::timestamp::Timestamp ts = calin::util::timestamp::Timestamp::now();
  calin::ix::provenance::chronicle::RNGCoreRecord* record = nullptr;
  {
    std::lock_guard<std::mutex> lock { chronicle_mutex };
    record = singleton_chronicle->add_rng_core_record();
  }
  ts.as_proto(record->mutable_timestamp());
  record->mutable_rng_core()->CopyFrom(rng_core_data);
  record->set_comment(comment);
  record->set_created_by(created_by);
}

void calin::provenance::chronicle::
register_vcl_rng_core(const calin::ix::math::rng::VCLRNGCoreData& vcl_rng_core_data,
  const std::string& created_by, const std::string& comment)
{
  calin::util::timestamp::Timestamp ts = calin::util::timestamp::Timestamp::now();
  calin::ix::provenance::chronicle::VCLRNGCoreRecord* record = nullptr;
  {
    std::lock_guard<std::mutex> lock { chronicle_mutex };
    record = singleton_chronicle->add_vcl_rng_core_record();
  }
  ts.as_proto(record->mutable_timestamp());
  record->mutable_vcl_rng_core()->CopyFrom(vcl_rng_core_data);
  record->set_comment(comment);
  record->set_created_by(created_by);
}

void calin::provenance::chronicle::
register_external_rng(uint64_t seed, const std::string& rng_type,
  const std::string& created_by, void* rng_data, unsigned rng_data_size,
  const std::string& comment)
{
  calin::util::timestamp::Timestamp ts = calin::util::timestamp::Timestamp::now();
  calin::ix::provenance::chronicle::RNGRecord* record = nullptr;
  {
    std::lock_guard<std::mutex> lock { chronicle_mutex };
    record = singleton_chronicle->add_rng_record();
  }
  ts.as_proto(record->mutable_timestamp());
  auto* rng_data_proto = record->mutable_external_rng();
  rng_data_proto->set_seed(seed);
  rng_data_proto->set_rng_type(rng_type);
  if(rng_data!=nullptr and rng_data_size>0)
    rng_data_proto->set_rng_data(rng_data, rng_data_size);
  record->set_comment(comment);
  record->set_created_by(created_by);
}
