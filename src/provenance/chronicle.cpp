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
