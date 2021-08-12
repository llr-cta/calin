/*

   calin/provenance/chronicle.hpp -- Stephen Fegan -- 2017-10-30

   Provenance system for recording noteworthy run-time events (io files etc)

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#pragma once

#include <string>

#include <provenance/chronicle.pb.h>

namespace calin { namespace provenance { namespace chronicle {

const calin::ix::provenance::chronicle::Chronicle* the_chronicle();

void prune_the_chronicle();

#ifndef SWIG
calin::ix::provenance::chronicle::Chronicle*
  copy_the_chronicle(calin::ix::provenance::chronicle::Chronicle* x = nullptr);
#else
calin::ix::provenance::chronicle::Chronicle* copy_the_chronicle();
void copy_the_chronicle(calin::ix::provenance::chronicle::Chronicle* x);
#endif

calin::ix::provenance::chronicle::FileIORecord*
register_file_open(const std::string& file_name,
  calin::ix::provenance::chronicle::AccessType access,
  const std::string& opened_by, const std::string& comment = "");
void register_file_close(calin::ix::provenance::chronicle::FileIORecord* record);

calin::ix::provenance::chronicle::NetworkIORecord*
register_network_open(const std::string& endpoint,
  const std::string& opened_by, const std::string& comment = "");
void register_network_close(calin::ix::provenance::chronicle::NetworkIORecord* record,
  int64_t nbytes_received = -1, int64_t nbytes_sent = -1);

calin::ix::provenance::chronicle::RNGRecord*
register_calin_rng_open(const calin::ix::math::rng::RNGData& rng_data,
  const std::string& created_by, const std::string& comment = "");
calin::ix::provenance::chronicle::RNGRecord*
register_rng_core_open(const calin::ix::math::rng::RNGCoreData& rng_data,
  const std::string& created_by, const std::string& comment = "");
calin::ix::provenance::chronicle::RNGRecord*
register_vcl_rng_core_open(const calin::ix::math::rng::VCLRNGCoreData& vcl_rng_data,
  const std::string& created_by, const std::string& comment = "");
calin::ix::provenance::chronicle::RNGRecord*
register_external_rng_open(uint64_t seed, const std::string& rng_type,
  const std::string& created_by, void* rng_data = nullptr, unsigned rng_data_size = 0,
  const std::string& comment = "");
void register_rng_close(calin::ix::provenance::chronicle::RNGRecord* record,
  uint64_t ncore_calls=0);

calin::ix::provenance::chronicle::CommandLineProcessingRecord*
register_command_line_processing(const std::string& processed_by, const std::string& comment = "");

calin::ix::provenance::chronicle::ProcessingRecord*
register_subprocessing_start(
  calin::ix::provenance::chronicle::ProcessingRecord* parent_processing_record,
  const std::string& type, const std::string& description,
  const std::string& created_by, const std::string& comment = "");
calin::ix::provenance::chronicle::ProcessingRecord*
register_processing_start(const std::string& type, const std::string& description,
  const std::string& created_by, const std::string& comment = "");
void register_processing_finish(calin::ix::provenance::chronicle::ProcessingRecord* record);

//void register_rng_creation();

} } } // namespace calin::provenance::system
