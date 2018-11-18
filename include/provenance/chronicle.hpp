/*

   calin/provenance/chronicle.hpp -- Stephen Fegan -- 2017-10-30

   Provenance system for recording noteworthy run-time events (io files etc)

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#pragma once

#include <string>

#include <provenance/chronicle.pb.h>

namespace calin { namespace provenance { namespace chronicle {

const calin::ix::provenance::chronicle::Chronicle* the_chronicle();

#ifndef SWIG
calin::ix::provenance::chronicle::Chronicle*
  copy_the_chronicle(calin::ix::provenance::chronicle::Chronicle* x = nullptr);
#else
calin::ix::provenance::chronicle::Chronicle* copy_the_chronicle();
void copy_the_chronicle(calin::ix::provenance::chronicle::Chronicle* x);
#endif

void register_file_open(const std::string& file_name,
  calin::ix::provenance::chronicle::AccessType access,
  const std::string& opened_by, const std::string& comment = "");
void register_network_open(const std::string& endpoint,
  const std::string& opened_by, const std::string& comment = "");

void register_calin_rng(const calin::ix::math::rng::RNGData& rng_data,
  const std::string& created_by, const std::string& comment = "");
void register_rng_core(const calin::ix::math::rng::RNGCoreData& rng_data,
  const std::string& created_by, const std::string& comment = "");
void register_vcl_rng_core(const calin::ix::math::rng::VCLRNGCoreData& vcl_rng_data,
  const std::string& created_by, const std::string& comment = "");
void register_external_rng(uint64_t seed, const std::string& rng_type,
  const std::string& created_by, void* rng_data = nullptr, unsigned rng_data_size = 0,
  const std::string& comment = "");

//void register_rng_creation();

} } } // namespace calin::provenance::system
