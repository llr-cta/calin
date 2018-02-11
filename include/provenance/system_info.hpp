/*

   calin/provenance/system_info.hpp -- Stephen Fegan -- 2017-11-06

   Provenance information about build-time and run-time system environment

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <util/log.hpp>
#include <provenance/system_info.pb.h>

namespace calin { namespace provenance { namespace system_info {

const calin::ix::provenance::system_info::BuildInfo* the_build_info();
const calin::ix::provenance::system_info::HostAndProcessInfo* the_host_info();

inline const calin::ix::provenance::system_info::BuildInfo* build_info() {
  // Backward compatibility function
  return the_build_info();
}

#ifndef SWIG
calin::ix::provenance::system_info::BuildInfo*
  copy_the_build_info(calin::ix::provenance::system_info::BuildInfo* x = nullptr);
calin::ix::provenance::system_info::HostAndProcessInfo*
  copy_the_host_info(calin::ix::provenance::system_info::HostAndProcessInfo* x = nullptr);
#else
calin::ix::provenance::system_info::BuildInfo* copy_the_build_info();
void copy_the_build_info(calin::ix::provenance::system_info::BuildInfo* x);
calin::ix::provenance::system_info::HostAndProcessInfo* copy_the_host_info();
void copy_the_host_info(calin::ix::provenance::system_info::HostAndProcessInfo* x);
#endif

void write_system_info_to_log(calin::util::log::Level level = calin::util::log::INFO,
  calin::util::log::Logger* logger = calin::util::log::default_logger());

} } } // namespace calin::provenance::system
