/*

   calin/bin/write_system_info.cpp -- Stephen Fegan -- 2018-02-11

   Write the system info to the log

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

#include <provenance/system_info.hpp>
#include <util/log.hpp>

int main(int argc, char** argv)
{
  calin::util::log::default_logger()->add_cout(/* apply_timestamp = */ false,
    /* use_colors = */ false);
  calin::provenance::system_info::write_system_info_to_log();
  if(argc>1) {
    const auto* host_info = calin::provenance::system_info::the_host_info();
    const auto* build_info = calin::provenance::system_info::the_build_info();
    auto L = calin::util::log::LOG(calin::util::log::INFO);
    L << "----------\nBUILD INFO\n----------\n" << build_info->DebugString();
    L << "\n---------\nHOST INFO\n---------\n" << host_info->DebugString();
  }
}
