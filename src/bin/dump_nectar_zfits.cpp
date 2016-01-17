/*

   calin/bin/dump_nectar_zfits.cpp -- Stephen Fegan -- 2016-01-14

   Dump telescope events from NectarCam data source.

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <iostream>
#include <string>
#include <cstdlib>
#include <cassert>

#include <google/protobuf/text_format.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include <iact_data/nectarcam_data_source.hpp>

using namespace calin::iact_data::nectarcam_data_source;

int main(int argc, char **argv)
{
  std::string progname(*argv);
  argv++,argc--;
  if(argc==0)
  {
    std::cout << "Usage: " << progname << " filename [nsample] [nevent_max]\n";
    exit(EXIT_FAILURE);
  }

  std::string filename(*argv);
  argv++,argc--;

  NectarCamZFITSDataSource source(filename);

  if(argc)
  {
    source.mutable_config()->set_demand_nsample(std::atoi(*argv));
    argv++,argc--;
  }

  if(argc)
  {
    source.mutable_config()->set_num_events_max(std::atoi(*argv));
    argv++,argc--;
  }

  google::protobuf::io::OstreamOutputStream stream(&std::cout);
  while(auto* event = source.getNext())
  {
    google::protobuf::TextFormat::Print(*event, &stream);
    delete event;
  }
}
