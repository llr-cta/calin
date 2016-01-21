/*

   calin/bin/write_nectar_zfits_as_raw.cpp -- Stephen Fegan -- 2016-01-17

   Write telescope events from NectarCam data source in a raw file.

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

#include <iact_data/telescope_data_source.hpp>
#include <iact_data/nectarcam_data_source.hpp>

using namespace calin::iact_data::telescope_data_source;
using namespace calin::iact_data::nectarcam_data_source;

int main(int argc, char **argv)
{
  std::string progname(*argv);
  argv++,argc--;
  if(argc==0)
  {
    std::cout << "Usage: " << progname
      << " in_filename out_filename [nsample] [nevent_max]\n";
    exit(EXIT_FAILURE);
  }

  std::string in_filename(*argv);
  argv++,argc--;

  std::string out_filename(*argv);
  argv++,argc--;

  NectarCamZFITSDataSource source(in_filename);

  if(argc)
  {
    source.mutable_decoder_config()->set_demand_nsample(std::atoi(*argv));
    argv++,argc--;
  }

  if(argc)
  {
    source.mutable_reader_config()->set_num_events_max(std::atoi(*argv));
    argv++,argc--;
  }

  RawFileTelescopeDataSink::config_type config;
  //config.set_use_compression(true);
  RawFileTelescopeDataSink sink(out_filename,false,config);

  unsigned nevent = 0;
  while(auto* event = source.get_next())
  {
    sink.put_next(event);
    nevent++;
  }
  std::cout << "Copied " << nevent << " events\n";
}
