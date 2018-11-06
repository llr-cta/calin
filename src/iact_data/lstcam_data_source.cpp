/*

   calin/iact_data/lstcam_data_source.cpp -- Stephen Fegan -- 2016-01-11

   A supplier of single telescope data from LSTCam DAQ data files

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

#include <stdexcept>
#include <string>
#include <memory>
#include <numeric>

#include <util/log.hpp>
#include <util/file.hpp>
#include <iact_data/lstcam_data_source.hpp>
#include <iact_data/zfits_data_source.hpp>
#include <iact_data/lstcam_layout.hpp>
#include <iact_data/lstcam_actl_event_decoder.hpp>
#include <math/simd.hpp>
#include <provenance/system_info.hpp>

using namespace calin::iact_data::lstcam_data_source;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::ix::iact_data::lstcam_data_source;
using namespace calin::util::log;
using namespace calin::util::file;
using namespace calin::iact_data::lstcam_actl_event_decoder;

#include <ProtobufIFits.h>
#include <IFits.h>
#include <R1.pb.h>

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

LSTCamZFITSDataSource_R1::
LSTCamZFITSDataSource_R1(const std::string& filename,
  const config_type& config, const decoder_config_type& decoder_config):
  calin::iact_data::zfits_data_source::ZFITSDataSource_R1(filename,
    decoder_ = new LSTCam_ACTL_R1_CameraEventDecoder(filename,
      calin::util::file::extract_run_number_from_filename(filename),
      decoder_config), false /* we delete it! */, config)
{
  // nothing to see here
}

LSTCamZFITSDataSource_R1::
LSTCamZFITSDataSource_R1(const std::string& filename,
  const decoder_config_type& decoder_config, const config_type& config):
    LSTCamZFITSDataSource_R1(filename, config, decoder_config)
{
  // nothing to see here
}

LSTCamZFITSDataSource_R1::~LSTCamZFITSDataSource_R1()
{
  delete decoder_;
}
