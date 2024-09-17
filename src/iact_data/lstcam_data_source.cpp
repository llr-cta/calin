/*

   calin/iact_data/lstcam_data_source.cpp -- Stephen Fegan -- 2016-01-11

   A supplier of single telescope data from LSTCam DAQ data files

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <stdexcept>
#include <string>
#include <memory>
#include <numeric>

#include <util/log.hpp>
#include <util/file.hpp>
#include <iact_data/lstcam_data_source.hpp>
#include <iact_data/zfits_data_source.hpp>
#include <iact_data/lstcam_layout.hpp>
#include <iact_data/lstcam_acada_event_decoder.hpp>
#include <provenance/system_info.hpp>

using namespace calin::iact_data::lstcam_data_source;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::ix::iact_data::lstcam_data_source;
using namespace calin::util::log;
using namespace calin::util::file;
using namespace calin::iact_data::lstcam_acada_event_decoder;

/*

    RRRRRRRRRRRRRRRRR     1111111                              000000000     
    R::::::::::::::::R   1::::::1                            00:::::::::00   
    R::::::RRRRRR:::::R 1:::::::1                          00:::::::::::::00 
    RR:::::R     R:::::R111:::::1                         0:::::::000:::::::0
      R::::R     R:::::R   1::::1vvvvvvv           vvvvvvv0::::::0   0::::::0
      R::::R     R:::::R   1::::1 v:::::v         v:::::v 0:::::0     0:::::0
      R::::RRRRRR:::::R    1::::1  v:::::v       v:::::v  0:::::0     0:::::0
      R:::::::::::::RR     1::::l   v:::::v     v:::::v   0:::::0 000 0:::::0
      R::::RRRRRR:::::R    1::::l    v:::::v   v:::::v    0:::::0 000 0:::::0
      R::::R     R:::::R   1::::l     v:::::v v:::::v     0:::::0     0:::::0
      R::::R     R:::::R   1::::l      v:::::v:::::v      0:::::0     0:::::0
      R::::R     R:::::R   1::::l       v:::::::::v       0::::::0   0::::::0
    RR:::::R     R:::::R111::::::111     v:::::::v        0:::::::000:::::::0
    R::::::R     R:::::R1::::::::::1      v:::::v          00:::::::::::::00 
    R::::::R     R:::::R1::::::::::1       v:::v             00:::::::::00   
    RRRRRRRR     RRRRRRR111111111111        vvv                000000000     

*/

LSTCamZFITSDataSource_R1v0::
LSTCamZFITSDataSource_R1v0(const std::string& filename,
    calin::iact_data::zfits_acada_data_source::
      ZFITSACADACameraEventDataSource<calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0>* acada_zfits,
    const decoder_config_type& decoder_config, bool adopt_acada_zfits):
  calin::iact_data::zfits_data_source::ZFITSDataSource_R1v0(
    acada_zfits, 
    decoder_ = new LSTCam_ACADACameraEventDecoder_R1v0(filename,
      calin::util::file::extract_run_number_from_filename(filename),
      decoder_config), 
    adopt_acada_zfits, /* adopt_decoder_= */ false)
{
  // nothing to see here
}

LSTCamZFITSDataSource_R1v0::
LSTCamZFITSDataSource_R1v0(const std::string& filename,
  const config_type& config, const decoder_config_type& decoder_config):
  calin::iact_data::zfits_data_source::ZFITSDataSource_R1v0(filename,
    decoder_ = new LSTCam_ACADACameraEventDecoder_R1v0(filename,
      calin::util::file::extract_run_number_from_filename(filename),
      decoder_config), false /* we delete it! */, config)
{
  // nothing to see here
}

LSTCamZFITSDataSource_R1v0::
LSTCamZFITSDataSource_R1v0(const std::string& filename,
  const decoder_config_type& decoder_config, const config_type& config):
    LSTCamZFITSDataSource_R1v0(filename, config, decoder_config)
{
  // nothing to see here
}

LSTCamZFITSDataSource_R1v0::~LSTCamZFITSDataSource_R1v0()
{
  delete decoder_;
}
