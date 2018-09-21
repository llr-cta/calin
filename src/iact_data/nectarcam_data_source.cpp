/*

   calin/iact_data/nectarcam_data_source.cpp -- Stephen Fegan -- 2016-01-11

   A supplier of single telescope data from NectarCam DAQ data files

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
#include <iact_data/nectarcam_data_source.hpp>
#include <iact_data/zfits_data_source.hpp>
#include <iact_data/nectarcam_layout.hpp>
#include <iact_data/nectarcam_actl_event_decoder.hpp>
#include <iact_data/nectarcam_module_configuration.hpp>
#include <math/simd.hpp>
#include <provenance/system_info.hpp>

using namespace calin::iact_data::nectarcam_data_source;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::ix::iact_data::nectarcam_data_source;
using namespace calin::util::log;
using namespace calin::iact_data::nectarcam_actl_event_decoder;

#include <ProtobufIFits.h>
#include <L0.pb.h>

// =============================================================================
// =============================================================================
//
// NectarCamZFITSDataSource
//
// =============================================================================
// =============================================================================

NectarCamZFITSDataSource::
NectarCamZFITSDataSource(const std::string& filename,
  const config_type& config, const decoder_config_type& decoder_config):
  calin::iact_data::zfits_data_source::ZFITSDataSource(filename,
    decoder_ = new NectarCAM_ACTL_L0_CameraEventDecoder(filename,
      calin::util::file::extract_first_number_from_filename(filename),
      decoder_config), false /* we delete it! */, config)
{
  // nothing to see here
}

NectarCamZFITSDataSource::
NectarCamZFITSDataSource(const std::string& filename,
  const decoder_config_type& decoder_config, const config_type& config):
    NectarCamZFITSDataSource(filename, config, decoder_config)
{
  // nothing to see here
}

NectarCamZFITSDataSource::~NectarCamZFITSDataSource()
{
  delete decoder_;
}
