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
using namespace calin::util::file;
using namespace calin::iact_data::nectarcam_actl_event_decoder;

#include <ProtobufIFits.h>
#include <IFits.h>
#include <L0.pb.h>
#include <R1.pb.h>

/*

              LLLLLLLLLLL                       000000000
              L:::::::::L                     00:::::::::00
              L:::::::::L                   00:::::::::::::00
              LL:::::::LL                  0:::::::000:::::::0
                L:::::L                    0::::::0   0::::::0
                L:::::L                    0:::::0     0:::::0
                L:::::L                    0:::::0     0:::::0
                L:::::L                    0:::::0 000 0:::::0
                L:::::L                    0:::::0 000 0:::::0
                L:::::L                    0:::::0     0:::::0
                L:::::L                    0:::::0     0:::::0
                L:::::L         LLLLLL     0::::::0   0::::::0
              LL:::::::LLLLLLLLL:::::L     0:::::::000:::::::0
              L::::::::::::::::::::::L      00:::::::::::::00
              L::::::::::::::::::::::L        00:::::::::00
              LLLLLLLLLLLLLLLLLLLLLLLL          000000000

*/

NectarCamZFITSDataSource_L0::
NectarCamZFITSDataSource_L0(const std::string& filename,
  const config_type& config, const decoder_config_type& decoder_config):
  calin::iact_data::zfits_data_source::ZFITSDataSource_L0(filename,
    decoder_ = new NectarCAM_ACTL_L0_CameraEventDecoder(filename,
      calin::util::file::extract_first_number_from_filename(filename),
      decoder_config), false /* we delete it! */, config)
{
  // nothing to see here
}

NectarCamZFITSDataSource_L0::
NectarCamZFITSDataSource_L0(const std::string& filename,
  const decoder_config_type& decoder_config, const config_type& config):
    NectarCamZFITSDataSource_L0(filename, config, decoder_config)
{
  // nothing to see here
}

NectarCamZFITSDataSource_L0::~NectarCamZFITSDataSource_L0()
{
  delete decoder_;
}

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

NectarCamZFITSDataSource_R1::
NectarCamZFITSDataSource_R1(const std::string& filename,
  const config_type& config, const decoder_config_type& decoder_config):
  calin::iact_data::zfits_data_source::ZFITSDataSource_R1(filename,
    decoder_ = new NectarCAM_ACTL_R1_CameraEventDecoder(filename,
      calin::util::file::extract_first_number_from_filename(filename),
      decoder_config), false /* we delete it! */, config)
{
  // nothing to see here
}

NectarCamZFITSDataSource_R1::
NectarCamZFITSDataSource_R1(const std::string& filename,
  const decoder_config_type& decoder_config, const config_type& config):
    NectarCamZFITSDataSource_R1(filename, config, decoder_config)
{
  // nothing to see here
}

NectarCamZFITSDataSource_R1::~NectarCamZFITSDataSource_R1()
{
  delete decoder_;
}

/*

               AAA                                     tttt
              A:::A                                 ttt:::t
             A:::::A                                t:::::t
            A:::::::A                               t:::::t
           A:::::::::A        uuuuuu    uuuuuuttttttt:::::ttttttt       ooooooooooo
          A:::::A:::::A       u::::u    u::::ut:::::::::::::::::t     oo:::::::::::oo
         A:::::A A:::::A      u::::u    u::::ut:::::::::::::::::t    o:::::::::::::::o
        A:::::A   A:::::A     u::::u    u::::utttttt:::::::tttttt    o:::::ooooo:::::o
       A:::::A     A:::::A    u::::u    u::::u      t:::::t          o::::o     o::::o
      A:::::AAAAAAAAA:::::A   u::::u    u::::u      t:::::t          o::::o     o::::o
     A:::::::::::::::::::::A  u::::u    u::::u      t:::::t          o::::o     o::::o
    A:::::AAAAAAAAAAAAA:::::A u:::::uuuu:::::u      t:::::t    tttttto::::o     o::::o
   A:::::A             A:::::Au:::::::::::::::uu    t::::::tttt:::::to:::::ooooo:::::o
  A:::::A               A:::::Au:::::::::::::::u    tt::::::::::::::to:::::::::::::::o
 A:::::A                 A:::::Auu::::::::uu:::u      tt:::::::::::tt oo:::::::::::oo
AAAAAAA                   AAAAAAA uuuuuuuu  uuuu        ttttttttttt     ooooooooooo

*/

NectarCamZFITSDataSource::
NectarCamZFITSDataSource(const std::string& filename,
    const config_type& config, const decoder_config_type& decoder_config):
  TelescopeRandomAccessDataSourceWithRunConfig(),
  Delegator<calin::iact_data::telescope_data_source::
    TelescopeRandomAccessDataSourceWithRunConfig>(
      construct_delegate(filename, config, decoder_config), /* adopt_delegate=*/ true)
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
  // nothing to see here
}

calin::ix::iact_data::telescope_event::TelescopeEvent*
NectarCamZFITSDataSource::get_next(uint64_t& seq_index_out, google::protobuf::Arena** arena)
{
  return delegate_->get_next(seq_index_out, arena);
}

uint64_t NectarCamZFITSDataSource::size()
{
  return delegate_->size();
}

void NectarCamZFITSDataSource::set_next_index(uint64_t next_index)
{
  return delegate_->set_next_index(next_index);
}

calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration*
NectarCamZFITSDataSource::get_run_configuration()
{
  return delegate_->get_run_configuration();
}

calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig*
NectarCamZFITSDataSource::construct_delegate(const std::string& filename,
  const config_type& config, const decoder_config_type& decoder_config)
{
  bool use_r1 = true;
  if(config.data_model() ==
      calin::ix::iact_data::zfits_data_source::ACTL_DATA_MODEL_AUTO_DETECT) {
    if(!is_file(filename))
      throw(std::runtime_error(
        "NectarCamZFITSDataSource::construct_delegate: File not found: " + filename));
    std::string events_table_name = config.events_table_name();
    if(events_table_name == "")events_table_name =
      zfits_actl_data_source::ZFITSSingleFileACTL_R1_CameraEventDataSource::
        default_config().events_table_name();
    IFits ifits(filename, events_table_name, /* force= */ true);
    const std::string& message_name = ifits.GetStr("PBFHEAD");
    if (message_name == "DataModel.CameraEvent") {
      use_r1 = false;
    } else if (message_name == "R1.CameraEvent") {
      use_r1 = true;
    } else {
      throw(std::runtime_error(
        "NectarCamZFITSDataSource::construct_delegate: Unknown message type: " + message_name));
    }
  } else if(config.data_model() ==
      calin::ix::iact_data::zfits_data_source::ACTL_DATA_MODEL_L0) {
    use_r1 = false;
  } else if(config.data_model() ==
      calin::ix::iact_data::zfits_data_source::ACTL_DATA_MODEL_R1) {
    use_r1 = true;
  } else {
    throw(std::runtime_error("NectarCamZFITSDataSource::construct_delegate: Requested data model not unknown"));
  }

  if(use_r1)return new NectarCamZFITSDataSource_R1(filename, config, decoder_config);
  else return new NectarCamZFITSDataSource_L0(filename, config, decoder_config);
}
