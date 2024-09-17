/*

   calin/iact_data/nectarcam_data_source.cpp -- Stephen Fegan -- 2016-01-11

   A supplier of single telescope data from NectarCam DAQ data files

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
#include <iact_data/nectarcam_data_source.hpp>
#include <iact_data/zfits_data_source.hpp>
#include <iact_data/zfits_acada_data_source.hpp>
#include <iact_data/nectarcam_layout.hpp>
#include <iact_data/nectarcam_acada_event_decoder.hpp>
#include <iact_data/nectarcam_configuration.hpp>
#include <provenance/system_info.hpp>

using namespace calin::iact_data::nectarcam_data_source;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::ix::iact_data::nectarcam_data_source;
using namespace calin::util::log;
using namespace calin::util::file;
using namespace calin::iact_data::nectarcam_acada_event_decoder;


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
    decoder_ = new NectarCam_ACADACameraEventDecoder_L0(filename,
      calin::util::file::extract_run_number_from_filename(filename),
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

NectarCamZFITSDataSource_R1v0::
NectarCamZFITSDataSource_R1v0(const std::string& filename,
    calin::iact_data::zfits_acada_data_source::
      ZFITSACADACameraEventDataSource<calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0>* acada_zfits,
    const decoder_config_type& decoder_config, bool adopt_acada_zfits):
  calin::iact_data::zfits_data_source::ZFITSDataSource_R1v0(
    acada_zfits, 
    decoder_ = new NectarCam_ACADACameraEventDecoder_R1v0(filename,
      calin::util::file::extract_run_number_from_filename(filename),
      decoder_config), 
    adopt_acada_zfits, /* adopt_decoder_= */ false)
{
  // nothing to see here
}

NectarCamZFITSDataSource_R1v0::
NectarCamZFITSDataSource_R1v0(const std::string& filename,
    const config_type& config, const decoder_config_type& decoder_config):
  calin::iact_data::zfits_data_source::ZFITSDataSource_R1v0(filename,
    decoder_ = new NectarCam_ACADACameraEventDecoder_R1v0(filename,
      calin::util::file::extract_run_number_from_filename(filename),
      decoder_config), false /* we delete it! */, config)
{
  // nothing to see here
}

NectarCamZFITSDataSource_R1v0::
NectarCamZFITSDataSource_R1v0(const std::string& filename,
  const decoder_config_type& decoder_config, const config_type& config):
    NectarCamZFITSDataSource_R1v0(filename, config, decoder_config)
{
  // nothing to see here
}

NectarCamZFITSDataSource_R1v0::~NectarCamZFITSDataSource_R1v0()
{
  delete decoder_;
}

/*

        RRRRRRRRRRRRRRRRR     1111111                        1111111   
        R::::::::::::::::R   1::::::1                       1::::::1   
        R::::::RRRRRR:::::R 1:::::::1                      1:::::::1   
        RR:::::R     R:::::R111:::::1                      111:::::1   
          R::::R     R:::::R   1::::1vvvvvvv           vvvvvvv1::::1   
          R::::R     R:::::R   1::::1 v:::::v         v:::::v 1::::1   
          R::::RRRRRR:::::R    1::::1  v:::::v       v:::::v  1::::1   
          R:::::::::::::RR     1::::l   v:::::v     v:::::v   1::::l   
          R::::RRRRRR:::::R    1::::l    v:::::v   v:::::v    1::::l   
          R::::R     R:::::R   1::::l     v:::::v v:::::v     1::::l   
          R::::R     R:::::R   1::::l      v:::::v:::::v      1::::l   
          R::::R     R:::::R   1::::l       v:::::::::v       1::::l   
        RR:::::R     R:::::R111::::::111     v:::::::v     111::::::111
        R::::::R     R:::::R1::::::::::1      v:::::v      1::::::::::1
        R::::::R     R:::::R1::::::::::1       v:::v       1::::::::::1
        RRRRRRRR     RRRRRRR111111111111        vvv        111111111111
                                                               

*/

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
NectarCamZFITSDataSource::construct_delegate(std::string filename,
  const config_type& config, const decoder_config_type& decoder_config)
{
  expand_filename_in_place(filename);

  auto model = config.data_model();
  if(model == calin::ix::iact_data::zfits_data_source::ACADA_DATA_MODEL_AUTO_DETECT) {
    model = calin::iact_data::zfits_acada_data_source::get_zfits_data_model(filename);
  }

  switch(model) {
  case calin::ix::iact_data::zfits_data_source::ACADA_DATA_MODEL_AUTO_DETECT:
    throw std::runtime_error("NectarCamZFITSDataSource::construct_delegate: Requested data model not known");
  case calin::ix::iact_data::zfits_data_source::ACADA_DATA_MODEL_L0:
    return new NectarCamZFITSDataSource_L0(filename, config, decoder_config);
  case calin::ix::iact_data::zfits_data_source::ACADA_DATA_MODEL_R1V0:
    return new NectarCamZFITSDataSource_R1v0(filename, config, decoder_config);
  case calin::ix::iact_data::zfits_data_source::ACADA_DATA_MODEL_R1V1:
    throw std::runtime_error("NectarCamZFITSDataSource::construct_delegate: R1v1 not supported yet");
  default:
    throw std::runtime_error("NectarCamZFITSDataSource::construct_delegate: Requested data model not known");
  }

  return nullptr;
}
