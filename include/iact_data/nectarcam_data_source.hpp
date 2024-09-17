/*

   calin/iact_data/nectarcam_data_source.hpp -- Stephen Fegan -- 2016-01-11

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

#pragma once

#include <string>

#include <calin_global_definitions.hpp>
#include <calin_global_config.hpp>
#include <iact_data/nectarcam_acada_event_decoder.hpp>
#include <iact_data/nectarcam_data_source.pb.h>
#include <iact_data/zfits_data_source.hpp>
#include <pattern/delegation.hpp>

namespace calin { namespace iact_data { namespace nectarcam_data_source {

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

class NectarCamZFITSDataSource_L0:
  public calin::iact_data::zfits_data_source::ZFITSDataSource<
    calin::iact_data::acada_data_source::ACADA_MessageSet_L0>
{
public:
  CALIN_TYPEALIAS(config_type, 
    ZFITSDataSource<calin::iact_data::acada_data_source::ACADA_MessageSet_L0>::config_type);
  CALIN_TYPEALIAS(decoder_config_type,
    calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig);

  calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig decoder_config() const {
    return decoder_->config(); }

  static calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig default_decoder_config() {
    return nectarcam_acada_event_decoder::NectarCam_ACADACameraEventDecoder_L0::default_config(); }

  NectarCamZFITSDataSource_L0(const std::string& filename,
    const config_type& config,
    const decoder_config_type& decoder_config = default_decoder_config());
  NectarCamZFITSDataSource_L0(const std::string& filename,
    const decoder_config_type& decoder_config = default_decoder_config(),
    const config_type& config = default_config());
  virtual ~NectarCamZFITSDataSource_L0();
private:
  nectarcam_acada_event_decoder::NectarCam_ACADACameraEventDecoder_L0* decoder_;
};

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

class NectarCamZFITSDataSource_R1v0:
  public calin::iact_data::zfits_data_source::ZFITSDataSource<
    calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0>
{
public:
  CALIN_TYPEALIAS(config_type, 
    ZFITSDataSource<calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0>::config_type);
  CALIN_TYPEALIAS(decoder_config_type,
    calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig);

  calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig decoder_config() const {
    return decoder_->config(); }

  static calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig default_decoder_config() {
    return nectarcam_acada_event_decoder::NectarCam_ACADACameraEventDecoder_R1v0::default_config(); }

  NectarCamZFITSDataSource_R1v0(const std::string& filename,
    calin::iact_data::zfits_acada_data_source::
      ZFITSACADACameraEventDataSource<calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0>* acada_zfits,
    const decoder_config_type& decoder_config = default_decoder_config(),
    bool adopt_acada_zfits = false);
  NectarCamZFITSDataSource_R1v0(const std::string& filename,
    const config_type& config,
    const decoder_config_type& decoder_config = default_decoder_config());
  NectarCamZFITSDataSource_R1v0(const std::string& filename,
    const decoder_config_type& decoder_config = default_decoder_config(),
    const config_type& config = default_config());
  virtual ~NectarCamZFITSDataSource_R1v0();
private:
  nectarcam_acada_event_decoder::NectarCam_ACADACameraEventDecoder_R1v0* decoder_;
};

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

class NectarCamZFITSDataSource:
  public calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig,
  public calin::pattern::delegation::Delegator<
    calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig>
{
public:
  CALIN_TYPEALIAS(decoder_config_type,
    calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig);
  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig);

  static calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig default_config() {
    config_type config = NectarCamZFITSDataSource_R1v0::default_config();
    config.set_data_model(calin::ix::iact_data::zfits_data_source::ACADA_DATA_MODEL_AUTO_DETECT);
    config.set_run_header_table_name(""); // Differs between L0 and R1 so let downstream decode
    return config;
  }

  static calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig default_decoder_config() {
    return NectarCamZFITSDataSource_R1v0::default_decoder_config(); }

  NectarCamZFITSDataSource(const std::string& filename,
    const calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig& config,
    const calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig& decoder_config = default_decoder_config());
  NectarCamZFITSDataSource(const std::string& filename,
    const calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig& decoder_config = default_decoder_config(),
    const calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig& config = default_config());
  virtual ~NectarCamZFITSDataSource();

  calin::ix::iact_data::telescope_event::TelescopeEvent* get_next(
    uint64_t& seq_index_out, google::protobuf::Arena** arena = nullptr) override;
  uint64_t size() override;
  void set_next_index(uint64_t next_index) override;

  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* get_run_configuration() override;

private:
  static TelescopeRandomAccessDataSourceWithRunConfig* construct_delegate(
    std::string filename, const config_type& config,
    const decoder_config_type& decoder_config);
};

} } } // namespace calin::iact_data::nectarcam_data_source
