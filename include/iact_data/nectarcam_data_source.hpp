/*

   calin/iact_data/nectarcam_data_source.hpp -- Stephen Fegan -- 2016-01-11

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

#pragma once

#include <string>

#include <calin_global_definitions.hpp>
#include <calin_global_config.hpp>
#include <iact_data/nectarcam_actl_event_decoder.hpp>
#include <iact_data/nectarcam_data_source.pb.h>
#include <iact_data/zfits_data_source.hpp>
#include <pattern/delegation.hpp>

namespace calin { namespace iact_data { namespace nectarcam_data_source {

#ifdef CALIN_HAVE_CTA_CAMERASTOACTL

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
  public calin::iact_data::zfits_data_source::ZFITSDataSource_L0
{
public:
  CALIN_TYPEALIAS(decoder_config_type,
    nectarcam_actl_event_decoder::NectarCam_ACTL_L0_CameraEventDecoder::config_type);

#if 0
  void set_decoder_config(const decoder_config_type& config) {
    decoder_->set_config(config); }
#endif
  decoder_config_type decoder_config() const {
    return decoder_->config(); }
#if 0
  decoder_config_type* mutable_decoder_config() {
    return decoder_->mutable_config(); }
#endif
  static decoder_config_type default_decoder_config() {
    return nectarcam_actl_event_decoder::NectarCam_ACTL_L0_CameraEventDecoder::default_config(); }

  NectarCamZFITSDataSource_L0(const std::string& filename,
    const config_type& config,
    const decoder_config_type& decoder_config = default_decoder_config());
  NectarCamZFITSDataSource_L0(const std::string& filename,
    const decoder_config_type& decoder_config = default_decoder_config(),
    const config_type& config = default_config());
  virtual ~NectarCamZFITSDataSource_L0();
private:
  nectarcam_actl_event_decoder::NectarCam_ACTL_L0_CameraEventDecoder* decoder_;
};

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

class NectarCamZFITSDataSource_R1:
  public calin::iact_data::zfits_data_source::ZFITSDataSource_R1
{
public:
  CALIN_TYPEALIAS(decoder_config_type,
    nectarcam_actl_event_decoder::NectarCam_ACTL_R1_CameraEventDecoder::config_type);

#if 0
  void set_decoder_config(const decoder_config_type& config) {
    decoder_->set_config(config); }
#endif
  decoder_config_type decoder_config() const {
    return decoder_->config(); }
#if 0
  decoder_config_type* mutable_decoder_config() {
    return decoder_->mutable_config(); }
#endif
  static decoder_config_type default_decoder_config() {
    return nectarcam_actl_event_decoder::NectarCam_ACTL_R1_CameraEventDecoder::default_config(); }

  NectarCamZFITSDataSource_R1(const std::string& filename,
    const config_type& config,
    const decoder_config_type& decoder_config = default_decoder_config());
  NectarCamZFITSDataSource_R1(const std::string& filename,
    const decoder_config_type& decoder_config = default_decoder_config(),
    const config_type& config = default_config());
  virtual ~NectarCamZFITSDataSource_R1();
private:
  nectarcam_actl_event_decoder::NectarCam_ACTL_R1_CameraEventDecoder* decoder_;
};

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

  static config_type default_config() {
    config_type config = NectarCamZFITSDataSource_R1::default_config();
    config.set_data_model(calin::ix::iact_data::zfits_data_source::ACTL_DATA_MODEL_AUTO_DETECT);
    config.set_run_header_table_name(""); // Differs between L0 and R1 so let downstream decode
    return config;
  }

  static decoder_config_type default_decoder_config() {
    return NectarCamZFITSDataSource_R1::default_decoder_config(); }

  NectarCamZFITSDataSource(const std::string& filename,
    const config_type& config,
    const decoder_config_type& decoder_config = default_decoder_config());
  NectarCamZFITSDataSource(const std::string& filename,
    const decoder_config_type& decoder_config = default_decoder_config(),
    const config_type& config = default_config());
  virtual ~NectarCamZFITSDataSource();

  calin::ix::iact_data::telescope_event::TelescopeEvent* get_next(
    uint64_t& seq_index_out, google::protobuf::Arena** arena = nullptr) override;
  uint64_t size() override;
  void set_next_index(uint64_t next_index) override;

  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* get_run_configuration() override;

private:
  static TelescopeRandomAccessDataSourceWithRunConfig* construct_delegate(
    const std::string& filename, const config_type& config,
    const decoder_config_type& decoder_config);
};

#endif

} } } // namespace calin::iact_data::nectarcam_data_source
