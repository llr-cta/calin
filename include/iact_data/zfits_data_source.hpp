/*

   calin/iact_data/zfits_data_source.hpp -- Stephen Fegan -- 2016-01-21

   A supplier of single telescope data from CTA ACTL ZFits data files

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
#include <io/chained_data_source.hpp>
#include <iact_data/actl_event_decoder.hpp>
#include <iact_data/zfits_actl_data_source.hpp>
#include <iact_data/zfits_data_source.pb.h>
#include <iact_data/telescope_data_source.hpp>

#ifdef CALIN_HAVE_CTA_CAMERASTOACTL
#include <CoreMessages.pb.h>
#endif

namespace calin { namespace iact_data { namespace zfits_data_source {

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

class ZFITSSingleFileDataSource_L0:
  public calin::iact_data::telescope_data_source::
    TelescopeRandomAccessDataSourceWithRunConfig
{
public:
  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig);

  ZFITSSingleFileDataSource_L0(calin::iact_data::zfits_actl_data_source::
      ZFITSSingleFileACTL_L0_CameraEventDataSource* actl_zfits,
    bool dont_decode_run_configuration,
    calin::iact_data::actl_event_decoder::ACTL_L0_CameraEventDecoder* decoder,
    bool adopt_decoder = false,
    bool adopt_actl_zfits = false);

  ZFITSSingleFileDataSource_L0(const std::string& filename,
    calin::iact_data::actl_event_decoder::ACTL_L0_CameraEventDecoder* decoder,
    bool adopt_decoder = false,
    const config_type& config = default_config());
  virtual ~ZFITSSingleFileDataSource_L0();

  calin::ix::iact_data::telescope_event::TelescopeEvent* get_next(
    uint64_t& seq_index_out, google::protobuf::Arena** arena = nullptr) override;
  uint64_t size() override;
  void set_next_index(uint64_t next_index) override;

  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* get_run_configuration() override;

  static config_type default_config() {
    return zfits_actl_data_source::ZFITSSingleFileACTL_L0_CameraEventDataSource::default_config(); }

private:
  calin::iact_data::actl_event_decoder::ACTL_L0_CameraEventDecoder* decoder_;
  bool adopt_decoder_ = false;
  calin::iact_data::zfits_actl_data_source::
    ZFITSSingleFileACTL_L0_CameraEventDataSource* actl_zfits_ = nullptr;
  bool adopt_actl_zfits_ = false;
  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* run_config_ = nullptr;
};

class ZFITSDataSource_L0:
  public calin::io::data_source::BasicChainedRandomAccessDataSource<
    calin::iact_data::telescope_data_source::
      TelescopeRandomAccessDataSourceWithRunConfig>
{
public:
  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig);

  ZFITSDataSource_L0(const std::string& filename,
    calin::iact_data::actl_event_decoder::ACTL_L0_CameraEventDecoder* decoder,
    bool adopt_decoder = false,
    const config_type& config = default_config());
  virtual ~ZFITSDataSource_L0();

  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* get_run_configuration() override;

  static config_type default_config() {
    return zfits_actl_data_source::ZFITSACTL_L0_CameraEventDataSource::default_config(); }

protected:
  calin::iact_data::actl_event_decoder::ACTL_L0_CameraEventDecoder* decoder_;
  bool adopt_decoder_ = false;
  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* run_config_ = nullptr;
};

class ZFITSDataSourceOpener_L0:
  public calin::io::data_source::DataSourceOpener<
    calin::iact_data::telescope_data_source::
      TelescopeRandomAccessDataSourceWithRunConfig>
{
public:
  CALIN_TYPEALIAS(data_source_type, calin::iact_data::telescope_data_source::
    TelescopeRandomAccessDataSourceWithRunConfig);
  ZFITSDataSourceOpener_L0(std::string filename,
    calin::iact_data::actl_event_decoder::ACTL_L0_CameraEventDecoder* decoder,
    const ZFITSDataSource_L0::config_type& config =
      ZFITSDataSource_L0::default_config());
  virtual ~ZFITSDataSourceOpener_L0();
  unsigned num_sources() override;
  std::string source_name(unsigned isource) override;
  ZFITSSingleFileDataSource_L0* open(unsigned isource) override;
private:
  calin::iact_data::zfits_actl_data_source::
    ZFITSACTL_L0_CameraEventDataSourceOpener* zfits_actl_opener_ = nullptr;
  calin::iact_data::actl_event_decoder::ACTL_L0_CameraEventDecoder* decoder_ = nullptr;
  ZFITSDataSource_L0::config_type config_;
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

class ZFITSSingleFileDataSource_R1:
  public calin::iact_data::telescope_data_source::
    TelescopeRandomAccessDataSourceWithRunConfig
{
public:
  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig);

  ZFITSSingleFileDataSource_R1(calin::iact_data::zfits_actl_data_source::
      ZFITSSingleFileACTL_R1_CameraEventDataSource* actl_zfits,
    bool dont_decode_run_configuration,
    calin::iact_data::actl_event_decoder::ACTL_R1_CameraEventDecoder* decoder,
    bool adopt_decoder = false,
    bool adopt_actl_zfits = false);

  ZFITSSingleFileDataSource_R1(const std::string& filename,
    calin::iact_data::actl_event_decoder::ACTL_R1_CameraEventDecoder* decoder,
    bool adopt_decoder = false,
    const config_type& config = default_config());

  virtual ~ZFITSSingleFileDataSource_R1();

  calin::ix::iact_data::telescope_event::TelescopeEvent* get_next(
    uint64_t& seq_index_out, google::protobuf::Arena** arena = nullptr) override;
  uint64_t size() override;
  void set_next_index(uint64_t next_index) override;

  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* get_run_configuration() override;

  static config_type default_config() {
    return zfits_actl_data_source::ZFITSSingleFileACTL_R1_CameraEventDataSource::default_config(); }

private:
  calin::iact_data::actl_event_decoder::ACTL_R1_CameraEventDecoder* decoder_;
  bool adopt_decoder_ = false;
  calin::iact_data::zfits_actl_data_source::
    ZFITSSingleFileACTL_R1_CameraEventDataSource* actl_zfits_ = nullptr;
  bool adopt_actl_zfits_ = false;
  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* run_config_ = nullptr;
};

class ZFITSDataSource_R1:
  public calin::io::data_source::BasicChainedRandomAccessDataSource<
    calin::iact_data::telescope_data_source::
      TelescopeRandomAccessDataSourceWithRunConfig>
{
public:
  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig);

  ZFITSDataSource_R1(const std::string& filename,
    calin::iact_data::actl_event_decoder::ACTL_R1_CameraEventDecoder* decoder,
    bool adopt_decoder = false,
    const config_type& config = default_config());
  virtual ~ZFITSDataSource_R1();

  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* get_run_configuration() override;

  static config_type default_config() {
    return zfits_actl_data_source::ZFITSACTL_R1_CameraEventDataSource::default_config(); }

protected:
  calin::iact_data::actl_event_decoder::ACTL_R1_CameraEventDecoder* decoder_;
  bool adopt_decoder_ = false;
  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* run_config_ = nullptr;
};

class ZFITSDataSourceOpener_R1:
  public calin::io::data_source::DataSourceOpener<
    calin::iact_data::telescope_data_source::
      TelescopeRandomAccessDataSourceWithRunConfig>
{
public:
  CALIN_TYPEALIAS(data_source_type, calin::iact_data::telescope_data_source::
    TelescopeRandomAccessDataSourceWithRunConfig);
  ZFITSDataSourceOpener_R1(std::string filename,
    calin::iact_data::actl_event_decoder::ACTL_R1_CameraEventDecoder* decoder,
    const ZFITSDataSource_R1::config_type& config =
      ZFITSDataSource_R1::default_config());
  virtual ~ZFITSDataSourceOpener_R1();
  unsigned num_sources() override;
  std::string source_name(unsigned isource) override;
  ZFITSSingleFileDataSource_R1* open(unsigned isource) override;
private:
  calin::iact_data::zfits_actl_data_source::
    ZFITSACTL_R1_CameraEventDataSourceOpener* zfits_actl_opener_ = nullptr;
  calin::iact_data::actl_event_decoder::ACTL_R1_CameraEventDecoder* decoder_ = nullptr;
  ZFITSDataSource_R1::config_type config_;
};

#endif
} } } // namespace calin::iact_data::zfits_data_source
