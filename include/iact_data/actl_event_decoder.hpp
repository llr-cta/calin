/*

   calin/iact_data/actl_event_decoder.hpp -- Stephen Fegan -- 2018-09-21

   A decoder of ACTL event types

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
#include <io/chained_data_source.hpp>
#include <iact_data/zfits_actl_data_source.hpp>
#include <iact_data/zfits_data_source.pb.h>
#include <iact_data/telescope_data_source.hpp>

#ifdef CALIN_HAVE_CTA_CAMERASTOACTL
#include <CoreMessages.pb.h>
#include <L0.pb.h>
#include <R1.pb.h>
#endif

namespace calin { namespace iact_data { namespace actl_event_decoder {

#ifdef CALIN_HAVE_CTA_CAMERASTOACTL

void decode_cdts_data(calin::ix::iact_data::telescope_event::CDTSData* calin_cdts_data,
  const DataModel::AnyArray& cta_array);

void decode_tib_data(calin::ix::iact_data::telescope_event::TIBData* calin_tib_data,
  const DataModel::AnyArray& cta_array);

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

class ACTL_L0_CameraEventDecoder
{
public:
  virtual ~ACTL_L0_CameraEventDecoder();
  virtual bool decode(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    const DataModel::CameraEvent* cta_event) = 0;
  virtual bool decode_run_config(
    calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config,
    const DataModel::CameraRunHeader* cta_run_header,
    const DataModel::CameraEvent* cta_event) = 0;
};

class DecodedACTL_L0_CameraEventDataSource:
  public calin::iact_data::telescope_data_source::TelescopeDataSource
{
public:
  DecodedACTL_L0_CameraEventDataSource(
    calin::iact_data::zfits_actl_data_source::ACTL_L0_CameraEventDataSource* actl_src,
    ACTL_L0_CameraEventDecoder* decoder,
    bool adopt_actl_src = false, bool adopt_decoder = false);

  virtual ~DecodedACTL_L0_CameraEventDataSource();

  calin::ix::iact_data::telescope_event::TelescopeEvent* get_next(
    uint64_t& seq_index_out, google::protobuf::Arena** arena = nullptr) override;

private:
  ACTL_L0_CameraEventDecoder* decoder_;
  bool adopt_decoder_ = false;
  calin::iact_data::zfits_actl_data_source::ACTL_L0_CameraEventDataSource* actl_src_ = nullptr;
  bool adopt_actl_src_ = false;
};

class DecodedConstACTL_L0_CameraEventDataSource:
  public calin::iact_data::telescope_data_source::TelescopeDataSource
{
public:
  DecodedConstACTL_L0_CameraEventDataSource(
    calin::iact_data::zfits_actl_data_source::ConstACTL_L0_CameraEventDataSource* actl_src,
    calin::iact_data::zfits_actl_data_source::ConstACTL_L0_CameraEventDataSink* actl_sink,
    ACTL_L0_CameraEventDecoder* decoder,
    bool adopt_actl_src = false, bool adopt_actl_sink = false, bool adopt_decoder = false);

  DecodedConstACTL_L0_CameraEventDataSource(
    calin::iact_data::zfits_actl_data_source::ConstACTL_L0_CameraEventDataSource* actl_src,
    ACTL_L0_CameraEventDecoder* decoder,
    bool adopt_actl_src = false, bool adopt_decoder = false);

  virtual ~DecodedConstACTL_L0_CameraEventDataSource();

  calin::ix::iact_data::telescope_event::TelescopeEvent* get_next(
    uint64_t& seq_index_out, google::protobuf::Arena** arena = nullptr) override;

private:
  ACTL_L0_CameraEventDecoder* decoder_;
  bool adopt_decoder_ = false;
  calin::iact_data::zfits_actl_data_source::ConstACTL_L0_CameraEventDataSource* actl_src_ = nullptr;
  bool adopt_actl_src_ = false;
  calin::iact_data::zfits_actl_data_source::ConstACTL_L0_CameraEventDataSink* actl_sink_ = nullptr;
  bool adopt_actl_sink_ = false;
};

class DecodedConstACTL_L0_CameraEventDataSourceFactory:
  public calin::iact_data::telescope_data_source::TelescopeDataSourceFactory
{
public:
  DecodedConstACTL_L0_CameraEventDataSourceFactory(
    calin::io::data_source::BidirectionalBufferedDataSourcePump<
      const DataModel::CameraEvent>* pump, ACTL_L0_CameraEventDecoder* decoder,
    bool adopt_pump = false, bool adopt_decoder = false);
  virtual ~DecodedConstACTL_L0_CameraEventDataSourceFactory();
  DecodedConstACTL_L0_CameraEventDataSource* new_data_source() override;

private:
  ACTL_L0_CameraEventDecoder* decoder_;
  bool adopt_decoder_ = false;
  calin::io::data_source::BidirectionalBufferedDataSourcePump<
    const DataModel::CameraEvent>* pump_ = nullptr;
  bool adopt_pump_ = false;
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

class ACTL_R1_CameraEventDecoder
{
public:
  virtual ~ACTL_R1_CameraEventDecoder();
  virtual bool decode(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    const R1::CameraEvent* cta_event) = 0;
  virtual bool decode_run_config(
    calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config,
    const R1::CameraConfiguration* cta_run_header,
    const R1::CameraEvent* cta_event) = 0;
};

class DecodedACTL_R1_CameraEventDataSource:
  public calin::iact_data::telescope_data_source::TelescopeDataSource
{
public:
  DecodedACTL_R1_CameraEventDataSource(
    calin::iact_data::zfits_actl_data_source::ACTL_R1_CameraEventDataSource* actl_src,
    ACTL_R1_CameraEventDecoder* decoder,
    bool adopt_actl_src = false, bool adopt_decoder = false);

  virtual ~DecodedACTL_R1_CameraEventDataSource();

  calin::ix::iact_data::telescope_event::TelescopeEvent* get_next(
    uint64_t& seq_index_out, google::protobuf::Arena** arena = nullptr) override;

protected:
  ACTL_R1_CameraEventDecoder* decoder_;
  virtual R1::CameraEvent* do_get_next(uint64_t& seq_index_out);

private:
  bool adopt_decoder_ = false;
  calin::iact_data::zfits_actl_data_source::ACTL_R1_CameraEventDataSource* actl_src_ = nullptr;
  bool adopt_actl_src_ = false;
};

class DecodedACTL_R1_CameraEventDataSourceWithRunConfig:
  public virtual calin::iact_data::telescope_data_source::TelescopeDataSourceWithRunConfig,
  public virtual DecodedACTL_R1_CameraEventDataSource
{
public:
  DecodedACTL_R1_CameraEventDataSourceWithRunConfig(
    calin::iact_data::zfits_actl_data_source::ACTL_R1_CameraEventDataSourceWithRunHeader* actl_src,
    ACTL_R1_CameraEventDecoder* decoder,
    bool adopt_actl_src = false, bool adopt_decoder = false);

  virtual ~DecodedACTL_R1_CameraEventDataSourceWithRunConfig();

  calin::ix::iact_data::telescope_event::TelescopeEvent* get_next(
    uint64_t& seq_index_out, google::protobuf::Arena** arena = nullptr) override;

  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* get_run_configuration() override;

protected:
  R1::CameraEvent* do_get_next(uint64_t& seq_index_out) override;
  void ensure_run_config();

private:
  calin::iact_data::zfits_actl_data_source::ACTL_R1_CameraEventDataSourceWithRunHeader* actl_src_ = nullptr;
  R1::CameraEvent* saved_event_ = nullptr;
  uint64_t saved_seq_index_ = 0;
  calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config_ = nullptr;
};

class DecodedConstACTL_R1_CameraEventDataSource:
  public calin::iact_data::telescope_data_source::TelescopeDataSource
{
public:
  DecodedConstACTL_R1_CameraEventDataSource(
    calin::iact_data::zfits_actl_data_source::ConstACTL_R1_CameraEventDataSource* actl_src,
    calin::iact_data::zfits_actl_data_source::ConstACTL_R1_CameraEventDataSink* actl_sink,
    ACTL_R1_CameraEventDecoder* decoder,
    bool adopt_actl_src = false, bool adopt_actl_sink = false, bool adopt_decoder = false);

  DecodedConstACTL_R1_CameraEventDataSource(
    calin::iact_data::zfits_actl_data_source::ConstACTL_R1_CameraEventDataSource* actl_src,
    ACTL_R1_CameraEventDecoder* decoder,
    bool adopt_actl_src = false, bool adopt_decoder = false);

  virtual ~DecodedConstACTL_R1_CameraEventDataSource();

  calin::ix::iact_data::telescope_event::TelescopeEvent* get_next(
    uint64_t& seq_index_out, google::protobuf::Arena** arena = nullptr) override;

private:
  ACTL_R1_CameraEventDecoder* decoder_;
  bool adopt_decoder_ = false;
  calin::iact_data::zfits_actl_data_source::ConstACTL_R1_CameraEventDataSource* actl_src_ = nullptr;
  bool adopt_actl_src_ = false;
  calin::iact_data::zfits_actl_data_source::ConstACTL_R1_CameraEventDataSink* actl_sink_ = nullptr;
  bool adopt_actl_sink_ = false;
};

class DecodedConstACTL_R1_CameraEventDataSourceFactory:
  public calin::iact_data::telescope_data_source::TelescopeDataSourceFactory
{
public:
  DecodedConstACTL_R1_CameraEventDataSourceFactory(
    calin::io::data_source::BidirectionalBufferedDataSourcePump<
      const R1::CameraEvent>* pump, ACTL_R1_CameraEventDecoder* decoder,
    bool adopt_pump = false, bool adopt_decoder = false);
  virtual ~DecodedConstACTL_R1_CameraEventDataSourceFactory();
  DecodedConstACTL_R1_CameraEventDataSource* new_data_source() override;

private:
  ACTL_R1_CameraEventDecoder* decoder_;
  bool adopt_decoder_ = false;
  calin::io::data_source::BidirectionalBufferedDataSourcePump<
    const R1::CameraEvent>* pump_ = nullptr;
  bool adopt_pump_ = false;
};

#endif
} } } // namespace calin::iact_data::zfits_data_source
