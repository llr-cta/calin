/*

   calin/iact_data/zfits_data_source.hpp -- Stephen Fegan -- 2016-01-21

   A supplier of single telescope data from CTA ACTL ZFits data files

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

#include <iact_data/acada_data_source.hpp>
#include <iact_data/acada_event_decoder.hpp>
#include <iact_data/zfits_acada_data_source.hpp>

#include <iact_data/zfits_data_source.pb.h>
#include <iact_data/telescope_data_source.hpp>

#include <CoreMessages.pb.h>

namespace calin { namespace iact_data { namespace zfits_data_source {

#if 0 // UNUSED
template<typename EventMessage, typename HeaderMessage>
class ZFITSSingleFileDataSource:
  public calin::iact_data::telescope_data_source::
    TelescopeRandomAccessDataSourceWithRunConfig
{
public:
  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig);

  ZFITSSingleFileDataSource(
    calin::iact_data::zfits_acada_data_source::
      ZFITSSingleFileACADACameraEventDataSource<EventMessage,HeaderMessage>* actl_zfits,
    bool dont_decode_run_configuration,
    calin::iact_data::acada_event_decoder::
      ACADACameraEventDecoder<EventMessage,HeaderMessage>* decoder,
    bool adopt_decoder = false,
    bool adopt_actl_zfits = false);

  ZFITSSingleFileDataSource(const std::string& filename,
    calin::iact_data::acada_event_decoder::
      ACADACameraEventDecoder<EventMessage,HeaderMessage>* decoder,
    bool adopt_decoder = false,
    const config_type& config = default_config());

  virtual ~ZFITSSingleFileDataSource();

  calin::ix::iact_data::telescope_event::TelescopeEvent* get_next(
    uint64_t& seq_index_out, google::protobuf::Arena** arena = nullptr) override;
  uint64_t size() override;
  void set_next_index(uint64_t next_index) override;

  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* get_run_configuration() override;

  static config_type default_config();

private:
  calin::iact_data::acada_event_decoder::
    ACADACameraEventDecoder<EventMessage,HeaderMessage>* decoder_;
  bool adopt_decoder_ = false;
  calin::iact_data::zfits_acada_data_source::
    ZFITSSingleFileACADACameraEventDataSource<EventMessage,HeaderMessage>* acada_zfits_ = nullptr;
  bool adopt_acada_zfits_ = false;
  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* run_config_ = nullptr;
};
#endif // UNUSED


template<typename MessageSet>
class ZFITSDataSource:
  public calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSourceWithRunConfig,
  public calin::io::data_source::FragmentList
{
public:
  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig);

  ZFITSDataSource(const std::string& filename,
    calin::iact_data::acada_event_decoder::
      ACADACameraEventDecoder<MessageSet>* decoder,
    bool adopt_decoder = false,
    const config_type& config = default_config());

  ZFITSDataSource(const std::string& filename,
    ZFITSDataSource<MessageSet>* base_datasource, 
    const config_type& config = default_config());

  virtual ~ZFITSDataSource();

  calin::ix::iact_data::telescope_event::TelescopeEvent* get_next(
    uint64_t& seq_index_out, google::protobuf::Arena** arena = nullptr) override;

  uint64_t size() override;
  void set_next_index(uint64_t next_index) override;

  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* get_run_configuration() override;

  unsigned current_fragment_index() const override;
  unsigned num_fragments() const override;
  std::string fragment_name(unsigned index) const override;

  static config_type default_config();

protected:

  calin::iact_data::acada_event_decoder::
    ACADACameraEventDecoder<MessageSet>* decoder_;
  bool adopt_decoder_ = false;
  calin::iact_data::zfits_acada_data_source::
    ZFITSACADACameraEventDataSource<MessageSet>* acada_zfits_ = nullptr;
  bool adopt_acada_zfits_ = false;
  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* run_config_ = nullptr;
};

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

CALIN_TYPEALIAS(ZFITSDataSource_L0,
  ZFITSDataSource<calin::iact_data::acada_data_source::ACADA_MessageSet_L0>);

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

CALIN_TYPEALIAS(ZFITSDataSource_R1v0,
  ZFITSDataSource<calin::iact_data::acada_data_source::ACADA_MessageSet_R1v0>);

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

CALIN_TYPEALIAS(ZFITSDataSource_R1v1,
  ZFITSDataSource<calin::iact_data::acada_data_source::ACADA_MessageSet_R1v1>);

} } } // namespace calin::iact_data::zfits_data_source
