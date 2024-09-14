/*

   calin/iact_data/zfits_acada_data_source.hpp -- Stephen Fegan -- 2022-09-05

   Source of "raw" ACADA data types from ZFITS files

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

#include <iact_data/acada_data_source.hpp>
#include <iact_data/zfits_data_source.pb.h>

#include <ProtobufIFits.h>

namespace calin { namespace iact_data { namespace zfits_acada_data_source {

std::vector<std::string> get_zfits_table_names(std::string filename);
std::vector<std::string> get_zfits_table_column_names(std::string filename, std::string tablename);
std::vector<std::string> get_zfits_table_keys(std::string filename, std::string tablename);
std::map<std::string,std::string> get_zfits_table_key_values(std::string filename, std::string tablename);

template<typename Message>
class ZFITSSingleFileSingleMessageDataSource:
  public virtual calin::io::data_source::RandomAccessDataSource<Message>
{
public:
  CALIN_TYPEALIAS(message_type, Message);

  ZFITSSingleFileSingleMessageDataSource(const std::string& filename, const std::string& tablename = {},
    bool suppress_file_record = false);
  virtual ~ZFITSSingleFileSingleMessageDataSource();

  const message_type* borrow_next_message(uint64_t& seq_index_out);
  void release_borrowed_message(const message_type* message);

  message_type* get_next(uint64_t& seq_index_out,
    google::protobuf::Arena** arena = nullptr) override;

  uint64_t size() override;
  void set_next_index(uint64_t next_index) override;
  uint64_t get_next_index() const { return next_message_index_; };

private:
  std::string filename_;
  std::string tablename_;
  ADH::IO::ProtobufIFits* zfits_ = nullptr;
  calin::ix::provenance::chronicle::FileIORecord* file_record_ = nullptr;
  uint64_t next_message_index_ = 0;
};

template<typename MessageSet>
class ZFITSSingleFileACADACameraEventDataSource:
  public calin::iact_data::acada_data_source::
    ACADACameraEventRandomAccessDataSourceWithRunHeader<MessageSet>
{
public:
  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig);
  CALIN_TYPEALIAS(message_set_type, MessageSet);
  CALIN_TYPEALIAS(event_type, typename MessageSet::event_type);
  CALIN_TYPEALIAS(header_type, typename MessageSet::header_type);
  CALIN_TYPEALIAS(data_stream_type, typename MessageSet::data_stream_type);

  ZFITSSingleFileACADACameraEventDataSource(const std::string& filename,
    const config_type& config = default_config());
  virtual ~ZFITSSingleFileACADACameraEventDataSource();

  const event_type* borrow_next_event(uint64_t& seq_index_out) override;
  void release_borrowed_event(const event_type* event) override;

  event_type* get_next(uint64_t& seq_index_out,
    google::protobuf::Arena** arena = nullptr) override;
  uint64_t size() override;
  void set_next_index(uint64_t next_index) override;

  header_type* get_run_header() override;
  data_stream_type* get_data_stream() override;

  static config_type default_config();
  const config_type& config() const { return config_; }

private:
  std::string filename_;
  ZFITSSingleFileSingleMessageDataSource<event_type>* zfits_;
  header_type* run_header_ = nullptr;
  data_stream_type* data_stream_ = nullptr;
  config_type config_;
};

template<typename MessageSet>
class ZFITSACADACameraEventDataSource:
  public calin::io::data_source::BasicChainedRandomAccessDataSource<
    calin::iact_data::acada_data_source::
      ACADACameraEventRandomAccessDataSourceWithRunHeader<MessageSet> >
{
public:
  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig);
  CALIN_TYPEALIAS(message_set_type, MessageSet);
  CALIN_TYPEALIAS(event_type, typename MessageSet::event_type);
  CALIN_TYPEALIAS(header_type, typename MessageSet::header_type);
  CALIN_TYPEALIAS(data_stream_type, typename MessageSet::data_stream_type);
  CALIN_TYPEALIAS(BaseDataSource, calin::io::data_source::BasicChainedRandomAccessDataSource<
    calin::iact_data::acada_data_source::
      ACADACameraEventRandomAccessDataSourceWithRunHeader<MessageSet> >);

  ZFITSACADACameraEventDataSource(const std::string& filename,
    const config_type& config = default_config());
  virtual ~ZFITSACADACameraEventDataSource();

  header_type* get_run_header() override;
  data_stream_type* get_data_stream() override;

  const event_type* borrow_next_event(uint64_t& seq_index_out) override;
  void release_borrowed_event(const event_type* event) override;

  event_type* get_next(uint64_t& seq_index_out,
    google::protobuf::Arena** arena = nullptr) override;
  uint64_t size() override;
  void set_next_index(uint64_t next_index) override;

  static config_type default_config();
  const config_type& config() const { return config_; }

protected:
  using BaseDataSource::source_;
  using BaseDataSource::opener_;
  using BaseDataSource::isource_;
  using BaseDataSource::seq_index_;
  using BaseDataSource::open_file;

private:
  config_type config_;
  header_type* run_header_ = nullptr;
  data_stream_type* data_stream_ = nullptr;
};

template<typename MessageSet>
class ZFITSACADACameraEventDataSourceOpener:
  public calin::io::data_source::DataSourceOpener<
    calin::iact_data::acada_data_source::
      ACADACameraEventRandomAccessDataSourceWithRunHeader<MessageSet> >
{
public:
  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig);
  CALIN_TYPEALIAS(message_set_type, MessageSet);
  CALIN_TYPEALIAS(event_type, typename MessageSet::event_type);
  CALIN_TYPEALIAS(header_type, typename MessageSet::header_type);
  CALIN_TYPEALIAS(data_stream_type, typename MessageSet::data_stream_type);

  ZFITSACADACameraEventDataSourceOpener(std::string filename,
    const config_type& config = default_config());
  virtual ~ZFITSACADACameraEventDataSourceOpener();
  unsigned num_sources() const override;
  std::string source_name(unsigned isource) const override;
  ZFITSSingleFileACADACameraEventDataSource<MessageSet>* open(unsigned isource) override;
  bool has_opened_file() { return has_opened_file_; }

  static config_type default_config();
  const config_type& config() const { return config_; }
private:
  std::vector<std::string> filenames_;
  config_type config_;
  bool has_opened_file_ = false;
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

// bool is_zfits_l0(std::string filename, std::string events_table_name = "");

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

// bool is_zfits_r1v0(std::string filename, std::string events_table_name = "");

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


} } } // namespace calin::iact_data::zfits_acada_data_source

