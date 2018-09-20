/*

   calin/iact_data/zfits_acts_data_source.hpp -- Stephen Fegan -- 2016-05-04

   A supplier of ACTL telescope data types from CTA ACTL ZFits data files

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
#include <io/data_source.hpp>
#include <io/chained_data_source.hpp>
#include <iact_data/zfits_data_source.pb.h>

#ifdef CALIN_HAVE_CTA_CAMERASTOACTL

#include <ProtobufIFits.h>
#include <L0.pb.h>
#include <R1.pb.h>

namespace calin { namespace iact_data { namespace zfits_actl_data_source {

CALIN_TYPEALIAS(ACTL_L0_CameraEventDataSource,
  calin::io::data_source::DataSource<DataModel::CameraEvent>);
CALIN_TYPEALIAS(ACTL_L0_CameraEventRandomAccessDataSource,
  calin::io::data_source::RandomAccessDataSource<DataModel::CameraEvent>);
CALIN_TYPEALIAS(ConstACTL_L0_CameraEventDataSource,
  calin::io::data_source::DataSource<const DataModel::CameraEvent>);
CALIN_TYPEALIAS(ConstACTL_L0_CameraEventDataSink,
  calin::io::data_source::DataSink<const DataModel::CameraEvent>);

CALIN_TYPEALIAS(ACTL_R1_CameraEventDataSource,
  calin::io::data_source::DataSource<R1::CameraEvent>);
CALIN_TYPEALIAS(ACTL_R1_CameraEventRandomAccessDataSource,
  calin::io::data_source::RandomAccessDataSource<R1::CameraEvent>);
CALIN_TYPEALIAS(ConstACTL_R1_CameraEventDataSource,
  calin::io::data_source::DataSource<const R1::CameraEvent>);
CALIN_TYPEALIAS(ConstACTL_R1_CameraEventDataSink,
  calin::io::data_source::DataSink<const R1::CameraEvent>);

} } } // namespace calin::iact_data::zfits_actl_data_source

#ifdef SWIG
%template(ACTL_L0_CameraEventRandomAccessDataSource)
  calin::io::data_source::RandomAccessDataSource<DataModel::CameraEvent>;
%template(ACTL_R1_CameraEventRandomAccessDataSource)
  calin::io::data_source::RandomAccessDataSource<R1::CameraEvent>;
#endif

namespace calin { namespace iact_data { namespace zfits_actl_data_source {

class ACTL_L0_CameraEventRandomAccessDataSourceWithRunHeader:
  public calin::io::data_source::RandomAccessDataSource<DataModel::CameraEvent>
{
public:
  ACTL_L0_CameraEventRandomAccessDataSourceWithRunHeader():
      calin::io::data_source::RandomAccessDataSource<DataModel::CameraEvent>() {
    /* nothing to see here */ }
  virtual ~ACTL_L0_CameraEventRandomAccessDataSourceWithRunHeader();
  virtual DataModel::CameraRunHeader* get_run_header() = 0;

  virtual const DataModel::CameraEvent* borrow_next_event(uint64_t& seq_index_out) = 0;
  virtual void release_borrowed_event(const DataModel::CameraEvent* event) = 0;
};

class ACTL_R1_CameraEventRandomAccessDataSourceWithRunHeader:
  public calin::io::data_source::RandomAccessDataSource<R1::CameraEvent>
{
public:
  ACTL_R1_CameraEventRandomAccessDataSourceWithRunHeader():
      calin::io::data_source::RandomAccessDataSource<R1::CameraEvent>() {
    /* nothing to see here */ }
  virtual ~ACTL_R1_CameraEventRandomAccessDataSourceWithRunHeader();
  virtual R1::CameraConfiguration* get_run_header() = 0;

  virtual const R1::CameraEvent* borrow_next_event(uint64_t& seq_index_out) = 0;
  virtual void release_borrowed_event(const R1::CameraEvent* event) = 0;
};

CALIN_TYPEALIAS(ACTL_L0_CameraEventChainedRandomAccessDataSourceWithRunHeader,
  calin::io::data_source::BasicChainedRandomAccessDataSource<
    ACTL_L0_CameraEventRandomAccessDataSourceWithRunHeader>);

CALIN_TYPEALIAS(ACTL_L0_CameraEventRandomAccessDataSourceOpener,
  calin::io::data_source::DataSourceOpener<
    ACTL_L0_CameraEventRandomAccessDataSourceWithRunHeader>);

CALIN_TYPEALIAS(ACTL_R1_CameraEventChainedRandomAccessDataSourceWithRunHeader,
  calin::io::data_source::BasicChainedRandomAccessDataSource<
    ACTL_R1_CameraEventRandomAccessDataSourceWithRunHeader>);

CALIN_TYPEALIAS(ACTL_R1_CameraEventRandomAccessDataSourceOpener,
  calin::io::data_source::DataSourceOpener<
    ACTL_R1_CameraEventRandomAccessDataSourceWithRunHeader>);

} } } // namespace calin::iact_data::zfits_actl_data_source

#ifdef SWIG
%template(ACTL_L0_CameraEventChainedDataSourceWithRunHeader)
  calin::io::data_source::BasicChainedDataSource<
    calin::iact_data::zfits_actl_data_source::
      ACTL_L0_CameraEventRandomAccessDataSourceWithRunHeader>;
%template(ACTL_L0_CameraEventChainedRandomAccessDataSourceWithRunHeader)
  calin::io::data_source::BasicChainedRandomAccessDataSource<
    calin::iact_data::zfits_actl_data_source::
      ACTL_L0_CameraEventRandomAccessDataSourceWithRunHeader>;
%template(ACTL_L0_CameraEventRandomAccessDataSourceOpener)
  calin::io::data_source::DataSourceOpener<
    calin::iact_data::zfits_actl_data_source::
      ACTL_L0_CameraEventRandomAccessDataSourceWithRunHeader>;

%template(ACTL_R1_CameraEventChainedDataSourceWithRunHeader)
  calin::io::data_source::BasicChainedDataSource<
    calin::iact_data::zfits_actl_data_source::
      ACTL_R1_CameraEventRandomAccessDataSourceWithRunHeader>;
%template(ACTL_R1_CameraEventChainedRandomAccessDataSourceWithRunHeader)
  calin::io::data_source::BasicChainedRandomAccessDataSource<
    calin::iact_data::zfits_actl_data_source::
      ACTL_R1_CameraEventRandomAccessDataSourceWithRunHeader>;
%template(ACTL_R1_CameraEventRandomAccessDataSourceOpener)
  calin::io::data_source::DataSourceOpener<
    calin::iact_data::zfits_actl_data_source::
      ACTL_R1_CameraEventRandomAccessDataSourceWithRunHeader>;
#endif

namespace calin { namespace iact_data { namespace zfits_actl_data_source {

class ZFITSSingleFileACTL_L0_CameraEventDataSource:
  public ACTL_L0_CameraEventRandomAccessDataSourceWithRunHeader
{
public:
  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig);

  ZFITSSingleFileACTL_L0_CameraEventDataSource(const std::string& filename,
    config_type config = default_config());
  virtual ~ZFITSSingleFileACTL_L0_CameraEventDataSource();

  const DataModel::CameraEvent* borrow_next_event(uint64_t& seq_index_out) override;
  void release_borrowed_event(const DataModel::CameraEvent* event) override;

  DataModel::CameraEvent* get_next(uint64_t& seq_index_out,
    google::protobuf::Arena** arena = nullptr) override;
  uint64_t size() override;
  void set_next_index(uint64_t next_index) override;

  DataModel::CameraRunHeader* get_run_header() override;

  static config_type default_config();
  const config_type& config() const { return config_; }

private:
  std::string filename_;
  ACTL::IO::ProtobufIFits* zfits_ = nullptr;
  uint64_t next_event_index_ = 0;
  DataModel::CameraRunHeader* run_header_ = nullptr;
  config_type config_;
};

class ZFITSACTL_L0_CameraEventDataSource:
  public calin::io::data_source::BasicChainedRandomAccessDataSource<
    ACTL_L0_CameraEventRandomAccessDataSourceWithRunHeader>
{
public:
  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig);

  ZFITSACTL_L0_CameraEventDataSource(const std::string& filename,
    const config_type& config = default_config());
  virtual ~ZFITSACTL_L0_CameraEventDataSource();

  DataModel::CameraRunHeader* get_run_header() override;

  static config_type default_config() {
    return ZFITSSingleFileACTL_L0_CameraEventDataSource::default_config(); }
  const config_type& config() const { return config_; }

  const DataModel::CameraEvent* borrow_next_event(uint64_t& seq_index_out) override;
  void release_borrowed_event(const DataModel::CameraEvent* event) override;

  DataModel::CameraEvent* get_next(uint64_t& seq_index_out,
    google::protobuf::Arena** arena = nullptr) override;
  uint64_t size() override;
  void set_next_index(uint64_t next_index) override;

private:
  config_type config_;
  DataModel::CameraRunHeader* run_header_ = nullptr;
};

class ZFITSACTL_L0_CameraEventDataSourceOpener:
  public calin::io::data_source::DataSourceOpener<
    ACTL_L0_CameraEventRandomAccessDataSourceWithRunHeader>
{
public:
  CALIN_TYPEALIAS(data_source_type, ACTL_L0_CameraEventRandomAccessDataSourceWithRunHeader);
  ZFITSACTL_L0_CameraEventDataSourceOpener(std::string filename,
    const ZFITSACTL_L0_CameraEventDataSource::config_type& config =
      ZFITSACTL_L0_CameraEventDataSource::default_config());
  virtual ~ZFITSACTL_L0_CameraEventDataSourceOpener();
  unsigned num_sources() override;
  std::string source_name(unsigned isource) override;
  ZFITSSingleFileACTL_L0_CameraEventDataSource* open(unsigned isource) override;
  bool has_opened_file() { return has_opened_file_; }
private:
  std::vector<std::string> filenames_;
  ZFITSACTL_L0_CameraEventDataSource::config_type config_;
  bool has_opened_file_ = false;
};

class ZFITSConstACTL_L0_CameraEventDataSourceBorrowAdapter:
  public calin::io::data_source::DataSource<const DataModel::CameraEvent>
{
public:
  ZFITSConstACTL_L0_CameraEventDataSourceBorrowAdapter(ZFITSACTL_L0_CameraEventDataSource* src);
  virtual ~ZFITSConstACTL_L0_CameraEventDataSourceBorrowAdapter();
  const DataModel::CameraEvent* get_next(uint64_t& seq_index_out,
    google::protobuf::Arena** arena = nullptr) override;
private:
  ZFITSACTL_L0_CameraEventDataSource* src_;
};

class ZFITSConstACTL_L0_CameraEventDataSourceReleaseAdapter:
  public calin::io::data_source::DataSink<const DataModel::CameraEvent>
{
public:
  ZFITSConstACTL_L0_CameraEventDataSourceReleaseAdapter(ZFITSACTL_L0_CameraEventDataSource* src);
  virtual ~ZFITSConstACTL_L0_CameraEventDataSourceReleaseAdapter();
  bool put_next(const DataModel::CameraEvent* data, uint64_t seq_index,
    google::protobuf::Arena* arena = nullptr, bool adopt_data = false) override;
private:
  ZFITSACTL_L0_CameraEventDataSource* src_;
};

} } } // namespace calin::iact_data::zfits_actl_data_source

#endif
