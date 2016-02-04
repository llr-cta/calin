/*

   calin/iact_data/zfits_data_source.hpp -- Stephen Fegan -- 2016-01-21

   A supplier of single telescope data from CTA ACTL ZFits data files

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

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
#include <iact_data/zfits_data_source.pb.h>
#include <iact_data/telescope_data_source.hpp>

// Forward declaration of ACTL::IO::ProtobufIFits

namespace ACTL { namespace IO {
  class ProtobufIFits;
} } // namespace ACTL::IO
namespace DataModel {
  class CameraEvent;
} // namespace DataModel

namespace calin { namespace iact_data { namespace zfits_data_source {

class CTACameraEventDecoder
{
public:
  virtual ~CTACameraEventDecoder();
  virtual calin::ix::iact_data::telescope_event::TelescopeEvent*
    decode(const DataModel::CameraEvent* cta_event) = 0;
};

class ZFITSSingleFileDataSource: public
  calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSource
{
public:
  ZFITSSingleFileDataSource(const std::string& filename,
    CTACameraEventDecoder* decoder, bool adopt_decoder = false);
  virtual ~ZFITSSingleFileDataSource();

  calin::ix::iact_data::telescope_event::TelescopeEvent* get_next() override;
  uint64_t size() override;
  void set_next_index(uint64_t next_index) override;

private:
  std::string filename_;
  ACTL::IO::ProtobufIFits* zfits_ = nullptr;
  CTACameraEventDecoder* decoder_ = nullptr;
  bool adopt_decoder_ = false;
  uint64_t next_event_index_ = 0;
};

class ZFitsDataSourceOpener:
  public calin::io::data_source::DataSourceOpener<
    calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSource>
{
public:
  CALIN_TYPEALIAS(data_source_type,
    calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSource);
  ZFitsDataSourceOpener(std::string filename,
    CTACameraEventDecoder* decoder, bool adopt_decoder = false,
    const std::string& extension = ".fits.fz");
  virtual ~ZFitsDataSourceOpener();
  unsigned num_sources() override;
  calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSource*
    open(unsigned isource) override;
private:
  std::vector<std::string> filenames_;
  CTACameraEventDecoder* decoder_ = nullptr;
  bool adopt_decoder_ = false;
};

class ZFITSDataSource:
  public calin::io::data_source::BasicChaninedRandomAccessDataSource<
    calin::iact_data::telescope_data_source::TelescopeRandomAccessDataSource>
{
public:
  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig);

  ZFITSDataSource(const std::string& filename,
    CTACameraEventDecoder* decoder, bool adopt_decoder = false,
    const config_type& config = default_config());
  virtual ~ZFITSDataSource();

  static const config_type& default_config() { return default_config_.x; }
  static config_type* new_default_config() {
    return new config_type(default_config()); }
  const config_type& config() const { return config_; }

private:
#ifndef SWIG
  // Purpose is to set the extension to our default value in the static member
  class config_helper
  {
  public:
    config_helper() { x.set_extension(".fits.fz"); }
    config_type x;
  };
  static config_helper default_config_;
#endif
  config_type config_;
};

} } } // namespace calin::iact_data::nectarcam_data_source
