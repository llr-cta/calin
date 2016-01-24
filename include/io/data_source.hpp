/*

   calin/io/data_source.hpp -- Stephen Fegan -- 2016-01-15

   A supplier of generic data, probably Protobufs

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

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>

#include <io/packet_stream.hpp>
#include <io/packet_stream.pb.h>

namespace calin { namespace io { namespace data_source {

template<typename T> class DataSource
{
public:
  CALIN_TYPEALIAS(data_type, T);
  virtual ~DataSource() { }
  virtual T* get_next() = 0;
};

template<typename T> class RandomAccessDataSource: public DataSource<T>
{
public:
  CALIN_TYPEALIAS(data_type, T);
  virtual ~RandomAccessDataSource() { }
  virtual T* get_next() = 0;
  virtual uint64_t size() = 0;
  virtual void set_next_index(uint64_t next_index) = 0;
};

template<typename T> class DataSink
{
public:
  CALIN_TYPEALIAS(data_type, T);
  virtual ~DataSink() { }
  virtual bool put_next(T* data, bool adopt_data = false) = 0;
};

// *****************************************************************************
//
// ProtobufPacketStreamDataSource/Sink : transform a PacketIn/OutStream into a
// DataSource/Sink for protobufs using protobuf (de)serealisation
//
// *****************************************************************************

template<typename T> class ProtobufPacketStreamDataSource: public DataSource<T>
{
public:
  ProtobufPacketStreamDataSource(
    calin::io::packet_stream::PacketInStream* stream, bool adopt_stream=false):
    DataSource<T>(), stream_(stream), adopt_stream_(adopt_stream)
    { /* nothing to see here */ }
  virtual ~ProtobufPacketStreamDataSource() {
    if(adopt_stream_)delete(stream_); }
  T* get_next() override {
    std::string serialized_d;
    if(!stream_->get_packet(serialized_d))return nullptr;
    T* d = new T;
    if(!d->ParseFromString(serialized_d)){ delete d; return nullptr; }
    return d;
  }
private:
  calin::io::packet_stream::PacketInStream* stream_ = nullptr;
  bool adopt_stream_;
};

template<typename T> class ProtobufPacketStreamDataSink: public DataSink<T>
{
public:
  ProtobufPacketStreamDataSink(
    calin::io::packet_stream::PacketOutStream* stream, bool adopt_stream=false):
    DataSink<T>(), stream_(stream), adopt_stream_(adopt_stream)
    { /* nothing to see here */ }
  virtual ~ProtobufPacketStreamDataSink() {
    if(adopt_stream_)delete(stream_); }
  bool put_next(T* data, bool adopt_data = false) override {
    std::unique_ptr<T> data_guard(adopt_data ? data : nullptr);
    return stream_->put_packet(data->SerializeAsString()); }
private:
  calin::io::packet_stream::PacketOutStream* stream_ = nullptr;
  bool adopt_stream_;
};

// *****************************************************************************
//
// ProtobufFileDataSource/Sink : specialize ProtobufPacketStreamDataSource/Sink
// using FramedFilePacketIn/OutStream, mainly for convenience of construction
//
// *****************************************************************************

template<typename T> class ProtobufFileDataSource: public DataSource<T>
{
public:
  CALIN_TYPEALIAS(config_type,
    ix::io::packet_stream::CompressedPacketInStreamOptions);
  ProtobufFileDataSource(const std::string& filename,
    const config_type& config = config_type::default_instance()):
    DataSource<T>()
  {
    auto* stream =
      new packet_stream::FramedFilePacketInStream(filename, config);
    source_ = new ProtobufPacketStreamDataSource<T>(stream, true);
  }
  virtual ~ProtobufFileDataSource() { delete source_; }
  T* get_next() override { return source_->get_next(); }
  static config_type default_options() {
    return config_type::default_instance(); }
private:
  ProtobufPacketStreamDataSource<T>* source_ = nullptr;
};

template<typename T> class ProtobufFileDataSink: public DataSink<T>
{
public:
  CALIN_TYPEALIAS(config_type,
    ix::io::packet_stream::CompressedPacketOutStreamOptions);
  ProtobufFileDataSink(const std::string& filename, bool append = false,
    const config_type& config = config_type::default_instance()):
    DataSink<T>()
  {
    auto* stream =
      new packet_stream::FramedFilePacketOutStream(filename, append, config);
    sink_ = new ProtobufPacketStreamDataSink<T>(stream, true);
  }
  virtual ~ProtobufFileDataSink() { delete sink_; }
  bool put_next(T* data, bool adopt_data = false) override {
    return sink_->put_next(data, adopt_data); }
  static config_type default_options() {
    return config_type::default_instance(); }
private:
  ProtobufPacketStreamDataSink<T>* sink_ = nullptr;
};

// *****************************************************************************
//
// ChaninedFile(RandomAccess)DataSource - chain a number of file-based
// (RandomAccess)DataSource together
//
// *****************************************************************************

template<typename DST> class DataSourceOpener
{
public:
  CALIN_TYPEALIAS(data_source_type, DST);
  virtual ~DataSourceOpener() { }
  virtual unsigned num_sources() = 0;
  virtual DST* open(unsigned isource) = 0;
};

template<typename DST> class FileOpener: public DataSourceOpener<DST>
{
public:
  using DataSourceOpener<DST>::data_source_type;

  FileOpener(const std::vector<std::string>& filenames):
    DataSourceOpener<DST>(), filenames_(filenames) { /* nothing to see here */ }
  virtual ~FileOpener() { /* nothing to see here */ }
  unsigned num_sources() override { return filenames_.size(); }
  DST* open(unsigned isource) override {
    if(isource<filenames_.size())return open_filename(filenames_[isource]);
    return nullptr;
  }
  std::string filename(unsigned isource) {
    if(isource<filenames_.size())return filenames_[isource];
    return {};
  }
  virtual DST* open_filename(const std::string& filename) = 0;
protected:
  std::vector<std::string> filenames_;
};

template<typename DST> class BasicChainedDataSource: public DST
{
public:
  CALIN_TYPEALIAS(data_type, typename DST::data_type);

  BasicChainedDataSource(DataSourceOpener<DST>* opener,
    bool adopt_opener = false):
    DST(), opener_(opener), adopt_opener_(adopt_opener)
    {
      open_file();
    }

  ~BasicChainedDataSource()
  {
    if(adopt_opener_)delete opener_;
    delete source_;
  }

  data_type* get_next() override
  {
    while(isource_ < opener_->num_sources())
    {
      if(data_type* next = source_->get_next())return next;
      ++isource_;
      open_file();
    }
    return nullptr;
  }

  unsigned source_index() const { return isource_; }

protected:
  virtual void open_file()
  {
    if(source_)delete source_;
    source_ = nullptr;
    if(isource_ < opener_->num_sources())
    {
      source_ = opener_->open(isource_);
      // Throw exception if opener doesn't throw its own
      if(source_ == nullptr)throw std::runtime_error("Could not open source " +
        std::to_string(isource_));
    }
  }

  DataSourceOpener<DST>* opener_ = nullptr;
  bool adopt_opener_ = false;
  unsigned isource_ = 0;
  DST* source_ = nullptr;
};

template<typename T> using ChainedDataSource =
  BasicChainedDataSource<DataSource<T> >;

template<typename RADST> class BasicChaninedRandomAccessDataSource:
  public BasicChainedDataSource<RADST>
{
public:
  CALIN_TYPEALIAS(data_type, typename RADST::data_type);

  BasicChaninedRandomAccessDataSource(DataSourceOpener<RADST>* opener,
    bool adopt_opener = false):
    BasicChainedDataSource<RADST>(opener, adopt_opener)
  {
    // Base class can't call virtual open_file function so we complete it
    if(source_)add_source_index(source_);
  }

  virtual ~BasicChaninedRandomAccessDataSource() {
    /* nothing to see here */ }

  virtual uint64_t size() override
  {
    for(unsigned i = chained_file_index_.size(); i<opener_->num_sources(); i++)
    {
      std::unique_ptr<RADST> source(opener_->open(i));
      // Throw exception in case opener doesn't do so itself
      if(!source)throw std::runtime_error("Could not open file " +
        std::to_string(i));
      add_source_index(source.get());
    }

    if(opener_->num_sources())return chained_file_index_.back();
    else return 0;
  }

  void set_next_index(uint64_t next_index) override
  {
    if(opener_->num_sources() == 0)return;

    if(next_index < chained_file_index_.back())
    {
      // next_index belongs to one of the files we have already indexed
      auto file_index = std::upper_bound(chained_file_index_.begin(),
        chained_file_index_.end(), next_index);
      unsigned new_isource = int(file_index-chained_file_index_.begin());
      if(new_isource != isource_)
      {
        isource_ = new_isource;
        BasicChainedDataSource<RADST>::open_file();
      }
    }
    else
    {
      while(isource_<opener_->num_sources() and
        chained_file_index_.back()<=next_index)
      {
        ++isource_;
        open_file();
      }
    }

    if(source_)
    {
      if(isource_ == 0)
        source_->set_next_index(next_index);
      else
        source_->set_next_index(next_index - chained_file_index_[isource_-1]);
    }
  }

private:
  using BasicChainedDataSource<RADST>::source_;
  using BasicChainedDataSource<RADST>::isource_;
  using BasicChainedDataSource<RADST>::opener_;

  void add_source_index(RADST* source)
  {
    if(chained_file_index_.empty())
      chained_file_index_.push_back(source->size());
    else
      chained_file_index_.push_back(source->size()+chained_file_index_.back());
  }

  void open_file() override
  {
    BasicChainedDataSource<RADST>::open_file();
    if(source_ and isource_ >= chained_file_index_.size())
    {
      add_source_index(source_);
    }
  }

  std::vector<uint64_t> chained_file_index_;
};

} } } // namespace calin::io::data_source
