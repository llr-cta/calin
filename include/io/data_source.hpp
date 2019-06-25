/*

   calin/io/data_source.hpp -- Stephen Fegan -- 2016-01-15

   A supplier of generic data, probably Protobufs

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

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>

#include <pattern/delegation.hpp>
#include <io/packet_stream.hpp>
#include <io/packet_stream.pb.h>

namespace calin { namespace io { namespace data_source {

template<typename T> class DataSource
{
public:
  CALIN_TYPEALIAS(data_type, T);
  virtual ~DataSource() { }
  virtual T* get_next(uint64_t& seq_index_out,
    google::protobuf::Arena** arena = nullptr) = 0;
};

template<typename T> class DelegatedDataSource:
  public virtual DataSource<T>, public calin::pattern::delegation::Delegator<DataSource<T> >
{
public:
  CALIN_TYPEALIAS(data_type, T);
  DelegatedDataSource(DataSource<T>* delegate, bool adopt_delegate = false):
    DataSource<T>(),
    pattern::delegation::Delegator<DataSource<T> >(delegate, adopt_delegate) { }
  virtual ~DelegatedDataSource() { }
  T* get_next(uint64_t& seq_index_out,
    google::protobuf::Arena** arena = nullptr) override
  {
    return this->delegate_->get_next(seq_index_out,arena);
  }
};


template<typename T> class RandomAccessDataSource: public virtual DataSource<T>
{
public:
  CALIN_TYPEALIAS(data_type, T);
  virtual ~RandomAccessDataSource() { }
  virtual T* get_next(uint64_t& seq_index_out,
    google::protobuf::Arena** arena = nullptr) override = 0;
  virtual uint64_t size() = 0;
  virtual void set_next_index(uint64_t next_index) = 0;
};

template<typename T> class DataSink
{
public:
  CALIN_TYPEALIAS(data_type, T);
  virtual ~DataSink() { }
  virtual bool put_next(T* data, uint64_t seq_index,
    google::protobuf::Arena* arena = nullptr, bool adopt_data = false) = 0;
};

template<typename T> class DataSourceFactory
{
public:
  virtual ~DataSourceFactory() { }
  virtual DataSource<T>* new_data_source() = 0;
};

template<typename T> class DataSinkFactory
{
public:
  virtual ~DataSinkFactory() { }
  virtual DataSink<T>* new_data_sink() = 0;
};

template<typename T> class VectorDataSourceFactory: public DataSourceFactory<T>
{
public:
  VectorDataSourceFactory(const std::vector<DataSource<T>*>& sources,
      bool adopt_sources = false):
    DataSourceFactory<T>(), sources_(sources), adopt_sources_(adopt_sources) { }
  virtual ~VectorDataSourceFactory() {
    if(adopt_sources_)for(auto* src: sources_)delete src;
  }
  DataSource<T>* new_data_source() override {
    unsigned isource = isource_.fetch_add(1,std::memory_order_relaxed);
    if(isource<sources_.size()) {
      return new DelegatedDataSource<T>(sources_[isource], false);
    } else {
      return nullptr;
    }
  }

private:
  std::vector<DataSource<T>*> sources_;
  bool adopt_sources_ = false;
  std::atomic<uint_fast32_t> isource_ { 0 };
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
  T* get_next(uint64_t& seq_index_out,
      google::protobuf::Arena** arena = nullptr) override {
    std::string serialized_d;
    if(!stream_->get_packet(serialized_d, seq_index_out))return nullptr;
    T* d = nullptr;
    T* delete_d = nullptr;
    google::protobuf::Arena* delete_arena = nullptr;
    if(arena) {
      if(!*arena)*arena = delete_arena = new google::protobuf::Arena;
      d = google::protobuf::Arena::CreateMessage<T>(*arena);
    }
    else d = delete_d = new T;
    if(!d->ParseFromString(serialized_d)){
      delete delete_d; delete delete_arena; return nullptr; }
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
  bool put_next(T* data, uint64_t seq_index,
    google::protobuf::Arena* arena = nullptr, bool adopt_data = false) override
  {
    bool retval = stream_->put_packet(data->SerializeAsString(), seq_index);
    if(adopt_data) {
      if(arena)delete arena;
      else delete data;
    }
    return retval;
  }
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
  T* get_next(uint64_t& seq_index_out,
      google::protobuf::Arena** arena = nullptr) override {
    return source_->get_next(seq_index_out, arena); }
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
  bool put_next(T* data, uint64_t seq_index,
    google::protobuf::Arena* arena = nullptr, bool adopt_data = false) override
  {
    return sink_->put_next(data, seq_index, arena, adopt_data);
  }
  static config_type default_options() {
    return config_type::default_instance(); }
private:
  ProtobufPacketStreamDataSink<T>* sink_ = nullptr;
};

} } } // namespace calin::io::data_source
