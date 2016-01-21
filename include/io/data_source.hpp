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

#include <io/packet_stream.hpp>
#include <io/packet_stream.pb.h>

namespace calin { namespace io { namespace data_source {

template<typename T> class DataSource
{
public:
  virtual ~DataSource() { }
  virtual T* get_next() = 0;
};

template<typename T> class RandomAccessDataSource: public DataSource<T>
{
public:
  virtual ~RandomAccessDataSource() { }
  virtual T* get_next() = 0;
  virtual uint64_t size() = 0;
  virtual void set_next_index(uint64_t next_index) = 0;
};

template<typename T> class DataSink
{
public:
  virtual ~DataSink() { }
  virtual bool put_next(T* d) = 0;
};

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
  bool put_next(T* d) override {
    return stream_->put_packet(d->SerializeAsString()); }
private:
  calin::io::packet_stream::PacketOutStream* stream_ = nullptr;
  bool adopt_stream_;
};

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
  bool put_next(T* d) override { return sink_->put_next(d); }
  static config_type default_options() {
    return config_type::default_instance(); }
private:
  ProtobufPacketStreamDataSink<T>* sink_ = nullptr;
};

} } } // namespace calin::io::data_source
