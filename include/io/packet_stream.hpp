/*

   calin/io/packet_stream.hpp -- Stephen Fegan -- 2016-01-15

   A supplier of binary packtets, probably serealized Protobufs, possibly
   compressed

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

#include <provenance/chronicle.pb.h>
#include <calin_global_definitions.hpp>
#include <io/packet_stream.pb.h>

namespace google { namespace protobuf { namespace io {
  class ZeroCopyInputStream;
  class CodedInputStream;
  class ZeroCopyOutputStream;
  class CodedOutputStream;
} } } // namespace google::protobuf::io;

namespace calin { namespace io { namespace packet_stream {

class PacketInStream
{
public:
  virtual ~PacketInStream();
  virtual bool get_packet(std::string& packet_out, uint64_t& seq_index_out) = 0;
};

class PacketOutStream
{
public:
  virtual ~PacketOutStream();
  virtual bool put_packet(const std::string& packet, uint64_t seq_index) = 0;
};

class FramedZeroCopyPacketInStream final: public PacketInStream
{
public:
  FramedZeroCopyPacketInStream(
    google::protobuf::io::ZeroCopyInputStream* instream,
    bool adopt_instream = false);
  virtual ~FramedZeroCopyPacketInStream();
  bool get_packet(std::string& packet_out, uint64_t& seq_index_out) override;
private:
  google::protobuf::io::ZeroCopyInputStream* instream_ = nullptr;
  bool adopt_instream_ = false;
  google::protobuf::io::CodedInputStream* coded_instream_ = nullptr;
};

class FramedZeroCopyPacketOutStream final : public PacketOutStream
{
public:
  FramedZeroCopyPacketOutStream(
    google::protobuf::io::ZeroCopyOutputStream* outstream,
    bool adopt_outstream = false);
  ~FramedZeroCopyPacketOutStream();
  bool put_packet(const std::string& packet, uint64_t seq_index) override;
private:
  google::protobuf::io::ZeroCopyOutputStream* outstream_ = nullptr;
  bool adopt_outstream_ = false;
  google::protobuf::io::CodedOutputStream* coded_outstream_ = nullptr;
};

class CompressedPacketInStream final: public PacketInStream
{
public:
  CALIN_TYPEALIAS(config_type,
    ix::io::packet_stream::CompressedPacketInStreamOptions);
  CompressedPacketInStream(PacketInStream* upstream,
    const config_type& config = config_type::default_instance(),
    bool adopt_upstream = false);
  virtual ~CompressedPacketInStream();
  bool get_packet(std::string& packet_out, uint64_t& seq_index_out) override;
private:
  PacketInStream* upstream_;
  bool adopt_upstream_ = false;
  config_type config_;
};

class CompressedPacketOutStream final: public PacketOutStream
{
public:
  CALIN_TYPEALIAS(config_type,
    ix::io::packet_stream::CompressedPacketOutStreamOptions);
  CompressedPacketOutStream(PacketOutStream* downstream,
    const config_type& config = config_type::default_instance(),
    bool adopt_downstream = false);
  virtual ~CompressedPacketOutStream();
  bool put_packet(const std::string& packet, uint64_t seq_index) override;
private:
  PacketOutStream* downstream_;
  bool adopt_downstream_ = false;
  config_type config_;
};

class FramedFilePacketInStream final: public PacketInStream
{
public:
  CALIN_TYPEALIAS(config_type,
    ix::io::packet_stream::CompressedPacketInStreamOptions);
  FramedFilePacketInStream(const std::string& filename,
    const config_type& config = config_type::default_instance());
  virtual ~FramedFilePacketInStream();
  bool get_packet(std::string& packet_out, uint64_t& seq_index_out) override;
private:
  PacketInStream* upstream_ = nullptr;
  calin::ix::provenance::chronicle::FileIORecord* file_record_ = nullptr;
};

class FramedFilePacketOutStream final: public PacketOutStream
{
public:
  CALIN_TYPEALIAS(config_type,
    ix::io::packet_stream::CompressedPacketOutStreamOptions);
  FramedFilePacketOutStream(const std::string& filename, bool append = false,
    const config_type& config = config_type::default_instance());
  virtual ~FramedFilePacketOutStream();
  bool put_packet(const std::string& packet, uint64_t seq_index) override;
private:
  PacketOutStream* downstream_ = nullptr;
  calin::ix::provenance::chronicle::FileIORecord* file_record_ = nullptr;
};

} } } // namespace calin::io::packet_stream
