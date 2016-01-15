/*

   calin/io/packet_stream.hpp -- Stephen Fegan -- 2016-01-15

   A supplier of binary packtets, probably serealized Protobufs, possibly
   compressed

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
  virtual bool getPacket(std::string& packet_out) = 0;
};

class PacketOutStream
{
public:
  virtual ~PacketOutStream();
  virtual bool putPacket(const std::string& packet) = 0;
};

class ZeroCopyCodedPacketInStream final: public PacketInStream
{
public:
  ZeroCopyCodedPacketInStream(
    google::protobuf::io::ZeroCopyInputStream* instream,
    bool adopt_instream = false);
  virtual ~ZeroCopyCodedPacketInStream();
  bool getPacket(std::string& packet_out) override;
private:
  google::protobuf::io::ZeroCopyInputStream* instream_ = nullptr;
  bool adopt_instream_ = false;
  google::protobuf::io::CodedInputStream* coded_instream_ = nullptr;
};

class ZeroCopyCodedPacketOutStream final : public PacketOutStream
{
public:
  ZeroCopyCodedPacketOutStream(
    google::protobuf::io::ZeroCopyOutputStream* outstream,
    bool adopt_outstream = false);
  ~ZeroCopyCodedPacketOutStream();
  bool putPacket(const std::string& packet) override;
private:
  google::protobuf::io::ZeroCopyOutputStream* outstream_ = nullptr;
  bool adopt_outstream_ = false;
  google::protobuf::io::CodedOutputStream* coded_outstream_ = nullptr;
};

class CompressedPacketInStream final: public PacketInStream
{
public:
  CompressedPacketInStream(PacketInStream* upstream,
    bool adopt_upstream = false);
  virtual ~CompressedPacketInStream();
  bool getPacket(std::string& packet_out) override;
private:
  PacketInStream* upstream_;
  bool adopt_upstream_ = false;
};

class CompressedPacketOutStream final: public PacketOutStream
{
public:
  CompressedPacketOutStream(PacketOutStream* downstream,
    bool adopt_downstream = false);
  virtual ~CompressedPacketOutStream();
  bool putPacket(const std::string& packet) override;
private:
  PacketOutStream* downstream_;
  bool adopt_downstream_ = false;
};

} } } // namespace calin::io::packet_stream
