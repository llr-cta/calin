/*

   calin/io/packet_stream.cpp -- Stephen Fegan -- 2016-01-15

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

#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/gzip_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/zero_copy_stream_impl_lite.h>

#include <io/packet_stream.hpp>

using namespace google::protobuf::io;
using namespace calin::io::packet_stream;

PacketInStream::~PacketInStream()
{
  // nothing to see here
}

PacketOutStream::~PacketOutStream()
{
  // nothing to see here
};

namespace {

bool copy_stream(ZeroCopyOutputStream& stream_out,
                 ZeroCopyInputStream& stream_in)
{
  int size_in = 0;
  const char* buffer_in = nullptr;
  while(stream_in.Next((const void**)&buffer_in, &size_in))
  {
    while(size_in)
    {
      int size_out = 0;
      char* buffer_out = nullptr;
      if(!stream_out.Next((void**)&buffer_out, &size_out))return false;
      if(size_out > size_in) {
        memcpy(buffer_out, buffer_in, size_in);
        stream_out.BackUp(size_out-size_in);
        size_in = 0;
      } else {
        memcpy(buffer_out, buffer_in, size_out);
        size_in -= size_out;
      }
    }
  }
  return true;
}

} // anonymous namespace

ZeroCopyCodedPacketInStream::
ZeroCopyCodedPacketInStream(ZeroCopyInputStream* instream, bool adopt_instream):
  PacketInStream(), instream_(instream), adopt_instream_(adopt_instream),
  coded_instream_(new CodedInputStream(instream))
{
  // nothing to see here
}

ZeroCopyCodedPacketInStream::~ZeroCopyCodedPacketInStream()
{
  delete coded_instream_;
  if(adopt_instream_)delete instream_;
}

bool ZeroCopyCodedPacketInStream::getPacket(std::string& packet_out)
{
  uint32_t packet_size = 0;
  packet_out.clear();
  return coded_instream_->ReadLittleEndian32(&packet_size) and
    coded_instream_->ReadString(&packet_out, packet_size);
}

ZeroCopyCodedPacketOutStream::
ZeroCopyCodedPacketOutStream(ZeroCopyOutputStream* outstream,
    bool adopt_outstream):
  PacketOutStream(), outstream_(outstream), adopt_outstream_(adopt_outstream),
  coded_outstream_(new CodedOutputStream(outstream))
{
  // nothing to see here
}

ZeroCopyCodedPacketOutStream::~ZeroCopyCodedPacketOutStream()
{
  delete coded_outstream_;
  if(adopt_outstream_)delete outstream_;
}

bool ZeroCopyCodedPacketOutStream::putPacket(const std::string& packet)
{
  coded_outstream_->WriteLittleEndian32(packet.size());
  coded_outstream_->WriteString(packet);
  return !coded_outstream_->HadError();
}

CompressedPacketInStream::
CompressedPacketInStream(PacketInStream* upstream, bool adopt_upstream):
  PacketInStream(), upstream_(upstream), adopt_upstream_(adopt_upstream)
{
  // nothing to see here
}

CompressedPacketInStream::~CompressedPacketInStream()
{
  if(adopt_upstream_)delete upstream_;
}

bool CompressedPacketInStream::getPacket(std::string& packet_out)
{
  packet_out.clear();
  std::string packet_in;
  if(!upstream_->getPacket(packet_in))return false;;
  ArrayInputStream zstream_in(packet_in.data(), packet_in.size());
  GzipInputStream stream_in(&zstream_in);
  StringOutputStream stream_out(&packet_out);
  if(!copy_stream(stream_out, stream_in))return false;
  return true;
}

CompressedPacketOutStream::
CompressedPacketOutStream(PacketOutStream* downstream, bool adopt_downstream)
{
  // nothing to see here
}

CompressedPacketOutStream::~CompressedPacketOutStream()
{
  if(adopt_downstream_)delete downstream_;
}

bool CompressedPacketOutStream::putPacket(const std::string& packet)
{
  ArrayInputStream stream_in(packet.data(), packet.size());
  std::string string_out;
  if(true)
  {
    StringOutputStream zstream_out(&string_out);
    GzipOutputStream stream_out(&zstream_out);
    if(!copy_stream(stream_out, stream_in))return false;
  }
  return downstream_->putPacket(string_out);
}
