/*

   calin/io/packet_stream.cpp -- Stephen Fegan -- 2016-01-15

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

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdexcept>

#include <provenance/chronicle.hpp>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/gzip_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/zero_copy_stream_impl_lite.h>

#include <util/log.hpp>
#include <io/packet_stream.hpp>

using namespace google::protobuf::io;
using namespace calin::util::log;
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

FramedZeroCopyPacketInStream::
FramedZeroCopyPacketInStream(ZeroCopyInputStream* instream, bool adopt_instream):
  PacketInStream(), instream_(instream), adopt_instream_(adopt_instream),
  coded_instream_(new CodedInputStream(instream))
{
  // nothing to see here
}

FramedZeroCopyPacketInStream::~FramedZeroCopyPacketInStream()
{
  delete coded_instream_;
  if(adopt_instream_)delete instream_;
}

bool FramedZeroCopyPacketInStream::
get_packet(std::string& packet_out, uint64_t& seq_index_out)
{
  uint32_t packet_size = 0;
  packet_out.clear();
  if(!coded_instream_->ReadLittleEndian32(&packet_size))
    return false; // Must assume EOF since API does not distinguish errors
  google::protobuf::uint64 google_seq_index_out;
  if(!coded_instream_->ReadLittleEndian64(&google_seq_index_out))
    throw std::runtime_error("Error reading sequence index from input stream.");
  seq_index_out = google_seq_index_out;
  if(!coded_instream_->ReadString(&packet_out, packet_size))
    throw std::runtime_error("Error reading packet data from input stream.");
  if(packet_out.size()!=packet_size)
    throw std::runtime_error("Packet framing error: expected "
      + std::to_string(packet_size) + " bytes, got "
      + std::to_string(packet_out.size()) + " bytes.");
  return true;
}

FramedZeroCopyPacketOutStream::
FramedZeroCopyPacketOutStream(ZeroCopyOutputStream* outstream,
    bool adopt_outstream):
  PacketOutStream(), outstream_(outstream), adopt_outstream_(adopt_outstream),
  coded_outstream_(new CodedOutputStream(outstream))
{
  // nothing to see here
}

FramedZeroCopyPacketOutStream::~FramedZeroCopyPacketOutStream()
{
  delete coded_outstream_;
  if(adopt_outstream_)delete outstream_;
}

bool FramedZeroCopyPacketOutStream::
put_packet(const std::string& packet, uint64_t seq_index)
{
  coded_outstream_->WriteLittleEndian32(packet.size());
  if(coded_outstream_->HadError())return false;
  coded_outstream_->WriteLittleEndian64(seq_index);
  if(coded_outstream_->HadError())
    throw std::runtime_error("Could not sequence index data, "
      "framing on output stream will be incorrect from this point.");
  coded_outstream_->WriteString(packet);
  if(coded_outstream_->HadError())
    throw std::runtime_error("Could not write packet data, "
      "framing on output stream will be incorrect from this point.");
  return true;
}

CompressedPacketInStream::
CompressedPacketInStream(PacketInStream* upstream, const config_type& config,
  bool adopt_upstream):
  PacketInStream(), upstream_(upstream), adopt_upstream_(adopt_upstream),
  config_(config)
{
  // nothing to see here
}

CompressedPacketInStream::~CompressedPacketInStream()
{
  if(adopt_upstream_)delete upstream_;
}

bool CompressedPacketInStream::
get_packet(std::string& packet_out, uint64_t& seq_index_out)
{
  if(!config_.use_compression())
    return upstream_->get_packet(packet_out, seq_index_out);
  packet_out.clear();
  std::string packet_in;
  if(!upstream_->get_packet(packet_in, seq_index_out))return false;;
  ArrayInputStream zstream_in(packet_in.data(), packet_in.size());
  GzipInputStream stream_in(&zstream_in);
  bool copy_good = true;
  if(true)
  {
    StringOutputStream stream_out(&packet_out);
    copy_good = copy_stream(stream_out, stream_in);
  }
  if(stream_in.ZlibErrorCode() < 0)
    throw std::runtime_error(std::string("Zlib decompression error: ") +
      stream_in.ZlibErrorMessage());
  if(!copy_good)
    throw std::runtime_error("Error writing to string stream.");
  return true;
}

CompressedPacketOutStream::
CompressedPacketOutStream(PacketOutStream* downstream,
  const config_type& config, bool adopt_downstream):
  PacketOutStream(), downstream_(downstream),
  adopt_downstream_(adopt_downstream), config_(config)
{
  // nothing to see here
}

CompressedPacketOutStream::~CompressedPacketOutStream()
{
  if(adopt_downstream_)delete downstream_;
}

bool CompressedPacketOutStream::
put_packet(const std::string& packet, uint64_t seq_index)
{
  if(!config_.use_compression())
    return downstream_->put_packet(packet, seq_index);
  ArrayInputStream stream_in(packet.data(), packet.size());
  std::string packet_out;
  StringOutputStream zstream_out(&packet_out);
  GzipOutputStream stream_out(&zstream_out);
  bool copy_good = copy_stream(stream_out, stream_in);
  stream_out.Flush();
  if(stream_out.ZlibErrorCode() < 0)
    throw std::runtime_error(std::string("Zlib decompression error: ") +
      stream_out.ZlibErrorMessage());
  if(!copy_good)
    throw std::runtime_error("Error writing to decompression stream.");
  return downstream_->put_packet(packet_out, seq_index);
}

FramedFilePacketInStream::FramedFilePacketInStream(const std::string& filename,
  const config_type& config): PacketInStream()
{
  int fd = open(filename.c_str(), O_RDONLY);
  if(fd==-1)throw std::runtime_error(std::string("Error opening \"") +
    filename + "\": " + std::strerror(errno));
  calin::provenance::chronicle::register_file_open(filename,
    calin::ix::provenance::chronicle::AT_READ, __PRETTY_FUNCTION__);
  auto* filestream = new FileInputStream(fd);
  filestream->SetCloseOnDelete(true);
  upstream_ = new FramedZeroCopyPacketInStream(filestream, true);
  if(config.use_compression())upstream_ =
    new CompressedPacketInStream(upstream_, config, true);
}

FramedFilePacketInStream::~FramedFilePacketInStream()
{
  delete upstream_;
}

bool FramedFilePacketInStream::
get_packet(std::string& packet_out, uint64_t& seq_index_out)
{
  if(!upstream_)throw std::runtime_error("File is not open");
  return upstream_->get_packet(packet_out, seq_index_out);
}

FramedFilePacketOutStream::FramedFilePacketOutStream(const std::string& filename,
  bool append, const config_type& config): PacketOutStream()
{
  constexpr int mode = S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH|S_IWOTH;
  int fd = open(filename.c_str(), O_WRONLY|O_CREAT|(append?O_APPEND:O_TRUNC), mode);
  if(fd==-1)throw std::runtime_error(std::string("Error opening \"") +
    filename + "\" for writing: " + std::strerror(errno));
  calin::provenance::chronicle::register_file_open(filename,
    append?calin::ix::provenance::chronicle::AT_WRITE:
      calin::ix::provenance::chronicle::AT_TRUNCATE, __PRETTY_FUNCTION__);
  auto* filestream = new FileOutputStream(fd);
  filestream->SetCloseOnDelete(true);
  downstream_ = new FramedZeroCopyPacketOutStream(filestream, true);
  if(config.use_compression())downstream_ =
    new CompressedPacketOutStream(downstream_, config, true);
}

FramedFilePacketOutStream::~FramedFilePacketOutStream()
{
  delete downstream_;
}

bool FramedFilePacketOutStream::
put_packet(const std::string& packet, uint64_t seq_index)
{
  if(!downstream_)throw std::runtime_error("File is not open");
  return downstream_->put_packet(packet, seq_index);
}
