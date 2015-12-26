//-*-mode:swig;-*-

/* 

   calin/google_protobuf.i -- Stephen Fegan -- 2015-12-10

   SWIG interface file for supported elements from protobuf Message
   and Descriptoy

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin") google_protobuf

%{
  #include<google/protobuf/message.h>
  #include<google/protobuf/descriptor.h>
  #define SWIG_FILE_WITH_INIT
%}

%import "package_wide_definitions.i"
%include "calin_typemaps.i"

%import<stdint.i>
%import<std_string.i>

#define int32 int32_t
#define uint32 uint32_t
#define int64 int64_t
#define uint64 uint64_t

namespace google { namespace protobuf {

%newobject Message::New() const;
%nodefaultctor Message;
%nodefaultctor Descriptor;

class Message
{
 public:
  //Message();
  ~Message();

  google::protobuf::Message* New() const;
  void CopyFrom(const Message & from);
  void MergeFrom(const Message & from);
  int SpaceUsed() const;

  std::string DebugString() const;
  std::string ShortDebugString() const;
  std::string GetTypeName() const;
  void Clear();
  bool IsInitialized();
  int ByteSize();

  bool ParseFromString(const std::string& CALIN_BYTES_IN);
  bool ParsePartialFromString(const std::string& CALIN_BYTES_IN);
  %extend {
    void SerializeAsString(std::string& CALIN_BYTES_OUT) {
      CALIN_BYTES_OUT = $self->SerializeAsString(); }
    void SerializePartialAsString(std::string& CALIN_BYTES_OUT) {
      CALIN_BYTES_OUT = $self->SerializePartialAsString(); }
  }
};

class Descriptor
{
 public:
  ~Descriptor();
  std::string name() const;
};
   
} }
