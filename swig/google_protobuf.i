//-*-mode:swig;-*-

/*

   calin/google_protobuf.i -- Stephen Fegan -- 2015-12-10

   SWIG interface file for supported elements from protobuf Message
   and Descriptor

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

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

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
%nodefaultctor Arena;

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

class Arena
{
public:
  ~Arena();
  uint64 SpaceUsed() const;
  uint64 Reset();
};

} }

%typemap(in, numinputs=0) google::protobuf::Arena** CALIN_NONNULL_ARENA_OUTPUT
  (google::protobuf::Arena* temp = nullptr) {
    // typemap(in) google::protobuf::Arena** OUTPUT
    $1 = &temp;
}

%typemap(argout)
  google::protobuf::Arena** CALIN_NONNULL_ARENA_OUTPUT {
    // typemap(argout) google::protobuf::Arena** OUTPUT
    if(*$1 == nullptr) {
      PyErr_Format(PyExc_TypeError,
                   "Memory management error: no Arena returned.");
      SWIG_fail;
    }
    %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(*$1), $*1_descriptor, SWIG_POINTER_OWN));
}
