%module (package="calin") google_protobuf

%{
  #include<google/protobuf/message.h>
  #include<google/protobuf/descriptor.h>
  #define SWIG_FILE_WITH_INIT
%}

%include "package_wide_definitions.i"

%import<stdint.i>

#define int32 int32_t
#define uint32 uint32_t
#define int64 int64_t
#define uint64 uint64_t

namespace google { namespace protobuf {

class Message
{
 public:
  Message();
  virtual ~Message();

  //google::protobuf::Message* New() const override;
  virtual void CopyFrom(const Message & from);
  virtual void MergeFrom(const Message & from);
  virtual int SpaceUsed() const;
};

class Descriptor
{
 public:

};
   
} }
