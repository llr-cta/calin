%module (package="calin") google_protobuf

%{
  #include<google/protobuf/message.h>
  #include<google/protobuf/descriptor.h>
  #define SWIG_FILE_WITH_INIT
%}

%include "package_wide_definitions.i"

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

  bool ParseFromString(const std::string& data);
  bool ParsePartialFromString(const std::string& data);
  std::string SerializeAsString();
  std::string SerializePartialAsString();
};

class Descriptor
{
 public:
  ~Descriptor();
  std::string name() const;
};
   
} }
