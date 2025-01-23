//-*-mode:swig;-*-

/*

   calin/google_protobuf.i -- Stephen Fegan -- 2015-12-10

   SWIG interface file for supported elements from protobuf Message
   and Descriptor

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin") google_protobuf

%{
  #include<string>
  #include<vector>
  #include<google/protobuf/message.h>
  #include<google/protobuf/descriptor.h>
  #include<google/protobuf/util/json_util.h>
  #include<google/protobuf/util/message_differencer.h>
  #define SWIG_FILE_WITH_INIT
%}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import"calin_stdint.i"
%import<std_string.i>

#define int32 int32_t
#define uint32 uint32_t
#define int64 int64_t
#define uint64 uint64_t

namespace google { namespace protobuf {

  //%newobject Message::New() const;
  %nodefaultctor Message;
  %nodefaultctor Descriptor;
  %nodefaultctor FieldDescriptor;
  //%nodefaultctor Arena;

  class Descriptor;

  class FieldDescriptor
  {
  public:
    enum Type {
      TYPE_DOUBLE = 1,
      TYPE_FLOAT = 2,
      TYPE_INT64 = 3,
      TYPE_UINT64 = 4,
      TYPE_INT32 = 5,
      TYPE_FIXED64 = 6,
      TYPE_FIXED32 = 7,
      TYPE_BOOL = 8,
      TYPE_STRING = 9,
      TYPE_GROUP = 10,
      TYPE_MESSAGE = 11,
      TYPE_BYTES = 12,
      TYPE_UINT32 = 13,
      TYPE_ENUM = 14,
      TYPE_SFIXED32 = 15,
      TYPE_SFIXED64 = 16,
      TYPE_SINT32 = 17,
      TYPE_SINT64 = 18,
      MAX_TYPE = 18
    };

    enum CppType {
      CPPTYPE_INT32 = 1,
      CPPTYPE_INT64 = 2,
      CPPTYPE_UINT32 = 3,
      CPPTYPE_UINT64 = 4,
      CPPTYPE_DOUBLE = 5,
      CPPTYPE_FLOAT = 6,
      CPPTYPE_BOOL = 7,
      CPPTYPE_ENUM = 8,
      CPPTYPE_STRING = 9,
      CPPTYPE_MESSAGE = 10,
      MAX_CPPTYPE = 10
    };

    enum Label {
      LABEL_OPTIONAL = 1,
      LABEL_REQUIRED = 2,
      LABEL_REPEATED = 3,
      MAX_LABEL = 3
    };

    std::string name() const;
    std::string full_name() const;
    std::string json_name() const;
    bool is_extension() const;
    int	number() const;
    std::string lowercase_name() const;
    std::string camelcase_name() const;

    Type type() const;
    std::string type_name() const;
    CppType cpp_type() const;
    std::string cpp_type_name() const;

    Label	label() const;
    bool is_required() const;
    bool is_optional() const;
    bool is_repeated() const;
    bool is_packable() const;
    bool is_packed() const;
    bool is_map() const;

    %extend {
      bool is_double() const { 
        return $self->cpp_type() == google::protobuf::FieldDescriptor::CPPTYPE_DOUBLE;
      }
      bool is_float() const {
        return $self->cpp_type() == google::protobuf::FieldDescriptor::CPPTYPE_FLOAT;
      }
      bool is_int64() const {
        return $self->cpp_type() == google::protobuf::FieldDescriptor::CPPTYPE_INT64;
      }
      bool is_uint64() const {
        return $self->cpp_type() == google::protobuf::FieldDescriptor::CPPTYPE_UINT64;
      }
      bool is_int32() const {
        return $self->cpp_type() == google::protobuf::FieldDescriptor::CPPTYPE_INT32;
      }
      bool is_uint32() const {
        return $self->cpp_type() == google::protobuf::FieldDescriptor::CPPTYPE_UINT32;
      }
      bool is_bool() const {
        return $self->cpp_type() == google::protobuf::FieldDescriptor::CPPTYPE_BOOL;
      }
      bool is_enum() const {
        return $self->cpp_type() == google::protobuf::FieldDescriptor::CPPTYPE_ENUM;
      }
      bool is_string() const {
        return $self->type() == google::protobuf::FieldDescriptor::TYPE_STRING;
      }
      bool is_bytes() const {
        return $self->type() == google::protobuf::FieldDescriptor::TYPE_BYTES;
      }
      bool is_message() const {
        return $self->cpp_type() == google::protobuf::FieldDescriptor::CPPTYPE_MESSAGE;
      }
      bool is_simple_message() const {
        return $self->cpp_type() == google::protobuf::FieldDescriptor::CPPTYPE_MESSAGE and not $self->is_repeated();
      }
      bool is_repeated_message() const {
        return $self->cpp_type() == google::protobuf::FieldDescriptor::CPPTYPE_MESSAGE and $self->is_repeated();
      }
      bool is_simple_pod() const {
        return $self->cpp_type() != google::protobuf::FieldDescriptor::CPPTYPE_MESSAGE and not $self->is_repeated();
      }
      bool is_repeated_pod() const {
        return $self->cpp_type() != google::protobuf::FieldDescriptor::CPPTYPE_MESSAGE and $self->is_repeated();
      }
    }

    bool has_optional_keyword() const;
    bool has_presence() const;
    int	index() const;
    bool has_default_value() const;
    bool has_json_name() const;

    const Descriptor * containing_type() const;

    const Descriptor * message_type() const;

    std::string	DebugString() const;

    static CppType TypeToCppType(Type type);
    static std::string TypeName(Type type);
    static std::string CppTypeName(CppType cpp_type);
    static bool IsTypePackable(Type field_type);
  };

  class Descriptor
  {
  public:
    std::string name() const;
    std::string full_name() const;
    int index() const;

    std::string	DebugString() const;

    const Descriptor* containing_type() const;
    
    int field_count() const;
    const FieldDescriptor* field(int index) const;
    const FieldDescriptor* FindFieldByName(const std::string& name) const;
    const FieldDescriptor* FindFieldByNumber(int number) const;
    const FieldDescriptor* FindFieldByLowercaseName(const std::string& name) const;
    const FieldDescriptor* FindFieldByCamelcaseName(const std::string& name) const;

    %extend {
      std::vector<std::string> GetSimpleMessages() {
        std::vector<std::string> OUTPUT;
        for(int i=0; i<$self->field_count(); i++) {
          const google::protobuf::FieldDescriptor* f = $self->field(i);
          if(f->type() == google::protobuf::FieldDescriptor::TYPE_MESSAGE and not f->is_repeated()) OUTPUT.push_back(f->name());
        }
        return OUTPUT;
      }
      std::vector<std::string> GetRepeatedMessages() {
        std::vector<std::string> OUTPUT;
        for(int i=0; i<$self->field_count(); i++) {
          const google::protobuf::FieldDescriptor* f = $self->field(i);
          if(f->type() == google::protobuf::FieldDescriptor::TYPE_MESSAGE and f->is_repeated()) OUTPUT.push_back(f->name());
        }
        return OUTPUT;
      }
      std::vector<std::string> GetSimplePODs() {
        std::vector<std::string> OUTPUT;
        for(int i=0; i<$self->field_count(); i++) {
          const google::protobuf::FieldDescriptor* f = $self->field(i);
          if(f->type() != google::protobuf::FieldDescriptor::TYPE_MESSAGE and not f->is_repeated()) OUTPUT.push_back(f->name());
        }
        return OUTPUT;
      }
      std::vector<std::string> GetRepeatedPODs() {
        std::vector<std::string> OUTPUT;
        for(int i=0; i<$self->field_count(); i++) {
          const google::protobuf::FieldDescriptor* f = $self->field(i);
          if(f->type() != google::protobuf::FieldDescriptor::TYPE_MESSAGE and f->is_repeated()) OUTPUT.push_back(f->name());
        }
        return OUTPUT;
      }
    }
  };

  class Message
  {
  public:
    //Message();
    ~Message();

    //google::protobuf::Message* New() const;
    void CopyFrom(const Message & from);
    void MergeFrom(const Message & from);
    uint64_t SpaceUsedLong() const;
    %extend {
      uint64_t SpaceUsed() {
        return $self->SpaceUsedLong();
      }
    }

    std::string DebugString() const;
    std::string ShortDebugString() const;
    std::string GetTypeName() const;
    void Clear();
    bool IsInitialized();
    uint64_t ByteSizeLong();
    %extend {
      uint64_t ByteSize() {
        return $self->ByteSizeLong();
      }
    }

    bool ParseFromString(const std::string& CALIN_BYTES_IN);
    bool ParsePartialFromString(const std::string& CALIN_BYTES_IN);
    %extend {
      void SerializeAsString(std::string& CALIN_BYTES_OUT) {
        CALIN_BYTES_OUT = $self->SerializeAsString(); }
      void SerializePartialAsString(std::string& CALIN_BYTES_OUT) {
        CALIN_BYTES_OUT = $self->SerializePartialAsString(); }
    }

    %extend {
      std::string SerializeAsJSON(bool include_empty_fields = false) {
        google::protobuf::util::JsonPrintOptions opt;
        opt.add_whitespace = true;
        opt.always_print_primitive_fields = include_empty_fields;
        std::string s;
        google::protobuf::util::MessageToJsonString(*$self, &s, opt);
        return s;
      }

      bool ParseFromJSON(const std::string& s,
          bool ignore_unknown_fields = false) {
        google::protobuf::util::JsonParseOptions opt;
        opt.ignore_unknown_fields = ignore_unknown_fields;
        return JsonStringToMessage(s, $self, opt).ok();
      }
    }
  };

  class Arena
  {
  public:
    ~Arena();
    uint64_t SpaceUsed() const;
    uint64_t Reset();
  };

} } // namespace google::protobuf

%typemap(in, numinputs=0) google::protobuf::Arena** CALIN_ARENA_OUTPUT
  (google::protobuf::Arena* temp = nullptr) {
    // typemap(in) google::protobuf::Arena** CALIN_ARENA_OUTPUT - google_protobuf.i
    $1 = &temp;
}

%typemap(argout)
  google::protobuf::Arena** CALIN_ARENA_OUTPUT {
    // typemap(argout) google::protobuf::Arena** CALIN_NONNULL_ARENA_OUTPUT - google_protobuf.i
    %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(*$1), $*1_descriptor, SWIG_POINTER_OWN));
}

%template(VectorGoogleProtobufFieldDescriptor) std::vector<const google::protobuf::FieldDescriptor*>;

namespace google { namespace protobuf { namespace util {

  class MessageDifferencer
  {
  public:
    MessageDifferencer();
    ~MessageDifferencer();

    static bool Equals(const google::protobuf::Message& message1, const google::protobuf::Message& message2);
    static bool Equivalent(const google::protobuf::Message& message1, const google::protobuf::Message& message2);
    static bool ApproximatelyEquals(const google::protobuf::Message& message1, const google::protobuf::Message& message2);
    static bool ApproximatelyEquivalent(const google::protobuf::Message& message1, const google::protobuf::Message& message2);

    void TreatAsSet(const FieldDescriptor * field);
    void TreatAsSmartSet(const FieldDescriptor * field);
    void TreatAsList(const FieldDescriptor * field);
    void TreatAsSmartList(const FieldDescriptor * field);
    void TreatAsMap(const FieldDescriptor * field, const FieldDescriptor * key);	
    void TreatAsMapWithMultipleFieldsAsKey(const FieldDescriptor * field, const std::vector< const FieldDescriptor * > & key_fields);
    // void TreatAsMapWithMultipleFieldPathsAsKey(const FieldDescriptor * field, const std::vector< std::vector< const FieldDescriptor * > > & key_field_paths);
    // void TreatAsMapUsingKeyComparator(const FieldDescriptor * field, const MapKeyComparator * key_comparator);
    // MapKeyComparator* CreateMultipleFieldsMapKeyComparator(const std::vector< std::vector< const FieldDescriptor * > > & key_field_paths);
    // void	AddIgnoreCriteria(IgnoreCriteria * ignore_criteria);
    void IgnoreField(const FieldDescriptor * field);

    void set_report_matches(bool report_matches);
    void set_report_moves(bool report_moves);
    void set_report_ignores(bool report_ignores);

    bool Compare(const google::protobuf::Message& message1, const google::protobuf::Message& message2);
    bool CompareWithFields(const Message & message1, const Message & message2, const std::vector< const FieldDescriptor * > & message1_fields, 
      const std::vector< const FieldDescriptor * > & message2_fields);
    
    %extend {
      std::string ReportDifferencesToString() {
        std::string s;
        $self->ReportDifferencesToString(&s);
        return s;
      }
    }
  };

} } } // namespace google::protobuf::util


