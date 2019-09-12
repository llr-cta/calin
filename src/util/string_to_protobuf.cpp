/*

   calin/util/string_to_protobuf.cpp -- Stephen Fegan -- 2016-08-10

   String to and from protobuf

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

   Based on original, copyright 2006, Stephen Fegan, see notice below

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <util/string.hpp>
#include <util/string_to_protobuf.hpp>

using namespace calin::util::string;

// =============================================================================
// =============================================================================
//
// String to Protobuf
//
// =============================================================================
// =============================================================================

namespace {

#define SET_SCALAR(HandleType, type)              \
{                                                 \
  type x;                                         \
  good_conversion = from_string(s, x);            \
  if(good_conversion) {                           \
    for(const auto* if_path : f_path) {           \
      m = r->MutableMessage(m, if_path);          \
      assert(m);                                  \
      r = m->GetReflection();                     \
      assert(r);                                  \
    }                                             \
    r->HandleType(m, f, x);                       \
  }                                               \
}

bool set_scalar_protobuf_field(const std::string& s,
  google::protobuf::Message* m, const google::protobuf::FieldDescriptor* f,
  const std::vector<const google::protobuf::FieldDescriptor*>& f_path = {})
{
  bool good_conversion = false;
  const google::protobuf::Reflection* r = m->GetReflection();
  assert(r);
  switch(f->type())
  {
  case google::protobuf::FieldDescriptor::TYPE_DOUBLE:
    SET_SCALAR(SetDouble, double);
    break;
  case google::protobuf::FieldDescriptor::TYPE_FLOAT:
    SET_SCALAR(SetFloat, float);
    break;
  case google::protobuf::FieldDescriptor::TYPE_SFIXED64: // fallthrough
  case google::protobuf::FieldDescriptor::TYPE_SINT64:   // fallthrough
  case google::protobuf::FieldDescriptor::TYPE_INT64:
    SET_SCALAR(SetInt64, int64_t);
    break;
  case google::protobuf::FieldDescriptor::TYPE_FIXED64:  // fallthrough
  case google::protobuf::FieldDescriptor::TYPE_UINT64:
    SET_SCALAR(SetUInt64, uint64_t);
    break;
  case google::protobuf::FieldDescriptor::TYPE_SFIXED32: // fallthrough
  case google::protobuf::FieldDescriptor::TYPE_SINT32:   // fallthrough
  case google::protobuf::FieldDescriptor::TYPE_INT32:
    SET_SCALAR(SetInt32, int32_t);
    break;
  case google::protobuf::FieldDescriptor::TYPE_FIXED32:
  case google::protobuf::FieldDescriptor::TYPE_UINT32:
    SET_SCALAR(SetUInt32, uint32_t);
    break;
  case google::protobuf::FieldDescriptor::TYPE_BOOL:
    SET_SCALAR(SetBool, bool);
    break;
  case google::protobuf::FieldDescriptor::TYPE_STRING:
    for(const auto* if_path : f_path) {
      m = r->MutableMessage(m, if_path);
      assert(m);
      r = m->GetReflection();
      assert(r);
    }
    r->SetString(m, f, s);
    good_conversion = true;
    break;
  case google::protobuf::FieldDescriptor::TYPE_BYTES:
    for(const auto* if_path : f_path) {
      m = r->MutableMessage(m, if_path);
      assert(m);
      r = m->GetReflection();
      assert(r);
    }
    r->SetString(m, f, s);
    good_conversion = true;
    break;
  case google::protobuf::FieldDescriptor::TYPE_ENUM:
    {
      const google::protobuf::EnumDescriptor* ed = f->enum_type();
      assert(ed);
      const google::protobuf::EnumValueDescriptor* evd = ed->FindValueByName(s);
      if(evd) {
        good_conversion = true;
        for(const auto* if_path : f_path) {
          m = r->MutableMessage(m, if_path);
          assert(m);
          r = m->GetReflection();
          assert(r);
        }
        r->SetEnum(m, f, evd);
      }
    }
  case google::protobuf::FieldDescriptor::TYPE_MESSAGE:
  case google::protobuf::FieldDescriptor::TYPE_GROUP:
  default:
    assert(0);
    break;
  }
  return good_conversion;
}

bool add_scalar_protobuf_field(const std::string& s,
  google::protobuf::Message* m, const google::protobuf::FieldDescriptor* f,
  const std::vector<const google::protobuf::FieldDescriptor*>& f_path = {})
{
  bool good_conversion = false;
  const google::protobuf::Reflection* r = m->GetReflection();
  assert(r);
  switch(f->type())
  {
  case google::protobuf::FieldDescriptor::TYPE_DOUBLE:
    SET_SCALAR(AddDouble, double);
    break;
  case google::protobuf::FieldDescriptor::TYPE_FLOAT:
    SET_SCALAR(AddFloat, float);
    break;
  case google::protobuf::FieldDescriptor::TYPE_SFIXED64: // fallthrough
  case google::protobuf::FieldDescriptor::TYPE_SINT64:   // fallthrough
  case google::protobuf::FieldDescriptor::TYPE_INT64:
    SET_SCALAR(AddInt64, int64_t);
    break;
  case google::protobuf::FieldDescriptor::TYPE_FIXED64:  // fallthrough
  case google::protobuf::FieldDescriptor::TYPE_UINT64:
    SET_SCALAR(AddUInt64, uint64_t);
    break;
  case google::protobuf::FieldDescriptor::TYPE_SFIXED32: // fallthrough
  case google::protobuf::FieldDescriptor::TYPE_SINT32:   // fallthrough
  case google::protobuf::FieldDescriptor::TYPE_INT32:
    SET_SCALAR(AddInt32, int32_t);
    break;
  case google::protobuf::FieldDescriptor::TYPE_FIXED32:
  case google::protobuf::FieldDescriptor::TYPE_UINT32:
    SET_SCALAR(AddUInt32, uint32_t);
    break;
  case google::protobuf::FieldDescriptor::TYPE_BOOL:
    SET_SCALAR(AddBool, bool);
    break;
  case google::protobuf::FieldDescriptor::TYPE_STRING:
    for(const auto* if_path : f_path) {
      m = r->MutableMessage(m, if_path);
      assert(m);
      r = m->GetReflection();
      assert(r);
    }
    r->AddString(m, f, s);
    good_conversion = true;
    break;
  case google::protobuf::FieldDescriptor::TYPE_BYTES:
    for(const auto* if_path : f_path) {
      m = r->MutableMessage(m, if_path);
      assert(m);
      r = m->GetReflection();
      assert(r);
    }
    r->AddString(m, f, s);
    good_conversion = true;
    break;
  case google::protobuf::FieldDescriptor::TYPE_ENUM:
    {
      const google::protobuf::EnumDescriptor* ed = f->enum_type();
      assert(ed);
      const google::protobuf::EnumValueDescriptor* evd = ed->FindValueByName(s);
      if(evd) {
        good_conversion = true;
        for(const auto* if_path : f_path) {
          m = r->MutableMessage(m, if_path);
          assert(m);
          r = m->GetReflection();
          assert(r);
        }
        r->AddEnum(m, f, evd);
      }
    }
  case google::protobuf::FieldDescriptor::TYPE_MESSAGE:  // fallthrough to assert(0)
  case google::protobuf::FieldDescriptor::TYPE_GROUP:    // fallthrough to assert(0)
  default:
    break;
  }
  return good_conversion;
}

} // anonmymous namespace

bool calin::util::string_to_protobuf::
string_to_protobuf_field(const std::string& s,
  google::protobuf::Message* m, const google::protobuf::FieldDescriptor* f,
  const std::vector<const google::protobuf::FieldDescriptor*>& f_path)
{
  if(f->is_repeated())
  {
    const google::protobuf::Reflection* r = m->GetReflection();
    for(const auto* if_path : f_path) {
      m = r->MutableMessage(m, if_path);
      assert(m);
      r = m->GetReflection();
      assert(r);
    }
    r->ClearField(m, f);
    bool its_all_good = true;
    std::vector<std::string> bits = split(s, ',');
    for(const auto& ibit : bits) {
      its_all_good &= add_scalar_protobuf_field(ibit,m,f);
      if(!its_all_good)break;
    }
    return its_all_good;
  }
  else return set_scalar_protobuf_field(s, m, f, f_path);
}

// =============================================================================
// =============================================================================
//
// Protobuf to String
//
// =============================================================================
// =============================================================================

namespace {

#define GET_SCALAR(HandleType, type)              \
{                                                 \
  for(const auto* if_path : f_path) {             \
    m = &r->GetMessage(*m, if_path);              \
    assert(m);                                    \
    r = m->GetReflection();                       \
    assert(r);                                    \
  }                                               \
  type x = r->HandleType(*m, f);                  \
  s = to_string(x);                               \
}

#define GET_VECTOR(HandleType, type)              \
{                                                 \
  for(const auto* if_path : f_path) {             \
    m = &r->GetMessage(*m, if_path);              \
    assert(m);                                    \
    r = m->GetReflection();                       \
    assert(r);                                    \
  }                                               \
  type x = r->HandleType(*m, f, ielement);        \
  s = to_string(x);                               \
}

void get_scalar_protobuf_field(std::string& s,
  const google::protobuf::Message* m,
  const google::protobuf::FieldDescriptor* f,
  const std::vector<const google::protobuf::FieldDescriptor*>& f_path = {})
{
  const google::protobuf::Reflection* r = m->GetReflection();
  assert(r);
  switch(f->type())
  {
  case google::protobuf::FieldDescriptor::TYPE_DOUBLE:
    GET_SCALAR(GetDouble, double);
    break;
  case google::protobuf::FieldDescriptor::TYPE_FLOAT:
    GET_SCALAR(GetFloat, float);
    break;
  case google::protobuf::FieldDescriptor::TYPE_SFIXED64: // fallthrough
  case google::protobuf::FieldDescriptor::TYPE_SINT64:   // fallthrough
  case google::protobuf::FieldDescriptor::TYPE_INT64:
    GET_SCALAR(GetInt64, int64_t);
    break;
  case google::protobuf::FieldDescriptor::TYPE_FIXED64:  // fallthrough
  case google::protobuf::FieldDescriptor::TYPE_UINT64:
    GET_SCALAR(GetUInt64, uint64_t);
    break;
  case google::protobuf::FieldDescriptor::TYPE_SFIXED32: // fallthrough
  case google::protobuf::FieldDescriptor::TYPE_SINT32:   // fallthrough
  case google::protobuf::FieldDescriptor::TYPE_INT32:
    GET_SCALAR(GetInt32, int32_t);
    break;
  case google::protobuf::FieldDescriptor::TYPE_FIXED32:
  case google::protobuf::FieldDescriptor::TYPE_UINT32:
    GET_SCALAR(GetUInt32, uint32_t);
    break;
  case google::protobuf::FieldDescriptor::TYPE_BOOL:
    GET_SCALAR(GetBool, bool);
    break;
  case google::protobuf::FieldDescriptor::TYPE_STRING:
    for(const auto* if_path : f_path) {
      m = &r->GetMessage(*m, if_path);
      assert(m);
      r = m->GetReflection();
      assert(r);
    }
    s = r->GetString(*m, f);
    break;
  case google::protobuf::FieldDescriptor::TYPE_BYTES:
    for(const auto* if_path : f_path) {
      m = &r->GetMessage(*m, if_path);
      assert(m);
      r = m->GetReflection();
      assert(r);
    }
    r->GetString(*m, f);
    break;
  case google::protobuf::FieldDescriptor::TYPE_ENUM:
    {
      for(const auto* if_path : f_path) {
        m = &r->GetMessage(*m, if_path);
        assert(m);
        r = m->GetReflection();
        assert(r);
      }
      const google::protobuf::EnumValueDescriptor* evd = r->GetEnum(*m, f);
      assert(evd);
      s = evd->name();
    }
  case google::protobuf::FieldDescriptor::TYPE_MESSAGE:
  case google::protobuf::FieldDescriptor::TYPE_GROUP:
  default:
    assert(0);
    break;
  }
}

void get_vector_protobuf_field(std::string& s,
  const google::protobuf::Message* m,
  const google::protobuf::FieldDescriptor* f, int ielement,
  const std::vector<const google::protobuf::FieldDescriptor*>& f_path = {})
{
  const google::protobuf::Reflection* r = m->GetReflection();
  assert(r);
  switch(f->type())
  {
  case google::protobuf::FieldDescriptor::TYPE_DOUBLE:
    GET_VECTOR(GetRepeatedDouble, double);
    break;
  case google::protobuf::FieldDescriptor::TYPE_FLOAT:
    GET_VECTOR(GetRepeatedFloat, float);
    break;
  case google::protobuf::FieldDescriptor::TYPE_SFIXED64: // fallthrough
  case google::protobuf::FieldDescriptor::TYPE_SINT64:   // fallthrough
  case google::protobuf::FieldDescriptor::TYPE_INT64:
    GET_VECTOR(GetRepeatedInt64, int64_t);
    break;
  case google::protobuf::FieldDescriptor::TYPE_FIXED64:  // fallthrough
  case google::protobuf::FieldDescriptor::TYPE_UINT64:
    GET_VECTOR(GetRepeatedUInt64, uint64_t);
    break;
  case google::protobuf::FieldDescriptor::TYPE_SFIXED32: // fallthrough
  case google::protobuf::FieldDescriptor::TYPE_SINT32:   // fallthrough
  case google::protobuf::FieldDescriptor::TYPE_INT32:
    GET_VECTOR(GetRepeatedInt32, int32_t);
    break;
  case google::protobuf::FieldDescriptor::TYPE_FIXED32:
  case google::protobuf::FieldDescriptor::TYPE_UINT32:
    GET_VECTOR(GetRepeatedUInt32, uint32_t);
    break;
  case google::protobuf::FieldDescriptor::TYPE_BOOL:
    GET_VECTOR(GetRepeatedBool, bool);
    break;
  case google::protobuf::FieldDescriptor::TYPE_STRING:
    for(const auto* if_path : f_path) {
      m = &r->GetMessage(*m, if_path);
      assert(m);
      r = m->GetReflection();
      assert(r);
    }
    s = r->GetRepeatedString(*m, f, ielement);
    break;
  case google::protobuf::FieldDescriptor::TYPE_BYTES:
    for(const auto* if_path : f_path) {
      m = &r->GetMessage(*m, if_path);
      assert(m);
      r = m->GetReflection();
      assert(r);
    }
    r->GetRepeatedString(*m, f, ielement);
    break;
  case google::protobuf::FieldDescriptor::TYPE_ENUM:
    {
      for(const auto* if_path : f_path) {
        m = &r->GetMessage(*m, if_path);
        assert(m);
        r = m->GetReflection();
        assert(r);
      }
      const google::protobuf::EnumValueDescriptor* evd =
        r->GetRepeatedEnum(*m, f, ielement);
      assert(evd);
      s = evd->name();
    }
  case google::protobuf::FieldDescriptor::TYPE_MESSAGE:
  case google::protobuf::FieldDescriptor::TYPE_GROUP:
  default:
    assert(0);
    break;
  }
}

} // anonmymous namespace

void calin::util::string_to_protobuf::
protobuf_field_to_string(std::string& s,
  const google::protobuf::Message* m,
  const google::protobuf::FieldDescriptor* f,
  const std::vector<const google::protobuf::FieldDescriptor*>& f_path)
{
  if(f->is_repeated())
  {
    const google::protobuf::Reflection* r = m->GetReflection();
    for(const auto* if_path : f_path) {
      m = &r->GetMessage(*m, if_path);
      assert(m);
      r = m->GetReflection();
      assert(r);
    }
    s.clear();
    int nelement = r->FieldSize(*m, f);
    for(int ielement=0; ielement<nelement; ielement++)
    {
      if(ielement)s += ",";
      std::string is;
      get_vector_protobuf_field(is, m, f, ielement);
      s += is;
    }
  }
  else return get_scalar_protobuf_field(s, m, f, f_path);
}
