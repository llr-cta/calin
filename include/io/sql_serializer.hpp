/*

   calin/io/sql_serializer.hpp -- Stephen Fegan -- 2020-03-31

   Base class for reading and writing protobuf structures to SQL databases

   Based on : calin/io/sql_transceiver.hpp -- Stephen Fegan -- 2015-09-08

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

#pragma once

#include <string>
#include <streambuf>
#include <iostream>
#include <sstream>
#include <vector>
#include <deque>

#include <google/protobuf/message.h>
#include <google/protobuf/descriptor.h>

//#include <io/sql_statement.hpp>

namespace calin { namespace io { namespace sql_serializer {

struct SQLStatement {};

struct SQLTable;

enum DataPointerType {
  PT_NULL, PT_CONST_UINT64, PT_CONST_MESSAGE, PT_UINT64, PT_MESSAGE  };

union DataPointer {
  DataPointer(): p_void(nullptr) { /* nothing to see here */ }
  void* p_void;
  const uint64_t* p_const_uint64;
  uint64_t* p_uint64;
  const google::protobuf::Message* p_const_message;
  google::protobuf::Message* p_message;
};

struct SQLTableField
{
  enum FieldType {
    KEY_INHERITED, KEY_PARENT_OID, KEY_LOOP_ID, KEY_MAP_KEY, POD };
  SQLTable*                                 table = nullptr;
  SQLTableField*                            field_origin = nullptr;
  FieldType                                 field_type;
  std::string                               field_name;
  const google::protobuf::FieldDescriptor*  field_d = nullptr;
  std::vector<const google::protobuf::FieldDescriptor*> field_d_path;

  DataPointerType                           data_type = PT_NULL;
  DataPointer                               data;
  const google::protobuf::OneofDescriptor*  oneof_d = nullptr;

  bool is_key() const { return field_type!=POD; }
  bool is_inherited() const { return field_type==KEY_INHERITED or
        field_type==KEY_PARENT_OID; }

  void set_data_null() { data_type = PT_NULL; data.p_void = nullptr; }
  void set_data_const_uint64(const uint64_t* p) {
    data_type = PT_CONST_UINT64; data.p_const_uint64 = p; }
  void set_data_uint64(uint64_t* p) {
    data_type = PT_UINT64; data.p_uint64 = p; }
  void set_data_const_message(const google::protobuf::Message* p) {
    data_type = PT_CONST_MESSAGE; data.p_const_message = p; }
  void set_data_message(google::protobuf::Message* p) {
    data_type = PT_MESSAGE; data.p_message = p; }

  bool is_data_null() const { return data_type == PT_NULL; }
  const uint64_t* data_const_uint64() const {
    if(data_type == PT_CONST_UINT64) return data.p_const_uint64;
    else if(data_type == PT_UINT64) return data.p_uint64;
    else return nullptr; }
  uint64_t* data_uint64() const {
    return data_type == PT_UINT64 ? data.p_uint64 : nullptr; }
  const google::protobuf::Message* data_const_message() const {
    if(data_type == PT_CONST_MESSAGE) return data.p_const_message;
    else if(data_type == PT_MESSAGE) return data.p_message;
    else return nullptr; }
  google::protobuf::Message* data_message() const {
    return data_type == PT_MESSAGE ? data.p_message : nullptr; }
};

#ifdef SWIG
} } }
%template(VectorSQLTablePtr) std::vector<calin::io::sql_serializer::SQLTable*>;
%template(VectorSQLTableFieldPtr) std::vector<calin::io::sql_serializer::SQLTableField*>;
namespace calin { namespace io { namespace sql_serializer {
#endif

struct SQLTable
{
  ~SQLTable() {
    for(auto ifield : fields)delete ifield;
    for(auto itable : sub_tables)delete itable;
    delete stmt;
  }

  std::string                               table_name;
  SQLTable*                                 parent_table = nullptr;
  const google::protobuf::FieldDescriptor*  parent_field_d = nullptr;
  std::vector<const google::protobuf::FieldDescriptor*> parent_field_d_path;
  std::vector<SQLTableField*>               fields;
  std::vector<SQLTable*>                    sub_tables;
  SQLStatement*                             stmt = nullptr;
  bool children_need_oid() const {
    for(auto t : sub_tables) { for(auto f : t->fields)
        if(f->field_type == SQLTableField::KEY_PARENT_OID)return true; }
    return false; }
};

class SQLSerializer
{
 public:
  SQLSerializer(bool write_sql_to_log = false):
      write_sql_to_log_(write_sql_to_log) { /* nothing to see here */ }
  virtual ~SQLSerializer();


  // ===========================================================================
  //
  // Tree-related functions - this should be brought private
  //
  // ===========================================================================

  SQLTable*
  make_sqltable_tree(const std::string& table_name,
                     const google::protobuf::Descriptor* d);

  void propagate_keys(SQLTable* t);

protected:

  bool write_sql_to_log_ = false;

  // ===========================================================================
  //
  // Private tree-related functions
  //
  // ===========================================================================

  SQLTable*
  r_make_sqltable_tree(const std::string& table_name,
                       const google::protobuf::Descriptor* d,
                       SQLTable* parent_table,
                       const google::protobuf::FieldDescriptor* parent_field_d);

  void r_propagate_keys(SQLTable* t, std::vector<const SQLTableField*> keys);

  // ===========================================================================
  //
  // Overridable member functions to create SQL strings
  //
  // ===========================================================================

  virtual std::string sub_name(const std::string& parent_name, const std::string& name);
  virtual std::string sql_oid_column_name();
};

} } } // namespace calin::io::sql_serializer
