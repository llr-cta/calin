/*

   calin/io/sql_transceiver.hpp -- Stephen Fegan -- 2015-09-08

   Base class for reading and writing protobuf structures to SQL databases

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

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

#include <io/sql_statement.hpp>

namespace calin { namespace io { namespace sql_transceiver {

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
    KEY_INHERITED, KEY_USER_SUPPLIED, KEY_PROTO_DEFINED, KEY_PARENT_OID,
    KEY_LOOP_ID, KEY_MAP_KEY, POD };
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

class SQLTransceiver
{
 public:
  SQLTransceiver(bool write_sql_to_log = false):
      write_sql_to_log_(write_sql_to_log) { /* nothing to see here */ }
  virtual ~SQLTransceiver();

  // ===========================================================================
  //
  // Primary interface
  //
  // ===========================================================================

  virtual bool
  create_tables(const std::string& table_name,
                const google::protobuf::Descriptor* d_data,
                const google::protobuf::Descriptor* d_key = nullptr,
                const std::string& instance_desc = "");

  bool insert(const std::string& table_name,
              const google::protobuf::Message* m_data,
              const google::protobuf::Message* m_key = nullptr)
  {
    uint64_t unused_oid;
    return insert(table_name, unused_oid, m_data, m_key);
  }

  bool create_tables_and_insert(const std::string& table_name,
              const google::protobuf::Message* m_data,
              const google::protobuf::Message* m_key = nullptr,
              const std::string& instance_desc = "")
  {
    if(!this->create_tables(table_name, m_data->GetDescriptor(),
      m_key?m_key->GetDescriptor():nullptr, instance_desc))return false;
    uint64_t unused_oid;
    return this->insert(table_name, unused_oid, m_data, m_key);
  }

  bool create_tables_and_insert(const std::string& table_name, uint64_t& oid,
              const google::protobuf::Message* m_data,
              const google::protobuf::Message* m_key = nullptr,
              const std::string& instance_desc = "")
  {
    if(!this->create_tables(table_name, m_data->GetDescriptor(),
      m_key?m_key->GetDescriptor():nullptr, instance_desc))return false;
    return this->insert(table_name, oid, m_data, m_key);
  }

  virtual bool
  insert(const std::string& table_name, uint64_t& oid,
         const google::protobuf::Message* m_data,
         const google::protobuf::Message* m_key = nullptr);

  virtual bool
  retrieve_by_oid(const std::string& table_name, uint64_t oid,
                  google::protobuf::Message* m_data,
                  google::protobuf::Message* m_key = nullptr);

  virtual uint64_t count_entries_in_table(const std::string& table_name);

  virtual std::vector<uint64_t>
  retrieve_all_oids(const std::string& table_name);

  // ===========================================================================
  //
  // Utility functions
  //
  // ===========================================================================

  SQLTable*
  make_keyed_sqltable_tree(const std::string& table_name,
                           const google::protobuf::Descriptor* d_data,
                           const google::protobuf::Descriptor* d_key = nullptr,
                           bool ignore_key_option = false);

  SQLTable*
  make_sqltable_tree(const std::string& table_name,
                     const google::protobuf::Descriptor* d,
                     bool ignore_key_option = false);

  SQLTable*
  make_extkey_tree(const std::string& table_name,
                   const google::protobuf::Descriptor* d);

  void prune_empty_tables(SQLTable* t);

  void propagate_keys(SQLTable* t,
                      std::vector<const SQLTableField*> keys = { });

  std::vector<std::pair<std::string,std::string> >
  list_all_table_columns(const SQLTable* t);

  template<typename FCN> void iterate_over_tables(SQLTable* t, FCN fcn)
  {
    std::deque<SQLTable*> all_t;
    all_t.push_back(t);
    while(!all_t.empty())
    {
      SQLTable* it = all_t.front();
      all_t.pop_front();
      all_t.insert(all_t.begin(), it->sub_tables.begin(), it->sub_tables.end());
      fcn(it);
    }
  }

  template<typename FCN> void iterate_over_tables(const SQLTable* t, FCN fcn)
  {
    std::deque<const SQLTable*> all_t;
    all_t.push_back(t);
    while(!all_t.empty())
    {
      const SQLTable* it = all_t.front();
      all_t.pop_front();
      fcn(it);
      all_t.insert(all_t.begin(), it->sub_tables.begin(), it->sub_tables.end());
    }
  }

  template<typename FCN, typename DATUM>
  void iterate_over_tables_with_data(SQLTable* t, const DATUM& d0, FCN fcn)
  {
    std::deque<SQLTable*> all_t;
    std::deque<DATUM> all_d;
    all_t.push_back(t);
    all_d.push_back(d0);
    while(!all_t.empty())
    {
      SQLTable* it = all_t.front();
      DATUM d = all_d.front();
      all_t.pop_front();
      all_d.pop_front();
      all_t.insert(all_t.begin(), it->sub_tables.begin(), it->sub_tables.end());
      fcn(it, d); // Modified d here
      all_d.insert(all_d.begin(), it->sub_tables.size(), d);
    }
  }

  template<typename FCN, typename DATUM>
  void iterate_over_tables_with_data(const SQLTable* t, const DATUM& d0, FCN fcn)
  {
    std::deque<const SQLTable*> all_t;
    std::deque<DATUM> all_d;
    all_t.push_back(t);
    all_d.push_back(d0);
    while(!all_t.empty())
    {
      const SQLTable* it = all_t.front();
      DATUM d = all_d.front();
      all_t.pop_front();
      all_d.pop_front();
      all_t.insert(all_t.begin(), it->sub_tables.begin(), it->sub_tables.end());
      fcn(it, d); // Modified d here
      all_d.insert(all_d.begin(), it->sub_tables.size(), d);
    }
  }

  template<typename FCN> void iterate_over_fields(SQLTable* t, FCN fcn)
  {
    iterate_over_tables(t, [&fcn](SQLTable* t) {
        for(auto f : t->fields) { fcn(t,f); } });
  }

  template<typename FCN> void iterate_over_fields(const SQLTable* t, FCN fcn)
  {
    iterate_over_tables(t, [&fcn](const SQLTable* t) {
        for(auto f : t->fields) { fcn(t,f); } });
  }

  template<typename FCN_f, typename FCN_cd, typename DATUM> void
  iterate_over_fields_with_data(SQLTable* t, const DATUM& d0, FCN_cd fcn_cd,
                                FCN_f fcn_f)
  {
    iterate_over_tables_with_data(t, d0,
                               [&fcn_cd,&fcn_f](SQLTable* t, DATUM& d) {
        fcn_cd(t, d); // possible change in d here
        for(auto f : t->fields) { DATUM sd = d; fcn_f(t,f,sd); } });
  }

  template<typename FCN_f, typename FCN_cd, typename DATUM> void
  iterate_over_fields_with_data(const SQLTable* t, const DATUM& d0,
                                FCN_cd fcn_cd, FCN_f fcn_f)
  {
    iterate_over_tables_with_data(t, d0,
                               [&fcn_cd, &fcn_f](const SQLTable* t, DATUM& d) {
        fcn_cd(t, d); // possible change in d here
        for(auto f : t->fields) { DATUM sd = d; fcn_f(t,f,sd); } });
  }

 protected:

  bool write_sql_to_log_ = false;
  bool internal_tables_created_ = false;

  // ===========================================================================
  //
  // Private tree-related functions
  //
  // ===========================================================================

  SQLTable*
  r_make_sqltable_tree(const std::string& table_name,
                       const google::protobuf::Descriptor* d,
                       SQLTable* parent_table,
                       const google::protobuf::FieldDescriptor* parent_field_d,
                       bool ignore_key_option);

  // ===========================================================================
  //
  // Overridable member functions to create SQL strings
  //
  // ===========================================================================

  virtual std::string user_key_name(const std::string& name);
  virtual std::string sub_name(const std::string& parent_name,
                               const std::string& name);
  virtual std::string sql_table_name(const std::string& name);
  virtual std::string sql_field_name(const std::string& name);

  virtual std::string sql_type(const google::protobuf::FieldDescriptor* d);
  virtual std::string sql_oid_column_name();

  virtual std::string sql_create_table(const SQLTable* t,
                           const std::map<std::string,std::string>& parent_dict,
                                       bool if_not_exists = false);
  virtual std::string sql_create_index(const SQLTable* t,
                                       bool if_not_exists = false);
  virtual std::string sql_insert(const SQLTable* t);
  virtual std::string sql_select(const SQLTable* t, bool select_oid = false);
  virtual std::string sql_where_oid_equals();
  virtual std::string sql_where_inherited_keys_match(const SQLTable* t);

  virtual std::string sql_comment(const std::string& comment,
                                  unsigned first_line_indent = 0,
                                  unsigned multi_line_indent = 0,
                                  bool newline_before_multi_line = false);

  virtual std::string sql_count_entries(const std::string& table_name);
  virtual std::string sql_retrieve_oids(const std::string& table_name);

  // ===========================================================================
  //
  // Overridable prepare and execute statements
  //
  // ===========================================================================

  virtual SQLStatement* prepare_statement(const std::string& sql);

  virtual bool begin_transaction();
  virtual bool commit_transaction();
  virtual bool rollback_transaction();

  virtual bool r_exec_simple(SQLTable* t, bool ignore_errors);

  void set_const_data_pointers(SQLTable* t,
                               const google::protobuf::Message* m_data,
                               const google::protobuf::Message* m_key,
                               const uint64_t* parent_oid,
                               const uint64_t* loop_id);

  void bind_fields_from_data_pointers(const SQLTable* t, uint64_t loop_id,
                                      SQLStatement* stmt,
                                      bool bind_inherited_keys_only = false);

  bool field_selected(const google::protobuf::FieldDescriptor* f,
                      const std::map<const google::protobuf::OneofDescriptor*,
                                     int>& oneof_map);

  virtual bool r_exec_insert(SQLTable* t,
                             const google::protobuf::Message* m_data,
                             const google::protobuf::Message* m_key,
                             uint64_t& oid, uint64_t parent_oid,
                             uint64_t loop_id, bool ignore_errors);

  virtual bool r_exec_select(SQLTable* t,
                             google::protobuf::Message* m_data,
                             google::protobuf::Message* m_key,
                             uint64_t& loop_id,
                             bool ignore_errors);

  virtual bool finalize_statements(SQLTable* t);

  virtual void create_internal_tables();
  virtual void insert_table_description(const SQLTable* t,
                                        const std::string& instance_desc);

  std::map<std::string, std::string>
  field_dict(const google::protobuf::FieldDescriptor *d,
             const std::map<std::string, std::string>& parent_dict);
};

} } } // namespace calin::io::sql_transceiver
