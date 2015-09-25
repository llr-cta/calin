/* 

   calin/io/sql_transceiver.hpp -- Stephen Fegan -- 2015-09-08

   Base class for reading and writing protobuf structures to SQL databases

*/

#pragma once

#include <string>
#include <streambuf>
#include <iostream>
#include <sstream>
#include <vector>
#include <deque>

#include <google/protobuf/descriptor.pb.h>

#include "io/sql_statement.hpp"

namespace calin { namespace io { namespace sql_transceiver {

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

  virtual bool
  insert(const std::string& table_name, uint64_t& oid,
         const google::protobuf::Message* m_data,
         const google::protobuf::Message* m_key = nullptr);

  // ===========================================================================
  //
  // Utility functions
  //
  // ===========================================================================

  struct SQLTable;
  
  struct SQLTableField
  {
    enum FieldType {
      KEY_INHERITED, KEY_USER_SUPPLIED, KEY_PROTO_DEFINED, KEY_OID,
      KEY_LOOP_ID, KEY_MAP_KEY, POD };
    SQLTable*                                 table = nullptr;
    SQLTableField*                            field_origin = nullptr;
    FieldType                                 field_type;
    std::string                               field_name;
    const google::protobuf::FieldDescriptor*  field_d = nullptr;
    std::vector<const google::protobuf::FieldDescriptor*> field_d_path;
    void*                                     data = nullptr;
    const google::protobuf::OneofDescriptor*  oneof_d = nullptr;
    bool is_key() const { return field_type!=POD; }
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
  };

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

  virtual std::string sql_create_table(const SQLTable* t,
                                       bool if_not_exists = false);
  virtual std::string sql_insert(const SQLTable* t);

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
  
  virtual bool r_exec_insert(SQLTable* t,
                             const google::protobuf::Message* m_data,
                             const google::protobuf::Message* m_key,
                             uint64_t& oid, uint64_t parent_oid,
                             uint64_t loop_id, bool ignore_errors);

  virtual bool finalize_statements(SQLTable* t);

  virtual void create_internal_tables();
  virtual void insert_table_description(const SQLTable* t,
                                        const std::string& instance_desc);
};

} } } // namespace calin::io::sql_transceiver
