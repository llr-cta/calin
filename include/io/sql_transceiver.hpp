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

namespace calin { namespace io { namespace sql_transceiver {

class SQLTransceiver
{
 protected:
  class Statement;
  
 public:
  virtual ~SQLTransceiver();

  // --------------------------------------------------------------------------
  // Primary interface
  // --------------------------------------------------------------------------

  virtual bool
  create_tables(const std::string& table_name,
                const google::protobuf::Descriptor* d_data,
                const google::protobuf::Descriptor* d_key = nullptr,
                const std::string& instance_desc = "",
                bool write_sql_to_log = false);

  virtual bool
  insert(const std::string& table_name,
         const google::protobuf::Message* m_data,
         const google::protobuf::Message* m_key = nullptr,
         bool write_sql_to_log = false);

  // --------------------------------------------------------------------------
  // Utility functions
  // --------------------------------------------------------------------------

  struct SQLTable;
  
  struct SQLTableField
  {
    enum FieldType {
      KEY_INHERITED, KEY_USER_SUPPLIED, KEY_PROTO_DEFINED, KEY_OID,
      KEY_LOOP_ID, KEY_MAP_KEY, POD };
    SQLTable* table = nullptr;      // Table the field belongs to
    SQLTableField* field_origin = nullptr; // Ptr to first occurance of field (if key)
    FieldType field_type;           // Type of field
    std::string field_name;         // Full name of field
    const google::protobuf::FieldDescriptor* field_d;
    std::vector<const google::protobuf::FieldDescriptor*> field_d_path;
    void* data = nullptr;           // Pointer to data (message or index)
    //google::protobuf::FieldDescriptor* message;
    //uint64_t* rowid_or_vecindex;
    bool is_key() const { return field_type!=POD; }
  };

  struct SQLTable
  {
    ~SQLTable() { for(auto ifield : fields)delete ifield;
      for(auto itable : sub_tables)delete itable; }

    std::string table_name;
    SQLTable* parent_table = nullptr;
    const google::protobuf::FieldDescriptor* parent_field_d = nullptr;
    std::vector<const google::protobuf::FieldDescriptor*> parent_field_d_path;
    std::vector<SQLTableField*> fields;
    std::vector<SQLTable*> sub_tables;
    Statement* stmt = nullptr;
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
  
  static bool is_single_table(const google::protobuf::Descriptor* d);

  std::vector<std::pair<std::string,std::string> >
  list_all_table_columns(SQLTable* t);
  
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
  
 protected:

  class Statement
  {
   public:
    Statement(const std::string& sql);
    virtual ~Statement();

    const std::string& sql() const { return sql_; }
    std::string bound_sql() const;
    
    virtual unsigned num_columns();

    virtual bool is_initialized();
    virtual void error_codes(int& error_num, std::string& error_msg);
    
    virtual void reset();
    virtual bool step();

    bool bind_field(unsigned ifield, const google::protobuf::Message* m,
                    const google::protobuf::FieldDescriptor* d);
    bool bind_repeated_field(unsigned ifield, uint64_t loop_id, 
                             const google::protobuf::Message* m,
                             const google::protobuf::FieldDescriptor* d);

    virtual bool bind_null(unsigned ifield);
    virtual bool bind_int64(unsigned ifield, int64_t value);
    virtual bool bind_int32(unsigned ifield, int32_t value);
    virtual bool bind_int16(unsigned ifield, int16_t value);
    virtual bool bind_int8(unsigned ifield, int8_t value);
    virtual bool bind_uint64(unsigned ifield, uint64_t value);
    virtual bool bind_uint32(unsigned ifield, uint32_t value);
    virtual bool bind_uint16(unsigned ifield, uint16_t value);
    virtual bool bind_uint8(unsigned ifield, uint8_t value);
    virtual bool bind_float(unsigned ifield, float value);
    virtual bool bind_double(unsigned ifield, double value);
    virtual bool bind_bool(unsigned ifield, bool value);
    virtual bool bind_string(unsigned ifield, const std::string& value);
    virtual bool bind_bytes(unsigned ifield, const std::string& value);

   protected:
    std::string sql_;
    std::vector<std::string> bound_values_;
  };

  SQLTable*
  r_make_sqltable_tree(const std::string& table_name,
                       const google::protobuf::Descriptor* d,
                       SQLTable* parent_table,
                       const google::protobuf::FieldDescriptor* parent_field_d,
                       bool ignore_key_option);
  
  virtual std::string user_key_name(const std::string& name);
  virtual std::string sub_name(const std::string& parent_name,
                               const std::string& name);
  virtual std::string sql_table_name(const std::string& name);
  virtual std::string sql_field_name(const std::string& name);

  virtual std::string sql_type(const google::protobuf::FieldDescriptor* d);

  virtual std::string sql_create_table(const SQLTable* t,
                                       bool if_not_exists = false);
  virtual std::string sql_insert(const SQLTable* t);


  virtual bool prepare_insert_statements(SQLTable* t);
  virtual bool insert_tables(SQLTable* t,
                             const google::protobuf::Message* m_data,
                             const google::protobuf::Message* m_key,
                             uint64_t parent_oid, uint64_t loop_id,
                             bool write_sql_to_log);
  virtual bool execute_one_insert_statement(SQLTable* t, uint64_t& oid);
  virtual bool finalize_insert_statements(SQLTable* t);
};

} } } // namespace calin::io::sql_transceiver
