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

#include <google/protobuf/descriptor.pb.h>

namespace calin { namespace io { namespace sql_transceiver {

class SQLTransceiver
{
 public:
  virtual ~SQLTransceiver();

  static bool is_single_table(const google::protobuf::Descriptor* d);
  
  virtual bool
  create_tables(const std::string& table_name,
                const google::protobuf::Descriptor* d_data,
                const google::protobuf::Descriptor* d_key = nullptr,
                const std::string& instance_desc = "",
                bool write_sql_to_log = false, bool commit = true) = 0;

  static std::vector<std::pair<std::string,std::string>>
      list_all_table_columns(const std::string& table_name,
                             const google::protobuf::Descriptor* d);
  
 protected:

  struct SQLTable;
  
  struct SQLTableField
  {
    enum FieldType {
      KEY_USER_SUPPLIED, KEY_PROTO_DEFINED, KEY_OID,
      KEY_LOOP_ID, KEY_MAP_KEY, POD };
    SQLTable* table;
    FieldType field_type;
    std::string field_name;
    const google::protobuf::FieldDescriptor* field_d;
    std::vector<const google::protobuf::FieldDescriptor*> field_d_path;
    //google::protobuf::FieldDescriptor* message;
    //uint64_t* rowid_or_vecindex;
    bool is_key() const { return field_type!=POD; }
  };

  struct SQLTable
  {
    ~SQLTable() { for(auto ifield : fields)delete ifield;
      for(auto itable : sub_tables)delete itable; }

    std::string table_name;
    SQLTable* parent_table;
    const google::protobuf::FieldDescriptor* parent_field_d;
    std::vector<const google::protobuf::FieldDescriptor*> parent_field_d_path;
    std::vector<SQLTableField*> fields;
    std::vector<SQLTable*> sub_tables;
  };

  static std::vector<std::pair<std::string,std::string> >
  extract_all_table_columns(SQLTable* t);
  
  static SQLTable*
  make_table(const std::string& table_name,
             const google::protobuf::Descriptor* d,
             SQLTable* parent_table = nullptr,
             const google::protobuf::FieldDescriptor* parent_field_d
             = nullptr);
  
  static std::string key_name(const std::string& name);
  static std::string sub_name(const std::string& parent_name,
                              const std::string& name);
  
};

} } } // namespace calin::io::sql_transceiver
