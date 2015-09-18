/* 

   calin/io/sqlite3_transceiver.hpp -- Stephen Fegan -- 2015-09-08

   Provides reading and writing of protobuf structures to SQLite3 databases

*/

#pragma once

#include <string>
#include <streambuf>
#include <iostream>
#include <sstream>
#include <vector>

#include <google/protobuf/descriptor.pb.h>

#include <sqlite3.h>

#include "io/sql_transceiver.hpp"

namespace calin { namespace io { namespace sql_transceiver {

class SQLite3Transceiver: public SQLTransceiver
{
 public:
  enum OpenMode { EXISTING_OR_NEW_RW, EXISTING_RW, TRUNCATE_RW, READ_ONLY };
  
  SQLite3Transceiver(const std::string& filename,
                     OpenMode open_mode = EXISTING_OR_NEW_RW);
  virtual ~SQLite3Transceiver();

  bool create_tables(const std::string& table_name,
                     const google::protobuf::Descriptor* d_data,
                     const google::protobuf::Descriptor* d_key = nullptr,
                     const std::string& instance_desc = "",
                     bool write_sql_to_log = false) override;

  bool insert(const std::string& table_name,
              const google::protobuf::Message* m_data,
              const google::protobuf::Message* m_key,
              bool write_sql_to_log) override;

 private:
  sqlite3* db_;
  bool inherit_db_ = false;
  OpenMode open_mode_ = EXISTING_OR_NEW_RW;
  
  int execute_simple_sql(const std::string& sql, bool write_sql_to_log = false,
                         bool ignore_errors = false);

  void create_calin_tables();
};

} } } // namespace calin::io::sql_transceiver
