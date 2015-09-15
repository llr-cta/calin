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
  SQLite3Transceiver(const std::string& filename);
  SQLite3Transceiver(sqlite3* db, bool inherit_db = false):
      SQLTransceiver(), db_(db), inherit_db_(inherit_db) { }
  virtual ~SQLite3Transceiver();

  bool create_tables(const std::string& table_name,
                     const google::protobuf::Descriptor* d_data,
                     const google::protobuf::Descriptor* d_key = nullptr,
                     const std::string& instance_desc = "",
                     bool write_sql_to_log = false) override;

 private:
  sqlite3* db_;
  bool inherit_db_;
};

} } } // namespace calin::io::sql_transceiver
