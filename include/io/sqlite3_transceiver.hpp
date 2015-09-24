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
#include "io/sqlite3_statement.hpp"

namespace calin { namespace io { namespace sql_transceiver {

class SQLite3Transceiver: public SQLTransceiver
{
 public:
  enum OpenMode { EXISTING_OR_NEW_RW, EXISTING_RW, TRUNCATE_RW, READ_ONLY };
  
  SQLite3Transceiver(const std::string& filename,
                     OpenMode open_mode = EXISTING_OR_NEW_RW,
                     bool write_sql_to_log = false);
  virtual ~SQLite3Transceiver();

 protected:

  sqlite3* db_;
  bool inherit_db_ = false;
  OpenMode open_mode_ = EXISTING_OR_NEW_RW;
  
  SQLStatement* prepare_statement(const std::string& sql) override;
  
  void create_calin_tables();
};

} } } // namespace calin::io::sql_transceiver
