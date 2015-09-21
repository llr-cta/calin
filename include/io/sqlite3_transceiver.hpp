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

 protected:

  class SQLite3Statement: public SQLTransceiver::Statement
  {
   public:
    SQLite3Statement(const std::string& sql, sqlite3* db,
                     bool make_bound_sql = false);
    virtual ~SQLite3Statement();

    unsigned num_columns() override;
    
    bool is_initialized() override;
    void error_codes(int& error_num, std::string& error_msg) override;

    void reset() override;
    bool step() override;

    bool bind_null(unsigned ifield) override;
    bool bind_int64(unsigned ifield, int64_t value) override;
    bool bind_int32(unsigned ifield, int32_t value) override;
    bool bind_int16(unsigned ifield, int16_t value) override;
    bool bind_int8(unsigned ifield, int8_t value) override;
    bool bind_uint64(unsigned ifield, uint64_t value) override;
    bool bind_uint32(unsigned ifield, uint32_t value) override;
    bool bind_uint16(unsigned ifield, uint16_t value) override;
    bool bind_uint8(unsigned ifield, uint8_t value) override;
    bool bind_float(unsigned ifield, float value) override;
    bool bind_double(unsigned ifield, double value) override;
    bool bind_bool(unsigned ifield, bool value) override;
    bool bind_string(unsigned ifield, const std::string& value) override;
    bool bind_bytes(unsigned ifield, const std::string& value) override;

   protected:
    sqlite3* db_ = nullptr;
    bool make_bound_sql_;
    sqlite3_stmt* stmt_ = nullptr;
  };

  sqlite3* db_;
  bool inherit_db_ = false;
  OpenMode open_mode_ = EXISTING_OR_NEW_RW;
  
  int execute_simple_sql(const std::string& sql, bool write_sql_to_log = false,
                         bool ignore_errors = false);

  void create_calin_tables();
};

} } } // namespace calin::io::sql_transceiver
