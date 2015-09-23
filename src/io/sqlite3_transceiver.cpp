/* 

   calin/io/sqlite3_transceiver.cpp -- Stephen Fegan -- 2015-09-08

   Provides reading and writing of protobuf structures to SQLite3 databases

*/

#include <unistd.h>
#include <cstring>

#include "io/log.hpp"
#include "io/sqlite3_transceiver.hpp"

using namespace calin::io::log;
using namespace calin::io::sql_transceiver;

SQLite3Transceiver::
SQLite3Transceiver(const std::string& filename, OpenMode open_mode):
    SQLTransceiver(), db_(), inherit_db_(), open_mode_(open_mode)
{
  if(open_mode == TRUNCATE_RW and filename.front() != ':' and
     filename.back() != ':')unlink(filename.c_str());
  int flags { 0 };
  switch(open_mode) {
    case EXISTING_OR_NEW_RW:
    case TRUNCATE_RW:
      flags = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE; break;
    case EXISTING_RW:
      flags = SQLITE_OPEN_READWRITE; break;
    case READ_ONLY:
      flags = SQLITE_OPEN_READONLY; break;
  }
  sqlite3_open_v2(filename.c_str(), &db_, flags, NULL);
}

SQLite3Transceiver::~SQLite3Transceiver()
{
  if(inherit_db_)sqlite3_close(db_);
}

bool SQLite3Transceiver::
create_tables(const std::string& table_name,
              const google::protobuf::Descriptor* d_data,
              const google::protobuf::Descriptor* d_key,
              const std::string& instance_desc,
              bool write_sql_to_log)
{
  SQLTable* t = make_keyed_sqltable_tree(table_name, d_data, d_key, false);
  try {
    iterate_over_tables(t, [this,write_sql_to_log](const SQLTable* it) {
        execute_simple_sql(sql_create_table(it), write_sql_to_log); });
  } catch(...) {
    delete t;
    throw;
  }
  delete t;
  return true;
}

#if 0
bool SQLite3Transceiver::
insert(const std::string& table_name,
       const google::protobuf::Message* m_data,
       const google::protobuf::Message* m_key, bool write_sql_to_log)
{
  
}
#endif

int SQLite3Transceiver::
execute_simple_sql(const std::string& sql, bool write_sql_to_log,
                   bool ignore_errors)
{
  sqlite3_stmt* stmt { nullptr };
  int retcode { 0 };  

  retcode = sqlite3_prepare_v2(db_, sql.c_str(), -1, &stmt, nullptr);
  if(retcode != SQLITE_OK)
  {
    if(ignore_errors)return retcode;
    std::ostringstream str;
    str << "SQLITE error preparing statement: " << sqlite3_errmsg(db_)
        << " (" << retcode << " / " << sqlite3_extended_errcode(db_) << ")\n"
        << "Statement:\n" << sql;
    sqlite3_finalize(stmt);
    throw std::runtime_error(str.str());
  }
  
  retcode = sqlite3_step(stmt);
  if(retcode == SQLITE_ERROR)
  {
    if(ignore_errors)return retcode;
    std::ostringstream str;
    str << "SQLITE error executing statement: " << sqlite3_errmsg(db_)
        << " (" << retcode << " / " << sqlite3_extended_errcode(db_) << ")\n"
        << "Statement:\n" << sql;
    sqlite3_finalize(stmt);
    throw std::runtime_error(str.str());
  }
  else if(retcode == SQLITE_MISUSE)
  {
    std::ostringstream str;
    str << "SQLITE misuse executing statement:\n" << sql;
    LOG(ERROR) << str.str();
    sqlite3_finalize(stmt);
    throw std::runtime_error(str.str());
  }

  sqlite3_finalize(stmt);      

  if(write_sql_to_log)
    LOG(INFO) << sql << '\n';

  return SQLITE_DONE;
}

bool SQLite3Transceiver::
prepare_insert_statements(SQLTable* t, bool write_sql_to_log)
{
  iterate_over_tables(t, [this,write_sql_to_log](SQLTable* t) {
      t->stmt = new SQLite3Statement(sql_insert(t), db_, write_sql_to_log); });
  return true;
}


// ============================================================================
// ============================================================================
//
// SQLite3Transceiver::SQLite3Statement
//
// ============================================================================
// ============================================================================

SQLite3Transceiver::SQLite3Statement::
SQLite3Statement(const std::string& sql, sqlite3* db, bool make_bound_sql):
    SQLTransceiver::Statement(sql, make_bound_sql),
    db_(db), stmt_(nullptr), make_bound_sql_(make_bound_sql)
{
  sqlite3_prepare_v2(db_, sql.c_str(), -1, &stmt_, nullptr);
}

SQLite3Transceiver::SQLite3Statement::~SQLite3Statement()
{
  sqlite3_finalize(stmt_);
}

unsigned SQLite3Transceiver::SQLite3Statement::num_columns()
{
  return sqlite3_column_count(stmt_);
}

bool SQLite3Transceiver::SQLite3Statement::is_initialized()
{
  return stmt_ != nullptr;
}

int SQLite3Transceiver::SQLite3Statement::error_code()
{
  return sqlite3_errcode(db_);
}

std::string SQLite3Transceiver::SQLite3Statement::error_message()
{
  return sqlite3_errmsg(db_);
}

SQLTransceiver::Statement::StepStatus
SQLite3Transceiver::SQLite3Statement::step()
{
  int retcode = sqlite3_step(stmt_);
  if(retcode == SQLITE_DONE)return Statement::OK_NO_DATA;
  else if(retcode == SQLITE_ROW)return Statement::OK_HAS_DATA;
  else return Statement::ERROR;
}

uint64_t SQLite3Transceiver::SQLite3Statement::get_oid()
{
  return sqlite3_last_insert_rowid(db_);
}

void SQLite3Transceiver::SQLite3Statement::reset()
{
  sqlite3_reset(stmt_);
  sqlite3_clear_bindings(stmt_);
  if(make_bound_sql_)SQLTransceiver::Statement::reset();
}

bool SQLite3Transceiver::SQLite3Statement::bind_null(unsigned ifield)
{
  if(make_bound_sql_)SQLTransceiver::Statement::bind_null(ifield);
  return sqlite3_bind_null(stmt_, ifield+1) == SQLITE_OK;
}

bool SQLite3Transceiver::SQLite3Statement::
bind_int64(unsigned ifield, int64_t value)
{
  if(make_bound_sql_)SQLTransceiver::Statement::bind_int64(ifield, value);
  return sqlite3_bind_int64(stmt_, ifield+1, value) == SQLITE_OK;
}

bool SQLite3Transceiver::SQLite3Statement::
bind_int32(unsigned ifield, int32_t value)
{
  if(make_bound_sql_)SQLTransceiver::Statement::bind_int32(ifield, value);
  return sqlite3_bind_int(stmt_, ifield+1, value) == SQLITE_OK;
}

bool SQLite3Transceiver::SQLite3Statement::
bind_int16(unsigned ifield, int16_t value)
{
  if(make_bound_sql_)SQLTransceiver::Statement::bind_int16(ifield, value);
  return sqlite3_bind_int(stmt_, ifield+1, value) == SQLITE_OK;
}

bool SQLite3Transceiver::SQLite3Statement::
bind_int8(unsigned ifield, int8_t value)
{
  if(make_bound_sql_)SQLTransceiver::Statement::bind_int8(ifield, value);
  return sqlite3_bind_int(stmt_, ifield+1, value) == SQLITE_OK;
}

bool SQLite3Transceiver::SQLite3Statement::
bind_uint64(unsigned ifield, uint64_t value)
{
  if(make_bound_sql_)SQLTransceiver::Statement::bind_uint64(ifield, value);
  return sqlite3_bind_int64(stmt_, ifield+1, value) == SQLITE_OK;
}

bool SQLite3Transceiver::SQLite3Statement::
bind_uint32(unsigned ifield, uint32_t value)
{
  if(make_bound_sql_)SQLTransceiver::Statement::bind_uint32(ifield, value);
  return sqlite3_bind_int(stmt_, ifield+1, value) == SQLITE_OK;
}

bool SQLite3Transceiver::SQLite3Statement::
bind_uint16(unsigned ifield, uint16_t value)
{
  if(make_bound_sql_)SQLTransceiver::Statement::bind_uint16(ifield, value);
  return sqlite3_bind_int(stmt_, ifield+1, value) == SQLITE_OK;
}

bool SQLite3Transceiver::SQLite3Statement::
bind_uint8(unsigned ifield, uint8_t value)
{
  if(make_bound_sql_)SQLTransceiver::Statement::bind_uint8(ifield, value);
  return sqlite3_bind_int(stmt_, ifield+1, value) == SQLITE_OK;
}

bool SQLite3Transceiver::SQLite3Statement::
bind_float(unsigned ifield, float value)
{
  if(make_bound_sql_)SQLTransceiver::Statement::bind_float(ifield, value);
  return sqlite3_bind_double(stmt_, ifield+1, value) == SQLITE_OK;
}

bool SQLite3Transceiver::SQLite3Statement::
bind_double(unsigned ifield, double value)
{
  if(make_bound_sql_)SQLTransceiver::Statement::bind_double(ifield, value);
  return sqlite3_bind_double(stmt_, ifield+1, value) == SQLITE_OK;
}

bool SQLite3Transceiver::SQLite3Statement::
bind_bool(unsigned ifield, bool value)
{
  if(make_bound_sql_)SQLTransceiver::Statement::bind_bool(ifield, value);
  return sqlite3_bind_int(stmt_, ifield+1, value?1:0) == SQLITE_OK;  
}

bool SQLite3Transceiver::SQLite3Statement::
bind_string(unsigned ifield, const std::string& value)
{
  if(make_bound_sql_)SQLTransceiver::Statement::bind_string(ifield, value);
  int nbyte = ::strlen(value.c_str())+1;
  char* string_data = static_cast<char*>(::malloc(nbyte));
  ::memcpy(string_data, value.c_str(), nbyte);
  return sqlite3_bind_text(stmt_, ifield+1, string_data, nbyte, ::free)
      == SQLITE_OK;
}

bool SQLite3Transceiver::SQLite3Statement::
bind_bytes(unsigned ifield, const std::string& value)
{
  if(make_bound_sql_)SQLTransceiver::Statement::bind_bytes(ifield, value);
  int nbyte = value.size();
  void* blob_data = ::malloc(nbyte);
  ::memcpy(blob_data, value.c_str(), nbyte);
  return sqlite3_bind_blob(stmt_, ifield+1, blob_data, nbyte, ::free)
      == SQLITE_OK;
}
