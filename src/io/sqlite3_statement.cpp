/* 

   calin/io/sqlite3_statement.cpp -- Stephen Fegan -- 2015-09-24

   Derived SQLStatement for SQLite3 database

*/

#include "io/log.hpp"
#include "io/sqlite3_statement.hpp"

using namespace calin::io::log;
using namespace calin::io::sql_transceiver;

SQLite3Statement::
SQLite3Statement(const std::string& sql, sqlite3* db, bool make_bound_sql):
    SQLStatement(sql), db_(db), stmt_(nullptr), make_bound_sql_(make_bound_sql)
{
  sqlite3_prepare_v2(db_, sql.c_str(), -1, &stmt_, nullptr);
}

SQLite3Statement::~SQLite3Statement()
{
  sqlite3_finalize(stmt_);
}

unsigned SQLite3Statement::num_columns()
{
  return sqlite3_column_count(stmt_);
}

bool SQLite3Statement::is_initialized()
{
  return stmt_ != nullptr;
}

int SQLite3Statement::error_code()
{
  return sqlite3_errcode(db_);
}

std::string SQLite3Statement::error_message()
{
  return sqlite3_errmsg(db_);
}

SQLStatement::StepStatus
SQLite3Statement::step()
{
  int retcode = sqlite3_step(stmt_);
  if(retcode == SQLITE_DONE)return SQLStatement::OK_NO_DATA;
  else if(retcode == SQLITE_ROW)return SQLStatement::OK_HAS_DATA;
  else return SQLStatement::ERROR;
}

uint64_t SQLite3Statement::get_oid()
{
  return sqlite3_last_insert_rowid(db_);
}

void SQLite3Statement::reset()
{
  sqlite3_reset(stmt_);
  sqlite3_clear_bindings(stmt_);
  if(make_bound_sql_)SQLStatement::reset();
}

bool SQLite3Statement::bind_null(unsigned ifield)
{
  if(make_bound_sql_)SQLStatement::bind_null(ifield);
  return sqlite3_bind_null(stmt_, ifield+1) == SQLITE_OK;
}

bool SQLite3Statement::
bind_int64(unsigned ifield, int64_t value)
{
  if(make_bound_sql_)SQLStatement::bind_int64(ifield, value);
  return sqlite3_bind_int64(stmt_, ifield+1, value) == SQLITE_OK;
}

bool SQLite3Statement::
bind_int32(unsigned ifield, int32_t value)
{
  if(make_bound_sql_)SQLStatement::bind_int32(ifield, value);
  return sqlite3_bind_int(stmt_, ifield+1, value) == SQLITE_OK;
}

bool SQLite3Statement::
bind_int16(unsigned ifield, int16_t value)
{
  if(make_bound_sql_)SQLStatement::bind_int16(ifield, value);
  return sqlite3_bind_int(stmt_, ifield+1, value) == SQLITE_OK;
}

bool SQLite3Statement::
bind_int8(unsigned ifield, int8_t value)
{
  if(make_bound_sql_)SQLStatement::bind_int8(ifield, value);
  return sqlite3_bind_int(stmt_, ifield+1, value) == SQLITE_OK;
}

bool SQLite3Statement::
bind_uint64(unsigned ifield, uint64_t value)
{
  if(make_bound_sql_)SQLStatement::bind_uint64(ifield, value);
  return sqlite3_bind_int64(stmt_, ifield+1, value) == SQLITE_OK;
}

bool SQLite3Statement::
bind_uint32(unsigned ifield, uint32_t value)
{
  if(make_bound_sql_)SQLStatement::bind_uint32(ifield, value);
  return sqlite3_bind_int(stmt_, ifield+1, value) == SQLITE_OK;
}

bool SQLite3Statement::
bind_uint16(unsigned ifield, uint16_t value)
{
  if(make_bound_sql_)SQLStatement::bind_uint16(ifield, value);
  return sqlite3_bind_int(stmt_, ifield+1, value) == SQLITE_OK;
}

bool SQLite3Statement::
bind_uint8(unsigned ifield, uint8_t value)
{
  if(make_bound_sql_)SQLStatement::bind_uint8(ifield, value);
  return sqlite3_bind_int(stmt_, ifield+1, value) == SQLITE_OK;
}

bool SQLite3Statement::
bind_float(unsigned ifield, float value)
{
  if(make_bound_sql_)SQLStatement::bind_float(ifield, value);
  return sqlite3_bind_double(stmt_, ifield+1, value) == SQLITE_OK;
}

bool SQLite3Statement::
bind_double(unsigned ifield, double value)
{
  if(make_bound_sql_)SQLStatement::bind_double(ifield, value);
  return sqlite3_bind_double(stmt_, ifield+1, value) == SQLITE_OK;
}

bool SQLite3Statement::
bind_bool(unsigned ifield, bool value)
{
  if(make_bound_sql_)SQLStatement::bind_bool(ifield, value);
  return sqlite3_bind_int(stmt_, ifield+1, value?1:0) == SQLITE_OK;  
}

bool SQLite3Statement::
bind_string(unsigned ifield, const std::string& value)
{
  if(make_bound_sql_)SQLStatement::bind_string(ifield, value);
  int nbyte = ::strlen(value.c_str())+1;
  char* string_data = static_cast<char*>(::malloc(nbyte));
  ::memcpy(string_data, value.c_str(), nbyte);
  return sqlite3_bind_text(stmt_, ifield+1, string_data, nbyte-1, ::free)
      == SQLITE_OK;
}

bool SQLite3Statement::
bind_bytes(unsigned ifield, const std::string& value)
{
  if(make_bound_sql_)SQLStatement::bind_bytes(ifield, value);
  int nbyte = value.size();
  void* blob_data = ::malloc(nbyte);
  ::memcpy(blob_data, value.c_str(), nbyte);
  return sqlite3_bind_blob(stmt_, ifield+1, blob_data, nbyte, ::free)
      == SQLITE_OK;
}

bool SQLite3Statement::column_is_null(unsigned icol, bool* good)
{
  if(good)*good = true;
  return sqlite3_column_type(stmt_, icol) == SQLITE_NULL;
}

int64_t SQLite3Statement::extract_int64(unsigned icol, bool* good)
{
  if(good)*good = true;
  return sqlite3_column_int64(stmt_, icol);
}

int32_t SQLite3Statement::extract_int32(unsigned icol, bool* good)
{
  if(good)*good = true;
  return sqlite3_column_int(stmt_, icol);  
}

int16_t SQLite3Statement::extract_int16(unsigned icol, bool* good)
{
  if(good)*good = true;
  return sqlite3_column_int(stmt_, icol);  
}

int8_t SQLite3Statement::extract_int8(unsigned icol, bool* good)
{
  if(good)*good = true;
  return sqlite3_column_int(stmt_, icol);  
}

uint64_t SQLite3Statement::extract_uint64(unsigned icol, bool* good)
{
  if(good)*good = true; return 0;
  return sqlite3_column_int64(stmt_, icol);
}

uint32_t SQLite3Statement::extract_uint32(unsigned icol, bool* good)
{
  if(good)*good = true;
  return sqlite3_column_int(stmt_, icol);  
}

uint16_t SQLite3Statement::extract_uint16(unsigned icol, bool* good)
{
  if(good)*good = true;
  return sqlite3_column_int(stmt_, icol);  
}

uint8_t SQLite3Statement::extract_uint8(unsigned icol, bool* good)
{
  if(good)*good = true;
  return sqlite3_column_int(stmt_, icol);  
}

float SQLite3Statement::extract_float(unsigned icol, bool* good)
{
  if(good)*good = true;
  return sqlite3_column_double(stmt_, icol);  
}

double SQLite3Statement::extract_double(unsigned icol, bool* good)
{
  if(good)*good = true;
  return sqlite3_column_double(stmt_, icol);  
}

bool SQLite3Statement::extract_bool(unsigned icol, bool* good)
{
  if(good)*good = true;
  return sqlite3_column_int(stmt_, icol) != 0;
}

std::string SQLite3Statement::extract_string(unsigned icol, bool* good)
{
  if(good)*good = true;
  const char* text = (const char*) sqlite3_column_text(stmt_, icol);
  return std::string(text, sqlite3_column_bytes(stmt_, icol));
}

std::string SQLite3Statement::extract_bytes(unsigned icol, bool* good)
{
  if(good)*good = true;
  const char* blob = (const char*) sqlite3_column_blob(stmt_, icol);
  if(blob == nullptr)return std::string();
  return std::string(blob, sqlite3_column_bytes(stmt_, icol));
}
