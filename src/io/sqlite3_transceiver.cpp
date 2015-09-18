/* 

   calin/io/sqlite3_transceiver.cpp -- Stephen Fegan -- 2015-09-08

   Provides reading and writing of protobuf structures to SQLite3 databases

*/

#include <unistd.h>

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

bool SQLite3Transceiver::
insert(const std::string& table_name,
       const google::protobuf::Message* m_data,
       const google::protobuf::Message* m_key, bool write_sql_to_log)
{
  
}

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

