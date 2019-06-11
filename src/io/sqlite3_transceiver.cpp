/*

   calin/io/sqlite3_transceiver.cpp -- Stephen Fegan -- 2015-09-08

   Provides reading and writing of protobuf structures to SQLite3 databases

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <unistd.h>
#include <cstring>
#include <stdexcept>

#include <provenance/chronicle.hpp>
#include <util/log.hpp>
#include <io/sqlite3_transceiver.hpp>
#include <util/file.hpp>

using namespace calin::util::log;
using namespace calin::io::sql_transceiver;

SQLite3Transceiver::
SQLite3Transceiver(const std::string& filename_in, OpenMode open_mode,
                   bool write_sql_to_log):
    SQLTransceiver(write_sql_to_log), adopt_db_(true), open_mode_(open_mode)
{
  std::string filename = calin::util::file::expand_filename(filename_in);
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
  if(sqlite3_open_v2(filename.c_str(), &db_, flags, NULL) != SQLITE_OK)
    throw std::runtime_error("Could not open: " + filename_in + "\n" +
      sqlite3_errmsg(db_));

  calin::ix::provenance::chronicle::AccessType access;
  switch(open_mode) {
  case EXISTING_OR_NEW_RW:
  case EXISTING_RW:
    access = calin::ix::provenance::chronicle::AT_READWRITE; break;
  case TRUNCATE_RW:
    access = calin::ix::provenance::chronicle::AT_TRUNCATE; break;
  case READ_ONLY:
  default:
    access = calin::ix::provenance::chronicle::AT_READ; break;
  }

  calin::provenance::chronicle::register_file_open(filename, access,
    __PRETTY_FUNCTION__);
}

SQLite3Transceiver::~SQLite3Transceiver()
{
  if(adopt_db_)sqlite3_close(db_);
}

#if 0
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
#endif

SQLStatement* SQLite3Transceiver::prepare_statement(const std::string& sql)
{
  return new SQLite3Statement(sql, db_, true);
}

bool SQLite3Transceiver::begin_transaction()
{
  return sqlite3_exec(db_, "BEGIN TRANSACTION", NULL, NULL, NULL) == SQLITE_OK;
}

bool SQLite3Transceiver::commit_transaction()
{
  return sqlite3_exec(db_, "COMMIT TRANSACTION", NULL, NULL, NULL) == SQLITE_OK;
}

bool SQLite3Transceiver::rollback_transaction()
{
  return
      sqlite3_exec(db_, "ROLLBACK TRANSACTION", NULL, NULL, NULL) == SQLITE_OK;
}

std::string SQLite3Transceiver::sql_select_field_spec(const SQLTableField* f)
{
  if(f->field_d) {
    const google::protobuf::FieldOptions* fopt { &f->field_d->options() };
    if(fopt->HasExtension(CFO)) {
      const auto& cfo = fopt->GetExtension(CFO);
      switch(cfo.sql().transform()) {
        case SQLFieldOptions::TRANSFORM_UNIXTIME_TOFROM_DATETIME:
          return "STRFTIME('%s'," + sql_field_name(f->field_name) + ")";
        case SQLFieldOptions::TRANSFORM_NONE:
        default:
          return sql_field_name(f->field_name);
      }
    }
  }
  return sql_field_name(f->field_name);
}

std::string SQLite3Transceiver::sql_insert_field_spec(const SQLTableField* f)
{
  if(f->field_d) {
    const google::protobuf::FieldOptions* fopt { &f->field_d->options() };
    if(fopt->HasExtension(CFO)) {
      const auto& cfo = fopt->GetExtension(CFO);
      switch(cfo.sql().transform()) {
        case SQLFieldOptions::TRANSFORM_UNIXTIME_TOFROM_DATETIME:
          return "DATETIME(?,'UNIXEPOCH')";
        case SQLFieldOptions::TRANSFORM_NONE:
        default:
          return "?";
      }
    }
  }
  return "?";
}

std::string SQLite3Transceiver::sql_type(const google::protobuf::FieldDescriptor* d)
{
  if(d) {
    const google::protobuf::FieldOptions* fopt { &d->options() };
    if(fopt->HasExtension(CFO)) {
      const auto& cfo = fopt->GetExtension(CFO);
      switch(cfo.sql().transform()) {
        case SQLFieldOptions::TRANSFORM_UNIXTIME_TOFROM_DATETIME:
          return "DATETIME";
        case SQLFieldOptions::TRANSFORM_NONE:
        default:
          return SQLTransceiver::sql_type(d);
      }
    }
  }
  return SQLTransceiver::sql_type(d);
}
