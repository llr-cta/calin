/*

   calin/io/sqlite3_serializer.cpp -- Stephen Fegan -- 2020-04-03

   Provides reading and writing of protobuf structures to SQLite3 databases

   Based on : calin/io/sqlite3_transceiver.cpp -- Stephen Fegan -- 2015-09-08

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

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
#include <io/sqlite3_serializer.hpp>
#include <util/file.hpp>

using namespace calin::util::log;
using namespace calin::io::sql_serializer;
using namespace calin::io::sql_transceiver;

SQLite3Serializer::
SQLite3Serializer(const std::string& filename_in, OpenMode open_mode,
                  bool write_sql_to_log):
    SQLSerializer(write_sql_to_log), adopt_db_(true), open_mode_(open_mode)
{
  std::string filename = calin::util::file::expand_filename(filename_in);
  if(open_mode == TRUNCATE_RW and filename.front() != ':' and
     filename.back() != ':')unlink(filename.c_str());

  bool exists = false;
  if(filename.front() != ':' and filename.back() != ':')
    exists = calin::util::file::is_file(filename);

  int flags { 0 };
  switch(open_mode) {
  case EXISTING_OR_NEW_RW:
  case TRUNCATE_RW:
    flags = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE; break;
  case EXISTING_RW:
    flags = SQLITE_OPEN_READWRITE; break;
  case READ_ONLY:
  case READ_ONLY_NON_CALIN_DB:
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
  case READ_ONLY_NON_CALIN_DB:
  default:
    access = calin::ix::provenance::chronicle::AT_READ; break;
  }

  file_record_ = calin::provenance::chronicle::register_file_open(filename, access,
    __PRETTY_FUNCTION__);

  if(exists) {
    if(open_mode != READ_ONLY_NON_CALIN_DB) {
      retrieve_db_tables_and_fields();
    }
  } else {
    create_db_tables_and_fields();
  }
}

SQLite3Serializer::~SQLite3Serializer()
{
  // Must delete all open prepared statements here before calling sqlite3_close
  // so just delete all the SQLTables in the schema_, which in turn delete their
  // open select & insert statements, then clear the schema_ so the base class
  // doesn't delete them again
  for(auto& ischema : schema_) {
    for(auto& it : ischema.second) {
      delete it.second;
    }
  }
  schema_.clear();
  if(adopt_db_)
  {
    if(sqlite3_close(db_) !=  SQLITE_OK) {
      LOG(ERROR) << "SQLite3Serializer::~SQLite3Serializer: database connection did not close cleanly";
    }
    if(file_record_) {
      calin::provenance::chronicle::register_file_close(file_record_);
    }
  }
}

SQLStatement* SQLite3Serializer::prepare_statement(const std::string& sql)
{
  return new SQLite3Statement(sql, db_, /* make_bound_sql= */ write_sql_to_log_);
}

bool SQLite3Serializer::begin_transaction()
{
  if(transaction_count_.fetch_add(1) == 0) {
    if(write_sql_to_log_) {
      LOG(INFO) << "BEGIN TRANSACTION";
    }
    return sqlite3_exec(db_, "BEGIN TRANSACTION", NULL, NULL, NULL) == SQLITE_OK;
  }
  return true;
}

bool SQLite3Serializer::commit_transaction()
{
  if(transaction_count_.fetch_sub(1) == 1) {
    if(write_sql_to_log_) {
      LOG(INFO) << "COMMIT TRANSACTION";
    }
    return sqlite3_exec(db_, "COMMIT TRANSACTION", NULL, NULL, NULL) == SQLITE_OK;
  }
  return true;
}

bool SQLite3Serializer::rollback_transaction()
{
  if(transaction_count_.fetch_sub(1) == 1) {
    if(write_sql_to_log_) {
      LOG(INFO) << "ROLLBACK TRANSACTION";
    }
    return sqlite3_exec(db_, "ROLLBACK TRANSACTION", NULL, NULL, NULL) == SQLITE_OK;
  }
  return true;
}

std::string SQLite3Serializer::sql_select_field_spec(const SQLTableField* f)
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

std::string SQLite3Serializer::sql_insert_field_spec(const SQLTableField* f)
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

std::string SQLite3Serializer::sql_type(const google::protobuf::FieldDescriptor* d)
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
          return SQLSerializer::sql_type(d);
      }
    }
  }
  return SQLSerializer::sql_type(d);
}

std::string SQLite3Serializer::sql_add_field_to_table(const SQLTableField* f)
{
  // SQLite3 does not support a "comment" in the SQL for ALTER TABLE
  std::ostringstream sql;
  sql << "ALTER TABLE " << sql_table_name(f->table->table_name)
    << " ADD COLUMN " << sql_field_name(f->field_name) << ' ' << sql_type(f->field_d);
  return sql.str();
}
