/*

   calin/io/sqlite3_serializer.hpp -- Stephen Fegan -- 2020-04-03

   Provides reading and writing of protobuf structures to SQLite3 databases

   Based on : calin/io/sqlite3_transceiver.hpp -- Stephen Fegan -- 2015-09-08

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

#pragma once

#include <string>
#include <streambuf>
#include <iostream>
#include <sstream>
#include <vector>

#include <google/protobuf/descriptor.pb.h>

#include <sqlite3.h>

#include "io/sql_serializer.hpp"
#include "io/sqlite3_statement.hpp"

namespace calin { namespace io { namespace sql_serializer {

class SQLite3Serializer: public SQLSerializer
{
public:
  enum OpenMode { EXISTING_OR_NEW_RW, EXISTING_RW, TRUNCATE_RW, READ_ONLY, READ_ONLY_NON_CALIN_DB };

  SQLite3Serializer(const std::string& filename,
                    OpenMode open_mode = EXISTING_OR_NEW_RW,
                    bool write_sql_to_log = false);
  virtual ~SQLite3Serializer();

protected:
  sqlite3* db_ = nullptr;
  bool adopt_db_ = false;
  OpenMode open_mode_ = EXISTING_OR_NEW_RW;

  SQLStatement* prepare_statement(const std::string& sql) override;

  bool begin_transaction() override;
  bool commit_transaction() override;
  bool rollback_transaction() override;

  std::string sql_select_field_spec(const SQLTableField* f) override;
  std::string sql_insert_field_spec(const SQLTableField* f) override;
  std::string sql_type(const google::protobuf::FieldDescriptor* d) override;
  std::string sql_add_field_to_table(const SQLTableField* f) override;
};

} } } // namespace calin::io::sql_transceiver
