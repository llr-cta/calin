/*

   calin/io/sqlite3_transceiver.hpp -- Stephen Fegan -- 2015-09-08

   Provides reading and writing of protobuf structures to SQLite3 databases

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

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
  sqlite3* db_ = nullptr;
  bool adopt_db_ = false;
  OpenMode open_mode_ = EXISTING_OR_NEW_RW;

  SQLStatement* prepare_statement(const std::string& sql) override;

  bool begin_transaction() override;
  bool commit_transaction() override;
  bool rollback_transaction() override;
};

} } } // namespace calin::io::sql_transceiver
