/*

   calin/io/sqlite3_statement.hpp -- Stephen Fegan -- 2015-09-24

   Derived SQLStatement for SQLite3 database

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

#pragma once

#include <sqlite3.h>

#include "io/sql_statement.hpp"

namespace calin { namespace io { namespace sql_transceiver {

class SQLite3Statement: public SQLStatement
{
public:
  SQLite3Statement(const std::string& sql, sqlite3* db,
                   bool make_bound_sql = false);
  virtual ~SQLite3Statement();

  unsigned num_columns() override;

  bool is_initialized() override;
  int error_code() override;
  std::string error_message() override;

  void reset() override;
  StepStatus step() override;
  uint64_t get_oid() override;

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

  bool column_is_null(unsigned icol, bool* good = nullptr) override;
  int64_t extract_int64(unsigned icol, bool* good = nullptr) override;
  int32_t extract_int32(unsigned icol, bool* good = nullptr) override;
  int16_t extract_int16(unsigned icol, bool* good = nullptr) override;
  int8_t extract_int8(unsigned icol, bool* good = nullptr) override;
  uint64_t extract_uint64(unsigned icol, bool* good = nullptr) override;
  uint32_t extract_uint32(unsigned icol, bool* good = nullptr) override;
  uint16_t extract_uint16(unsigned icol, bool* good = nullptr) override;
  uint8_t extract_uint8(unsigned icol, bool* good = nullptr) override;
  float extract_float(unsigned icol, bool* good = nullptr) override;
  double extract_double(unsigned icol, bool* good = nullptr) override;
  bool extract_bool(unsigned icol, bool* good = nullptr) override;
  std::string extract_string(unsigned icol, bool* good = nullptr) override;
  std::string extract_bytes(unsigned icol, bool* good = nullptr) override;

protected:
  void set_good(unsigned icol, bool* good);

  sqlite3* db_ = nullptr;
  sqlite3_stmt* stmt_ = nullptr;
  unsigned ndatacol_ = 0;
  bool make_bound_sql_ = false;
};

} } } // namespace calin::io::sql_transceiver
