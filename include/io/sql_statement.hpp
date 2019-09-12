/* 

   calin/io/sql_statement.hpp -- Stephen Fegan -- 2015-09-24

   Base class for SQL statements

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

#include <google/protobuf/descriptor.pb.h>

namespace calin { namespace io { namespace sql_transceiver {

class SQLStatement
{
 public:
  enum StepStatus { ERROR, OK_NO_DATA, OK_HAS_DATA };

  SQLStatement(const std::string& sql);
  virtual ~SQLStatement();

  const std::string& sql() const { return sql_; }
  std::string bound_sql() const;
    
  virtual unsigned num_columns();

  virtual bool is_initialized();
  virtual int error_code();
  virtual std::string error_message();
    
  virtual void reset();
  virtual StepStatus step();
  virtual uint64_t get_oid();
    
  bool bind_field(unsigned ifield, const google::protobuf::Message* m,
                  const google::protobuf::FieldDescriptor* d);
  bool bind_repeated_field(unsigned ifield, uint64_t loop_id, 
                           const google::protobuf::Message* m,
                           const google::protobuf::FieldDescriptor* d);

  virtual bool bind_null(unsigned ifield);
  virtual bool bind_int64(unsigned ifield, int64_t value);
  virtual bool bind_int32(unsigned ifield, int32_t value);
  virtual bool bind_int16(unsigned ifield, int16_t value);
  virtual bool bind_int8(unsigned ifield, int8_t value);
  virtual bool bind_uint64(unsigned ifield, uint64_t value);
  virtual bool bind_uint32(unsigned ifield, uint32_t value);
  virtual bool bind_uint16(unsigned ifield, uint16_t value);
  virtual bool bind_uint8(unsigned ifield, uint8_t value);
  virtual bool bind_float(unsigned ifield, float value);
  virtual bool bind_double(unsigned ifield, double value);
  virtual bool bind_bool(unsigned ifield, bool value);
  virtual bool bind_string(unsigned ifield, const std::string& value);
  virtual bool bind_bytes(unsigned ifield, const std::string& value);

  bool extract_field(unsigned icol, google::protobuf::Message* m,
                     const google::protobuf::FieldDescriptor* d);
  bool extract_repeated_field(unsigned icol, uint64_t loop_id, 
                              google::protobuf::Message* m,
                              const google::protobuf::FieldDescriptor* d);

  virtual bool column_is_null(unsigned icol, bool* good = nullptr);
  virtual int64_t extract_int64(unsigned icol, bool* good = nullptr);
  virtual int32_t extract_int32(unsigned icol, bool* good = nullptr);
  virtual int16_t extract_int16(unsigned icol, bool* good = nullptr);
  virtual int8_t extract_int8(unsigned icol, bool* good = nullptr);
  virtual uint64_t extract_uint64(unsigned icol, bool* good = nullptr);
  virtual uint32_t extract_uint32(unsigned icol, bool* good = nullptr);
  virtual uint16_t extract_uint16(unsigned icol, bool* good = nullptr);
  virtual uint8_t extract_uint8(unsigned icol, bool* good = nullptr);
  virtual float extract_float(unsigned icol, bool* good = nullptr);
  virtual double extract_double(unsigned icol, bool* good = nullptr);
  virtual bool extract_bool(unsigned icol, bool* good = nullptr);
  virtual std::string extract_string(unsigned icol, bool* good = nullptr);
  virtual std::string extract_bytes(unsigned icol, bool* good = nullptr);

 protected:
  std::string sql_;
  std::vector<std::string> bound_values_;
};

} } } // namespace calin::io::sql_transceiver
