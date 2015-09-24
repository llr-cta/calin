/* 

   calin/io/sql_statement.hpp -- Stephen Fegan -- 2015-09-24

   Base class for SQL statements

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

 protected:
  std::string sql_;
  std::vector<std::string> bound_values_;
};

} } } // namespace calin::io::sql_transceiver
