/* 

   calin/io/sql_statement.hpp -- Stephen Fegan -- 2015-09-24

   Base class for SQL statements

*/

#include <sstream>

#include "proto/calin.pb.h"
#include "io/log.hpp"
#include "io/sql_statement.hpp"

using namespace calin::io::log;
using namespace calin::io::sql_transceiver;
using namespace google::protobuf;

SQLStatement::SQLStatement(const std::string& sql):
    sql_(sql), bound_values_()
{
  // nothing to see here
}

SQLStatement::~SQLStatement()
{
  // nothing to see here
}

std::string SQLStatement::bound_sql() const
{
  std::ostringstream L;
  unsigned ifield = 0;
  for(auto c : sql_)
  {
    if(c=='?' and ifield<bound_values_.size())
      L << bound_values_[ifield++];
    else L << c;
  }
  return L.str();
}

unsigned SQLStatement::num_columns()
{
  return 0;
}

bool SQLStatement::is_initialized()
{
  return true;
}

int SQLStatement::error_code()
{
  return 0;
}

std::string SQLStatement::error_message()
{
  return "";
}

void SQLStatement::reset()
{
  bound_values_.clear();
}

SQLStatement::StepStatus
SQLStatement::step()
{
  return SQLStatement::OK_NO_DATA;
}

uint64_t SQLStatement::get_oid()
{
  return 0;
}

bool SQLStatement::
bind_field(unsigned ifield, const google::protobuf::Message* m,
           const google::protobuf::FieldDescriptor* d)
{
  calin::FieldOptions::Int32StorageType int32_type =
      calin::FieldOptions::INT_32;
  const google::protobuf::FieldOptions* fopt { &d->options() };
  if(fopt->HasExtension(CFO))
    int32_type = fopt->GetExtension(CFO).int32_type();

  switch(d->type())
  {
    case FieldDescriptor::TYPE_DOUBLE:
      return bind_double(ifield, m->GetReflection()->GetDouble(*m, d));
    case FieldDescriptor::TYPE_FLOAT:
      return bind_float(ifield, m->GetReflection()->GetFloat(*m, d));
    case FieldDescriptor::TYPE_SFIXED64: // fallthrough
    case FieldDescriptor::TYPE_SINT64:   // fallthrough
    case FieldDescriptor::TYPE_INT64:
      return bind_int64(ifield, m->GetReflection()->GetInt64(*m, d));
    case FieldDescriptor::TYPE_FIXED64:  // fallthrough
    case FieldDescriptor::TYPE_UINT64:
      return bind_uint64(ifield, m->GetReflection()->GetUInt64(*m, d));
    case FieldDescriptor::TYPE_SFIXED32: // fallthrough
    case FieldDescriptor::TYPE_SINT32:   // fallthrough
    case FieldDescriptor::TYPE_INT32:
      switch(int32_type) {
	case FieldOptions::INT_16:
          return bind_int16(ifield, m->GetReflection()->GetInt32(*m, d));
	case FieldOptions::INT_8:
          return bind_int8(ifield, m->GetReflection()->GetInt32(*m, d));
	case FieldOptions::INT_32:       // fallthrough
        default:
          return bind_int32(ifield, m->GetReflection()->GetInt32(*m, d));
	};
    case FieldDescriptor::TYPE_FIXED32:
    case FieldDescriptor::TYPE_UINT32:
      switch(int32_type) {
	case FieldOptions::INT_16:
          return bind_uint16(ifield, m->GetReflection()->GetUInt32(*m, d));
	case FieldOptions::INT_8:
          return bind_uint8(ifield, m->GetReflection()->GetUInt32(*m, d));
	case FieldOptions::INT_32:       // fall through
        default:
          return bind_uint32(ifield, m->GetReflection()->GetUInt32(*m, d));
      }
    case FieldDescriptor::TYPE_BOOL:
      return bind_bool(ifield, m->GetReflection()->GetBool(*m, d));
    case FieldDescriptor::TYPE_STRING:
      return bind_string(ifield, m->GetReflection()->GetString(*m, d));
    case FieldDescriptor::TYPE_BYTES:
      return bind_bytes(ifield, m->GetReflection()->GetString(*m, d));
    case FieldDescriptor::TYPE_ENUM:
      return bind_int32(ifield, m->GetReflection()->GetEnumValue(*m, d));
    case FieldDescriptor::TYPE_MESSAGE:  // fallthrough to assert(0)
    case FieldDescriptor::TYPE_GROUP:    // fallthrough to assert(0)
    default:
      break;
  }
  assert(0);
}

bool SQLStatement::
bind_repeated_field(unsigned ifield, uint64_t iloop,
                    const google::protobuf::Message* m,
                    const google::protobuf::FieldDescriptor* d)
{
  calin::FieldOptions::Int32StorageType int32_type =
      calin::FieldOptions::INT_32;
  const google::protobuf::FieldOptions* fopt { &d->options() };
  if(fopt->HasExtension(CFO))
    int32_type = fopt->GetExtension(CFO).int32_type();

  switch(d->type())
  {
    case FieldDescriptor::TYPE_DOUBLE:
      return bind_double(ifield, m->GetReflection()->
                         GetRepeatedDouble(*m, d, iloop));
    case FieldDescriptor::TYPE_FLOAT:
      return bind_float(ifield, m->GetReflection()->
                        GetRepeatedFloat(*m, d, iloop));
    case FieldDescriptor::TYPE_SFIXED64: // fallthrough
    case FieldDescriptor::TYPE_SINT64:   // fallthrough
    case FieldDescriptor::TYPE_INT64:
      return bind_int64(ifield, m->GetReflection()->
                        GetRepeatedInt64(*m, d, iloop));
    case FieldDescriptor::TYPE_FIXED64:  // fallthrough
    case FieldDescriptor::TYPE_UINT64:
      return bind_uint64(ifield,
                         m->GetReflection()->GetRepeatedUInt64(*m, d, iloop));
    case FieldDescriptor::TYPE_SFIXED32: // fallthrough
    case FieldDescriptor::TYPE_SINT32:   // fallthrough
    case FieldDescriptor::TYPE_INT32:
      switch(int32_type) {
	case FieldOptions::INT_16:
          return bind_int16(ifield, m->GetReflection()->
                            GetRepeatedInt32(*m,d,iloop));
	case FieldOptions::INT_8:
          return bind_int8(ifield, m->GetReflection()->
                           GetRepeatedInt32(*m,d,iloop));
	case FieldOptions::INT_32:       // fallthrough
        default:
          return bind_int32(ifield, m->GetReflection()->
                            GetRepeatedInt32(*m,d,iloop));
	};
    case FieldDescriptor::TYPE_FIXED32:
    case FieldDescriptor::TYPE_UINT32:
      switch(int32_type) {
	case FieldOptions::INT_16:
          return bind_uint16(ifield, m->GetReflection()->
                             GetRepeatedUInt32(*m, d, iloop));
	case FieldOptions::INT_8:
          return bind_uint8(ifield, m->GetReflection()->
                            GetRepeatedUInt32(*m, d, iloop));
	case FieldOptions::INT_32:       // fall through
        default:
          return bind_uint32(ifield, m->GetReflection()->
                             GetRepeatedUInt32(*m, d, iloop));
      }
    case FieldDescriptor::TYPE_BOOL:
      return bind_bool(ifield, m->GetReflection()->
                       GetRepeatedBool(*m, d, iloop));
    case FieldDescriptor::TYPE_STRING:
      return bind_string(ifield, m->GetReflection()->
                         GetRepeatedString(*m, d, iloop));
    case FieldDescriptor::TYPE_BYTES:
      return bind_bytes(ifield, m->GetReflection()->
                        GetRepeatedString(*m, d, iloop));
    case FieldDescriptor::TYPE_ENUM:
      return bind_int32(ifield, m->GetReflection()->
                        GetRepeatedEnumValue(*m,d,iloop));
    case FieldDescriptor::TYPE_MESSAGE:  // fallthrough to assert(0)
    case FieldDescriptor::TYPE_GROUP:    // fallthrough to assert(0)
    default:
      break;
  }
  assert(0);
}

bool SQLStatement::bind_null(unsigned ifield)
{
  return SQLStatement::bind_string(ifield, "NULL");
}

bool SQLStatement::bind_int64(unsigned ifield, int64_t value)
{
  return SQLStatement::bind_string(ifield, std::to_string(value));
}

bool SQLStatement::bind_int32(unsigned ifield, int32_t value)
{
  return SQLStatement::bind_string(ifield, std::to_string(value));
}

bool SQLStatement::bind_int16(unsigned ifield, int16_t value)
{
  return SQLStatement::bind_string(ifield, std::to_string(value));
}

bool SQLStatement::bind_int8(unsigned ifield, int8_t value)
{
  return SQLStatement::bind_string(ifield, std::to_string(value));
}

bool SQLStatement::bind_uint64(unsigned ifield, uint64_t value)
{
  return SQLStatement::bind_string(ifield, std::to_string(value));
}

bool SQLStatement::bind_uint32(unsigned ifield, uint32_t value)
{
  return SQLStatement::bind_string(ifield, std::to_string(value));
}

bool SQLStatement::bind_uint16(unsigned ifield, uint16_t value)
{
  return SQLStatement::bind_string(ifield, std::to_string(value));
}

bool SQLStatement::bind_uint8(unsigned ifield, uint8_t value)
{
  return SQLStatement::bind_string(ifield, std::to_string(value));
}

bool SQLStatement::bind_float(unsigned ifield, float value)
{
  return SQLStatement::bind_string(ifield, std::to_string(value));
}

bool SQLStatement::bind_double(unsigned ifield, double value)
{
  return SQLStatement::bind_string(ifield, std::to_string(value));
}

bool SQLStatement::bind_bool(unsigned ifield, bool value)
{
  return SQLStatement::bind_string(ifield, std::to_string(value));
}

bool SQLStatement::
bind_string(unsigned ifield, const std::string& value)
{
  if(ifield >= bound_values_.size())bound_values_.resize(ifield+1);
  bound_values_[ifield] = value;
  return true;
}

bool SQLStatement::
bind_bytes(unsigned ifield, const std::string& value)
{
  return SQLStatement::bind_string(ifield, value);
}
