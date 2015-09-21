/* 

   calin/io/sql_transceiver.cpp -- Stephen Fegan -- 2015-09-08

   Base class for reading and writing protobuf structures to SQL databases

*/

#include <deque>
#include <stdexcept>

#include "proto/calin.pb.h"
#include "io/sql_transceiver.hpp"
#include "io/log.hpp"

using namespace calin::io::log;
using namespace calin::io::sql_transceiver;
using namespace google::protobuf;

SQLTransceiver::~SQLTransceiver()
{
  // nothing to see here
}

bool SQLTransceiver::
create_tables(const std::string& table_name,
              const google::protobuf::Descriptor* d_data,
              const google::protobuf::Descriptor* d_key,
              const std::string& instance_desc,
              bool write_sql_to_log)
{
  SQLTable* t = make_keyed_sqltable_tree(table_name, d_data, d_key, false);
  if(write_sql_to_log)
    iterate_over_tables(t, [this](const SQLTable* it){
        LOG(VERBOSE) << sql_create_table(it) << '\n'; });
  delete t;
  return true;
}

bool SQLTransceiver::
insert(const std::string& table_name,
       const google::protobuf::Message* m_data,
       const google::protobuf::Message* m_key,
       bool write_sql_to_log)
{
  SQLTable* t =
      make_keyed_sqltable_tree(table_name, m_data->GetDescriptor(),
                               m_key?m_key->GetDescriptor():nullptr, false);
  prepare_insert_statements(t);
  insert_tables(t, m_data, m_key, 0, 0, write_sql_to_log);
  finalize_insert_statements(t);
  delete t;
  return true;
}

std::vector<std::pair<std::string,std::string> >
SQLTransceiver::list_all_table_columns(SQLTable* t)
{
  std::vector<std::pair<std::string,std::string>> fields;
  for(auto f : t->fields)
    fields.push_back(std::make_pair(t->table_name, f->field_name));
  for(auto st : t->sub_tables)
  {
    std::vector<std::pair<std::string,std::string>> sub_fields =
        list_all_table_columns(st);
    fields.insert(fields.end(), sub_fields.begin(), sub_fields.end());
  }
  return fields;
}
  
bool SQLTransceiver::is_single_table(const google::protobuf::Descriptor* d)
{
  for(int ifield = 0; ifield<d->field_count(); ifield++)
  {
    const FieldDescriptor* f { d->field(ifield) };
    const google::protobuf::FieldOptions* fopt { &f->options() };
    if(f->is_repeated())return false;
    if((f->type()==FieldDescriptor::TYPE_MESSAGE) and
       (!fopt->HasExtension(CFO) or
        !fopt->GetExtension(CFO).sql().inline_message() or
        !is_single_table(f->message_type())))return false;
  }
  return true;
}

calin::io::sql_transceiver::SQLTransceiver::SQLTable* SQLTransceiver::
make_keyed_sqltable_tree(const std::string& table_name,
                         const google::protobuf::Descriptor* d_data,
                         const google::protobuf::Descriptor* d_key,
                         bool ignore_key_option)
{
  SQLTable* t { make_sqltable_tree(table_name, d_data) };
  if(d_key)
  {
    SQLTable* key_t { make_extkey_tree("key", d_key) };
    t->fields.insert(t->fields.begin(),
                     key_t->fields.begin(), key_t->fields.end());
    key_t->fields.clear();
    delete key_t;
  }
  propagate_keys(t);
  return t;
}

calin::io::sql_transceiver::SQLTransceiver::SQLTable* SQLTransceiver::
make_sqltable_tree(const std::string& table_name,
                   const google::protobuf::Descriptor* d,
                   bool ignore_key_option)
{
  return r_make_sqltable_tree(table_name, d, nullptr, nullptr,
                              ignore_key_option);
}

calin::io::sql_transceiver::SQLTransceiver::SQLTable* SQLTransceiver::
make_extkey_tree(const std::string& table_name,
                 const google::protobuf::Descriptor* d)
{
  SQLTable* t =
      r_make_sqltable_tree(table_name, d, nullptr, nullptr, true);
  if(!t->sub_tables.empty())
    throw std::invalid_argument(std::string("SQL key message ") +
                                d->full_name() +
                                std::string(" cannot have any sub tables"));
  for(auto f : t->fields)
  {
    f->field_name = user_key_name(f->field_name);
    f->field_type = SQLTableField::KEY_USER_SUPPLIED;
  }
  return t;
}

calin::io::sql_transceiver::SQLTransceiver::SQLTable* SQLTransceiver::
r_make_sqltable_tree(const std::string& table_name,
                     const google::protobuf::Descriptor* d,
                     SQLTable* parent_table,
                     const google::protobuf::FieldDescriptor* parent_field_d,
                     bool ignore_key_option)
{
  SQLTable* t { new SQLTable };
  t->parent_table        = parent_table;
  t->table_name          = table_name;
  t->parent_field_d      = parent_field_d;

  for(int ifield = 0; ifield<d->field_count(); ifield++)
  {
    const FieldDescriptor* f { d->field(ifield) };
    const google::protobuf::FieldOptions* fopt { &f->options() };

    if(fopt->HasExtension(CFO) and
       (fopt->GetExtension(CFO).dont_store() or
        fopt->GetExtension(CFO).sql().dont_store()))
      continue;

    SQLTable* sub_table { nullptr };
    
    if(f->type()==FieldDescriptor::TYPE_MESSAGE)
    {
      sub_table = r_make_sqltable_tree(sub_name(table_name, f->name()),
                                       f->message_type(), t, f, true);

      if((parent_field_d and parent_field_d->is_map()) or
         (!f->is_repeated() and fopt->HasExtension(CFO) and
          fopt->GetExtension(CFO).sql().inline_message()))
      {
        // Message can be inlined - so move its fields into primary table

        for(auto ifield : sub_table->fields)
        {
          ifield->table           = t;
          ifield->field_name      = sub_name(f->name(), ifield->field_name);
          ifield->table           = t;
          ifield->field_d_path.insert(ifield->field_d_path.begin(), f);
          t->fields.push_back(ifield);
        }
        sub_table->fields.clear();

        for(auto itable : sub_table->sub_tables)
        {
          itable->parent_table = t;
          itable->parent_field_d_path.
              insert(itable->parent_field_d_path.begin(), f);
          t->sub_tables.push_back(itable);
        }
        sub_table->sub_tables.clear();

        delete(sub_table);
        continue;
      }
    }
    else if(f->is_repeated())
    {
      sub_table = new SQLTable;
      sub_table->parent_table    = t;
      sub_table->table_name      = sub_name(table_name,f->name());
      sub_table->parent_field_d  = f;      
      
      SQLTableField* tf { new SQLTableField };
      tf->table             = sub_table;
      tf->field_origin      = tf;
      tf->field_type        = SQLTableField::POD;
      tf->field_name        = f->name();
      tf->field_d           = f;

      sub_table->fields.push_back(tf);
    }
    else
    {
      SQLTableField* tf { new SQLTableField };
      tf->table             = t;
      tf->field_origin      = tf;
      tf->field_type        = SQLTableField::POD;
      tf->field_name        = f->name();
      tf->field_d           = f;
      
      if(parent_table == nullptr and ignore_key_option==false and
         fopt->HasExtension(CFO) and fopt->GetExtension(CFO).sql().is_key())
        tf->field_type      = SQLTableField::KEY_PROTO_DEFINED;

      t->fields.push_back(tf);
      continue;
    }

    assert(sub_table);
    
    if(f->is_map())
    {
      sub_table->fields.front()->field_name =
          sub_name(sub_table->table_name,
                   sub_table->fields.front()->field_name);
      sub_table->fields.front()->field_type   = SQLTableField::KEY_MAP_KEY;
    }
    else if(f->is_repeated())
    {
      SQLTableField* tf { new SQLTableField };
      tf->table             = sub_table;
      tf->field_origin      = tf;
      tf->field_type        = SQLTableField::KEY_LOOP_ID;
      tf->field_name        = sub_name(sub_table->table_name, "loop_id");
      tf->field_d           = nullptr;
      sub_table->fields.insert(sub_table->fields.begin(), tf);
    }

    t->sub_tables.insert(t->sub_tables.end(), sub_table);
  }
  return t;
}

void SQLTransceiver::prune_empty_tables(SQLTable* t)
{
  iterate_over_tables(t, [](SQLTable* t) {
      if(t->fields.empty()) {
        for(auto st : t->sub_tables) {
          st->parent_table = t->parent_table;
          st->parent_field_d_path.insert(st->parent_field_d_path.begin(),
                                         t->parent_field_d); }
        auto sti = std::find(t->parent_table->sub_tables.begin(),
                             t->parent_table->sub_tables.end(), t);
        assert(sti != t->parent_table->sub_tables.end());
        t->parent_table->sub_tables.insert(sti, t->sub_tables.begin(),
                                           t->sub_tables.end());        
        t->sub_tables.clear();
        sti = std::find(t->parent_table->sub_tables.begin(),
                        t->parent_table->sub_tables.end(), t);
        t->parent_table->sub_tables.erase(sti);
        delete t;
      } });  
}

void SQLTransceiver::propagate_keys(SQLTable* t,
                                    std::vector<const SQLTableField*> keys)
{
  for(unsigned ikey=0;ikey<keys.size();ikey++)
  {
    const SQLTableField* key = keys[keys.size()-ikey-1];
    SQLTableField* f = new SQLTableField(*key);
    if(f->field_type == SQLTableField::KEY_PROTO_DEFINED)
      f->field_name = sub_name(key->table->table_name, f->field_name);
    f->field_type = SQLTableField::KEY_INHERITED;
    t->fields.insert(t->fields.begin(), f);
  }
  for(unsigned ikey=keys.size();ikey<t->fields.size();ikey++)
    if(t->fields[ikey]->is_key())keys.push_back(t->fields[ikey]);
  if(keys.empty())
  {
    // No keys, add OID
    SQLTableField* f { new SQLTableField };
    enum FieldType {
      KEY_INHERITED, KEY_USER_SUPPLIED, KEY_PROTO_DEFINED, KEY_OID,
      KEY_LOOP_ID, KEY_MAP_KEY, POD };
    f->table          = t;
    f->field_origin   = f;
    f->field_type     = SQLTableField::KEY_OID;
    f->field_name     = sub_name(t->table_name,"oid");
    f->field_d        = nullptr;
    keys.push_back(f);
  }
  for(auto it : t->sub_tables)
    propagate_keys(it, keys);
}

std::string SQLTransceiver::user_key_name(const std::string& name)
{
  std::string new_name { "@" };
  new_name += name;
  return new_name;
}

std::string SQLTransceiver::
sub_name(const std::string& parent_name,const std::string& name)
{
  std::string new_name { parent_name };
  new_name += ".";
  new_name += name;
  return new_name;
}

std::string SQLTransceiver::sql_table_name(const std::string& name)
{
  std::string new_name { '`' };
  new_name += name;
  new_name += '`';
  return new_name;
}

std::string SQLTransceiver::sql_field_name(const std::string& name)
{
  std::string new_name { '`' };
  new_name += name;
  new_name += '`';
  return new_name;
}

std::string SQLTransceiver::
sql_type(const google::protobuf::FieldDescriptor* d)
{
  // The base class implements MySQL types since SQLite3 is very
  // forgiving about what it accepts as types
  
  if(!d)return "BIGINT UNSIGNED"; // OID or vector index

  calin::FieldOptions::Int32StorageType int32_type =
      calin::FieldOptions::INT_32;
  const google::protobuf::FieldOptions* fopt { &d->options() };
  if(fopt->HasExtension(CFO))
    int32_type = fopt->GetExtension(CFO).int32_type();

  switch(d->type())
  {
    case FieldDescriptor::TYPE_DOUBLE:   return "DOUBLE";
    case FieldDescriptor::TYPE_FLOAT:    return "FLOAT";
    case FieldDescriptor::TYPE_SFIXED64: // fallthrough
    case FieldDescriptor::TYPE_SINT64:   // fallthrough
    case FieldDescriptor::TYPE_INT64:    return "BIGINT";
    case FieldDescriptor::TYPE_FIXED64:  // fallthrough
    case FieldDescriptor::TYPE_UINT64:   return "BIGINT UNSIGNED";
    case FieldDescriptor::TYPE_SFIXED32: // fallthrough
    case FieldDescriptor::TYPE_SINT32:   // fallthrough
    case FieldDescriptor::TYPE_INT32:
      switch(int32_type) {
	case FieldOptions::INT_16:       return "SMALLINT";
	case FieldOptions::INT_8:        return "TINYINT";
	case FieldOptions::INT_32:       // fallthrough
        default:                         return "INT";
	};
    case FieldDescriptor::TYPE_BOOL:     return "BOOLEAN";
    case FieldDescriptor::TYPE_STRING:   return "TEXT";
    case FieldDescriptor::TYPE_BYTES:    return "BLOB";
    case FieldDescriptor::TYPE_FIXED32:  // fall through
    case FieldDescriptor::TYPE_UINT32:
      switch(int32_type) {
	case FieldOptions::INT_16:       return "SMALLINT UNSIGNED";
	case FieldOptions::INT_8:        return "TINYINT UNSIGNED";
	case FieldOptions::INT_32:       // fall through
        default:                         return "INT UNSIGNED";
      }
    case FieldDescriptor::TYPE_ENUM:     return "INT";

    case FieldDescriptor::TYPE_MESSAGE:  // fallthrough to assert(0)
    case FieldDescriptor::TYPE_GROUP:    // fallthrough to assert(0)
    default:
      break;
  }

  assert(0);
}

std::string SQLTransceiver::
sql_create_table(const SQLTable* t, bool if_not_exists) 
{
  std::ostringstream sql;
  std::vector<const SQLTableField*> keys;
  sql << "CREATE TABLE ";
  if(if_not_exists)sql << "IF NOT EXISTS ";
  sql << sql_table_name(t->table_name) << " ( \n";
  for ( auto f : t->fields )
  {
    if(f->is_key())keys.push_back(f);
    sql << "  " << sql_field_name(f->field_name) << ' ' << sql_type(f->field_d);
    if(f != t->fields.back() or !keys.empty())sql << ',';
    if(f->field_d)
    {
      const google::protobuf::FieldOptions* fopt { &f->field_d->options() };
      if(fopt->HasExtension(CFO) and
         (!fopt->GetExtension(CFO).desc().empty() or
          !fopt->GetExtension(CFO).units().empty()))
      {
        sql << " --";
        if(!fopt->GetExtension(CFO).desc().empty())
          sql << ' ' << fopt->GetExtension(CFO).desc();
        if(!fopt->GetExtension(CFO).units().empty())
          sql << " [" << fopt->GetExtension(CFO).units() << ']';
      }
    }
    sql << '\n';
  }
  if(!keys.empty())
  {
    sql << "  PRIMARY KEY (\n";
    for(auto f : keys)
    {
      sql << "    " << sql_field_name(f->field_name);
      if(f != keys.back())sql << ',';
      sql << '\n';
    }
    sql << "  )\n";
  }
  sql << ')';
  return sql.str();
}

std::string SQLTransceiver::
sql_insert(const SQLTable* t)
{
  std::ostringstream sql;
  std::vector<const SQLTableField*> keys;
  sql << "INSERT INTO " << sql_table_name(t->table_name) << " (\n";
  for ( auto f : t->fields )
  {
    sql << "  " << sql_field_name(f->field_name);
    if(f != t->fields.back() or !keys.empty())sql << ',';
    sql << '\n';
  }
  sql << ") VALUES (\n";
  for ( auto f : t->fields )
  {
    sql << "  ?";
    if(f != t->fields.back() or !keys.empty())sql << ',';
    sql << '\n';
  }
  sql << ')';
  return sql.str();
}

bool SQLTransceiver::prepare_insert_statements(SQLTable* t)
{
  iterate_over_tables(t, [this](SQLTable* t) {
      t->stmt = new Statement(sql_insert(t));
    });
  return true;
}

bool SQLTransceiver::
insert_tables(SQLTable* t, const google::protobuf::Message* m_data,
              const google::protobuf::Message* m_key,
              uint64_t parent_oid, uint64_t loop_id,
              bool write_sql_to_log)
{
  unsigned ifield = 0;
  for(auto f : t->fields)
  {
    const google::protobuf::Message* m = m_data;
    if(f->field_type == SQLTableField::KEY_USER_SUPPLIED)m = m_key;
    assert(m);
    const google::protobuf::Reflection* r = m->GetReflection();

    if(f->field_origin == f)
    {
      if(f->field_d != nullptr)
      {
        for(auto d : f->field_d_path) {
          if(!r->HasField(*m, d))goto next_field;
          m = &r->GetMessage(*m, d);
          r = m->GetReflection();
        }
        
        f->data =
            static_cast<void*>(const_cast<google::protobuf::Message*>(m));
      }
      else
      {
        switch(f->field_type)
        {
          case SQLTableField::KEY_OID:
            f->data = &parent_oid;
            break;
          case SQLTableField::KEY_LOOP_ID:
            f->data = &loop_id;
            break;
          case SQLTableField::KEY_PROTO_DEFINED:  // handled in if clause above
          case SQLTableField::KEY_USER_SUPPLIED:  // handled in if clause above
          case SQLTableField::KEY_INHERITED:      // handled earlier in tree
          case SQLTableField::KEY_MAP_KEY:        // handled in if clause above
          case SQLTableField::POD:                // handled in if clause above
            assert(0);
            break;
        }
      }
    }

    f = f->field_origin;
    assert(f->data);

    if(f->field_d == nullptr)
      t->stmt->bind_uint64(ifield, *static_cast<uint64_t*>(f->data));
    else if(f->field_d->is_repeated())
      t->stmt->bind_repeated_field(ifield, loop_id, 
               static_cast<google::protobuf::Message*>(f->data), f->field_d);
    else
      t->stmt->bind_field(ifield,
               static_cast<google::protobuf::Message*>(f->data), f->field_d);

 next_field:
    ifield++;
  }

  bool good = true;
  good &= execute_one_insert_statement(t, parent_oid);
  if(write_sql_to_log)
    LOG(INFO) << t->stmt->bound_sql() << " -> " << parent_oid;

  for(auto st : t->sub_tables)
  {
    const google::protobuf::Message* m = m_data;
    const google::protobuf::Reflection* r = m->GetReflection();
    for(auto d : st->parent_field_d_path) {
      if(!r->HasField(*m, d))goto next_sub_table;
      m = &r->GetMessage(*m, d);
      r = m->GetReflection();      
    }
    
    if(st->parent_field_d->is_repeated())
    {
      uint64_t nloop = r->FieldSize(*m, st->parent_field_d);
      for(uint64_t iloop = 0; iloop < nloop; iloop++)
      {
        const google::protobuf::Message* mi = m;
        if(st->parent_field_d->type() == FieldDescriptor::TYPE_MESSAGE) {
          mi = &r->GetRepeatedMessage(*m, st->parent_field_d, iloop);
          r = m->GetReflection(); }
        good &= insert_tables(st, mi, nullptr, parent_oid, iloop,
                              write_sql_to_log);
      }
    }
    else
    {
      assert(st->parent_field_d->type() == FieldDescriptor::TYPE_MESSAGE);
      if(!r->HasField(*m, st->parent_field_d))goto next_sub_table;
      m = &r->GetMessage(*m, st->parent_field_d);
      r = m->GetReflection();
      good &= insert_tables(st, m, nullptr, parent_oid, 0, write_sql_to_log);
    }
 next_sub_table:
    ;
  }
  
  return good;
}

bool SQLTransceiver::execute_one_insert_statement(SQLTable* t, uint64_t& oid)
{
  auto L = LOG(INFO);
  Statement* s { static_cast<Statement*>(t->stmt) };
  return true;
}

bool SQLTransceiver::finalize_insert_statements(SQLTable* t)
{
  iterate_over_tables(t, [this](SQLTable* t) {
      delete t->stmt; t->stmt = nullptr; });
  return true;
}

// ============================================================================
// ============================================================================
//
// SQLTransceiver::Statement
//
// ============================================================================
// ============================================================================

SQLTransceiver::Statement::Statement(const std::string& sql):
    sql_(sql), bound_values_()
{
  // nothing to see here
}

SQLTransceiver::Statement::~Statement()
{
  // nothing to see here
}

std::string SQLTransceiver::Statement::bound_sql() const
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

unsigned SQLTransceiver::Statement::num_columns()
{
  return 0;
}

bool SQLTransceiver::Statement::is_initialized()
{
  return true;
}

void SQLTransceiver::Statement::
error_codes(int& error_num, std::string& error_msg)
{
  error_num = 0;
  error_msg = "";
}

void SQLTransceiver::Statement::reset()
{
  bound_values_.clear();
}

bool SQLTransceiver::Statement::step()
{
  // nothing to see here
}

bool SQLTransceiver::Statement::
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

bool SQLTransceiver::Statement::
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

bool SQLTransceiver::Statement::bind_null(unsigned ifield)
{
  return SQLTransceiver::Statement::bind_string(ifield, "");
}

bool SQLTransceiver::Statement::bind_int64(unsigned ifield, int64_t value)
{
  return SQLTransceiver::Statement::bind_string(ifield, std::to_string(value));
}

bool SQLTransceiver::Statement::bind_int32(unsigned ifield, int32_t value)
{
  return SQLTransceiver::Statement::bind_string(ifield, std::to_string(value));
}

bool SQLTransceiver::Statement::bind_int16(unsigned ifield, int16_t value)
{
  return SQLTransceiver::Statement::bind_string(ifield, std::to_string(value));
}

bool SQLTransceiver::Statement::bind_int8(unsigned ifield, int8_t value)
{
  return SQLTransceiver::Statement::bind_string(ifield, std::to_string(value));
}

bool SQLTransceiver::Statement::bind_uint64(unsigned ifield, uint64_t value)
{
  return SQLTransceiver::Statement::bind_string(ifield, std::to_string(value));
}

bool SQLTransceiver::Statement::bind_uint32(unsigned ifield, uint32_t value)
{
  return SQLTransceiver::Statement::bind_string(ifield, std::to_string(value));
}

bool SQLTransceiver::Statement::bind_uint16(unsigned ifield, uint16_t value)
{
  return SQLTransceiver::Statement::bind_string(ifield, std::to_string(value));
}

bool SQLTransceiver::Statement::bind_uint8(unsigned ifield, uint8_t value)
{
  return SQLTransceiver::Statement::bind_string(ifield, std::to_string(value));
}

bool SQLTransceiver::Statement::bind_float(unsigned ifield, float value)
{
  return SQLTransceiver::Statement::bind_string(ifield, std::to_string(value));
}

bool SQLTransceiver::Statement::bind_double(unsigned ifield, double value)
{
  return SQLTransceiver::Statement::bind_string(ifield, std::to_string(value));
}

bool SQLTransceiver::Statement::bind_bool(unsigned ifield, bool value)
{
  return SQLTransceiver::Statement::bind_string(ifield, std::to_string(value));
}

bool SQLTransceiver::Statement::
bind_string(unsigned ifield, const std::string& value)
{
  if(ifield >= bound_values_.size())bound_values_.resize(ifield+1);
  bound_values_[ifield] = value;
  return true;
}

bool SQLTransceiver::Statement::
bind_bytes(unsigned ifield, const std::string& value)
{
  return SQLTransceiver::Statement::bind_string(ifield, value);
}
