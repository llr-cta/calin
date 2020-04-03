/*

   calin/io/sql_serializer.cpp -- Stephen Fegan -- 2020-03-31

   Base class for reading and writing protobuf structures to SQL databases

   Reimplementation of calin/io/sql_transceiver.cpp -- Stephen Fegan -- 2015-09-08

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

#include <memory>
#include <deque>
#include <algorithm>
#include <stdexcept>

#include <google/protobuf/io/printer.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include <calin.pb.h>
#include <io/sql_serializer.pb.h>
#include <io/sql_serializer.hpp>
#include <util/log.hpp>
#include <util/string.hpp>

using namespace calin::util::log;
using namespace calin::util::string;
using namespace calin::io::sql_serializer;
using namespace google::protobuf;

namespace {
  std::string desc_units_comment(const std::string& desc, const std::string& units)
  {
   std::string comment = desc;
   if(not units.empty()) {
     if(comment.empty()) {
       comment += "[";
     } else {
       comment += " [";
     }
     comment += units;
     comment += "]";
   }
   return comment;
  }
}

std::string SQLTableField::field_comment() const {
  return desc_units_comment(field_desc, field_units);
}

std::string SQLTable::table_comment() const {
  return desc_units_comment(table_desc, table_units);
}

SQLSerializer::~SQLSerializer()
{
  for(auto itable : db_tables_)delete itable.second;
  for(auto ifield : db_table_fields_)delete ifield.second;
  for(auto& ischema : schema_) {
    for(auto& it : ischema.second) {
      delete it.second;
    }
  }
}

calin::io::sql_serializer::SQLTable* SQLSerializer::
make_sqltable_tree(const std::string& table_name,
  const google::protobuf::Descriptor* d, const std::string& instance_desc,
  bool propagate_keys)
{
  auto* t = r_make_sqltable_tree(table_name, d,
    /* parent_table= */ nullptr, /* parent_field_d= */ nullptr,
    /* table_desc = */ instance_desc, /* table_units= */ "");
  if(propagate_keys) {
    r_propagate_keys(t, { });
  }
  test_sqltable_tree_db_presence(t);
  return t;
}

bool SQLSerializer::create_or_extend_tables(const std::string& table_name,
  const google::protobuf::Descriptor* d, const std::string& instance_desc)
{
  std::unique_ptr<SQLTable> t { make_sqltable_tree(table_name, d, instance_desc) };
  return do_create_or_extend_tables(table_name, t.get());
}

bool SQLSerializer::insert(const std::string& table_name, uint64_t& oid,
  const google::protobuf::Message* m)
{
  const google::protobuf::Descriptor* d = m->GetDescriptor();
  SQLTable* t = schema_[table_name][d];
  if(t == nullptr) {
    t = make_sqltable_tree(table_name, d);
    schema_[table_name][d] = t;
  }

  if(not t->db_all_tree_fields_present)
  {
    for(auto& it : schema_[table_name]) {
      if(it.second != t) {
        delete it.second;
        it.second = nullptr;
      }
    }

    do_create_or_extend_tables(table_name, t);

    test_sqltable_tree_db_presence(t);
    if(not t->db_all_tree_fields_present)
      throw std::runtime_error("SQLSerializer::insert: not all fields preset in database");
  }

  bool success = true;
#if 0
  if(t.insert_stmt == nullptr) {
    unsigned nfield_before = db_table_fields_.size();

  }

  bool success = true;
  t->iterate_over_tables([this,&success](SQLTable* it) {
    if(success)
    {
      if(it->insert_stmt == nullptr) {
        it->insert_stmt = prepare_statement(sql_insert(it));
      }
      if(not it->insert_stmt->is_initialized())
      {
        LOG(ERROR) << "SQL error preparing INSERT: " << it->insert_stmt->error_message() << '\n'
          << "SQL: " << it->insert_stmt->sql();
        success = false;
      }
    }
  });

  if(!success)return success;

  begin_transaction();
//  success = r_exec_insert(t.get(), m, oid, 0, 0, false);
  if(!success)
  {
    rollback_transaction();
    return success;
  }
  commit_transaction();
#endif
  return success;
}

calin::io::sql_serializer::SQLTable* SQLSerializer::
r_make_sqltable_tree(const std::string& table_name, const google::protobuf::Descriptor* d,
  SQLTable* parent_table, const google::protobuf::FieldDescriptor* parent_field_d,
  const std::string& table_desc, const std::string& table_units)
{
  SQLTable* t { new SQLTable };
  t->table_d                     = d;
  t->parent_table                = parent_table;
  t->table_name                  = table_name;
  t->parent_field_d              = parent_field_d;
  t->table_desc                  = table_desc;
  t->table_units                 = table_units;

  for(int ioneof = 0; ioneof<d->oneof_decl_count(); ioneof++)
  {
    const OneofDescriptor* oo { d->oneof_decl(ioneof) };
    SQLTableField* tf { new SQLTableField };
    tf->table             = t;
    tf->field_origin      = tf;
    tf->field_type        = SQLTableField::POD;
    tf->field_name        = oo->name();
    tf->oneof_d           = oo;
    tf->field_desc        = "One-of selector";
    t->fields.push_back(tf);
  }

  for(int ifield = 0; ifield<d->field_count(); ifield++)
  {
    const FieldDescriptor* f { d->field(ifield) };
    const google::protobuf::FieldOptions* fopt { &f->options() };

    std::string field_desc;
    std::string field_units;
    if(fopt->HasExtension(CFO))
    {
      if((fopt->GetExtension(CFO).dont_store() or fopt->GetExtension(CFO).sql().dont_store()))
        continue;
      field_desc = fopt->GetExtension(CFO).desc();
      field_units = fopt->GetExtension(CFO).units();
    }

    SQLTable* sub_table { nullptr };

    if(f->type()==FieldDescriptor::TYPE_MESSAGE)
    {
      sub_table = r_make_sqltable_tree(sub_name(table_name, f->name()),
        f->message_type(), t, f, field_desc, field_units);

      if(parent_field_d and parent_field_d->is_map())
      {
        // Message can be inlined - so move its fields into primary table

        for(auto ifield : sub_table->fields)
        {
          ifield->table           = t;
          ifield->field_name      = sub_name(f->name(), ifield->field_name);
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
      sub_table->parent_table                = t;
      sub_table->table_name                  = sub_name(table_name,f->name());
      sub_table->parent_field_d              = f;
      sub_table->table_desc                  = field_desc;
      sub_table->table_units                 = field_units;

      SQLTableField* tf { new SQLTableField };
      tf->table             = sub_table;
      tf->field_origin      = tf;
      tf->field_type        = SQLTableField::POD;
      tf->field_name        = f->name();
      tf->field_d           = f;
      tf->field_desc        = field_desc;
      tf->field_units       = field_units;

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
      tf->field_desc        = field_desc;
      tf->field_units       = field_units;

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
      tf->db_field_present  = db_table_fields_.count(sub_name(t->table_name, tf->field_name));

      sub_table->fields.insert(sub_table->fields.begin(), tf);
    }

    t->sub_tables.insert(t->sub_tables.end(), sub_table);
  }
  return t;
}

void SQLSerializer::r_propagate_keys(SQLTable* t, std::vector<const SQLTableField*> keys)
{
  for(unsigned ikey=0;ikey<keys.size();ikey++)
  {
    const SQLTableField* key = keys[keys.size()-ikey-1];
    SQLTableField* f = new SQLTableField(*key);
    f->field_type = SQLTableField::KEY_INHERITED;
    f->table      = t; // WARNING - is this correct ?
    t->fields.insert(t->fields.begin(), f);
  }

  if(t->parent_table != nullptr and keys.empty())
  {
    // We have a parent, but no inherited keys, so add parent OID as key
    SQLTableField* f { new SQLTableField };
    f->table          = t;
    f->field_origin   = f;
    f->field_type     = SQLTableField::KEY_PARENT_OID;
    f->field_name     = sub_name(t->parent_table->table_name,
                                 sql_oid_column_name());
    f->field_d        = nullptr;
    t->fields.insert(t->fields.begin(), f);
  }

  for(unsigned ikey=keys.size();ikey<t->fields.size();ikey++)
    if(t->fields[ikey]->is_key())keys.push_back(t->fields[ikey]);

  for(auto it : t->sub_tables)
    r_propagate_keys(it, keys);
}

void SQLSerializer::test_sqltable_tree_db_presence(SQLTable* t)
{
  t->iterate_over_tables([this](SQLTable* it) {
    it->db_all_table_fields_present = true;
    it->db_all_tree_fields_present = true;
    for(auto* f : it->fields) {
      f->db_field_present = this->db_table_fields_.count(sub_name(it->table_name, f->field_name));
      it->db_all_table_fields_present &= f->db_field_present;
    }
    it->db_all_tree_fields_present = it->db_all_table_fields_present;
  });
  t->reverse_iterate_over_tables([](SQLTable* it) {
    if(it->parent_table) {
      it->parent_table->db_all_tree_fields_present &= it->db_all_table_fields_present;
    }
  });
}

bool SQLSerializer::do_create_or_extend_tables(const std::string& table_name, SQLTable* t)
{
  bool success = true;

  std::vector<SQLTable*> new_tables;
  std::vector<SQLTableField*> new_table_fields;

  begin_transaction();
  t->iterate_over_tables([this, &success, &new_tables, &new_table_fields](SQLTable* itable) {
    if(success) {
      if(this->db_tables_.find(itable->table_name) == this->db_tables_.end()) {
        auto* stmt = prepare_statement(sql_create_table(itable));
        if(not stmt->is_initialized()) {
          LOG(ERROR) << "SQL error preparing CREATE TABLE: " << stmt->error_message() << '\n'
                     << "SQL: " << stmt->sql();
          success = false;
        } else {
          success = execute_one_no_data_statement(stmt);
          if(success) {
            new_tables.push_back(itable);
          }
        }
        delete stmt;
      } else {
        for(auto ifield : itable->fields) {
          std::string fqdn = sub_name(itable->table_name, ifield->field_name);
          if(this->db_table_fields_.find(fqdn) == this->db_table_fields_.end()) {
            auto* stmt = prepare_statement(sql_add_field_to_table(ifield));
            if(not stmt->is_initialized()) {
              LOG(ERROR) << "SQL error preparing EXTEND TABLE: " << stmt->error_message() << '\n'
                         << "SQL: " << stmt->sql();
              success = false;
            } else {
              success = execute_one_no_data_statement(stmt);
              if(success) {
                new_table_fields.push_back(ifield);
              }
            }
            delete stmt;
          }
        }
      }
    }
  });

  if(not success)
  {
    rollback_transaction();
    return success;
  }

  success = true;
  for(SQLTable* it : new_tables) {
    if(success) {
      std::string sql = sql_create_index(it);
      if(sql.empty()) continue;
      auto* stmt = prepare_statement(sql);
      if(not stmt->is_initialized()) {
        LOG(ERROR) << "SQL error preparing CREATE INDEX: " << stmt->error_message() << '\n'
                   << "SQL: " << stmt->sql();
        success = false;
      } else {
        success = execute_one_no_data_statement(stmt);
        delete stmt;
      }
    }
  }

  if(not success)
  {
    rollback_transaction();
    return success;
  }

  std::vector<calin::ix::io::sql_serializer::SQLTable*> new_table_protos;
  std::vector<calin::ix::io::sql_serializer::SQLTableField*> new_field_protos;
  for(auto t : new_tables) {
    new_table_protos.push_back(table_as_proto(t));
    for(auto f : t->fields) {
      new_field_protos.push_back(field_as_proto(f));
    }
  }
  for(auto f : new_table_fields) {
    new_field_protos.push_back(field_as_proto(f));
  }

  // Try to insert into the intenal fields here

  commit_transaction();

  for(auto pt : new_table_protos) {
    db_tables_[pt->table_name()] = pt;
  }
  for(auto pf : new_field_protos) {
    std::string fqdn = sub_name(pf->table_name(), pf->field_name());
    db_table_fields_[fqdn] = pf;
  }

  return true;
}

SQLStatement* SQLSerializer::prepare_statement(const std::string& sql)
{
  return new SQLStatement(sql);
}

bool SQLSerializer::begin_transaction()
{
  return true;
}

bool SQLSerializer::commit_transaction()
{
  return true;
}

bool SQLSerializer::rollback_transaction()
{
  return true;
}

bool SQLSerializer::execute_one_no_data_statement(SQLStatement* stmt, bool ignore_errors)
{
  bool good = true;
  if(stmt != nullptr)
  {
    SQLStatement::StepStatus status = stmt->step();
    if(write_sql_to_log_)LOG(INFO) << stmt->bound_sql();
    switch(status)
    {
    case SQLStatement::ERROR:
      good = false;
      if(not ignore_errors)
      {
        LOG(ERROR) << "SQL statement returned error: "
                   << stmt->error_message() << '\n'
                   << "SQL: " << stmt->sql();
        return good;
      }
      break;
    case SQLStatement::OK_NO_DATA:
      // this is what we expect
      break;
    case SQLStatement::OK_HAS_DATA:
      LOG(ERROR) << "Simple SQL statement returned data" << '\n'
                 << "SQL: " << stmt->sql();
      return false;
    }
    stmt->reset();
  }
  return good;
}

// =============================================================================
// =============================================================================
//
// Overridable member functions to create SQL strings
//
// =============================================================================
// =============================================================================

std::string SQLSerializer::
sub_name(const std::string& parent_name,const std::string& name)
{
  std::string new_name { parent_name };
  new_name += ".";
  new_name += name;
  return new_name;
}

std::string SQLSerializer::sql_oid_column_name()
{
  return "_ROWID_";
}

std::string SQLSerializer::sql_table_name(const std::string& name)
{
  std::string new_name { '`' };
  new_name += name;
  new_name += '`';
  return new_name;
}

std::string SQLSerializer::sql_field_name(const std::string& name)
{
  std::string new_name { '`' };
  new_name += name;
  new_name += '`';
  return new_name;
}

std::string SQLSerializer::sql_type(const google::protobuf::FieldDescriptor* d)
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
  return "UNKNOWN";
}

std::string SQLSerializer::sql_insert_field_spec(const SQLTableField* f)
{
  return "?";
}

std::string SQLSerializer::sql_comment(const std::string& comment,
  unsigned first_line_indent, unsigned multi_line_indent,
  bool newline_before_multi_line)
{
  if(comment.empty())return std::string();
  std::vector<std::string> bits = split(comment,'\n');
  if(bits.size() == 1)
    return std::string(first_line_indent, ' ')+std::string("-- ")+bits.front();
  std::string out;
  if(newline_before_multi_line)out = "\n";
  for(const auto& iline : bits)
  {
    if(&iline==&bits.front())
      out += std::string(newline_before_multi_line?multi_line_indent:
                         first_line_indent, ' ');
    else
      out += std::string("\n") + std::string(multi_line_indent, ' ');
    out += "-- ";
    out += iline;
  }
  return out;
}

std::string SQLSerializer::sql_create_table(const SQLTable* t)
{
  std::ostringstream sql;
  sql << "CREATE TABLE ";
  sql << sql_table_name(t->table_name) << " ( \n";
  std::string table_comment = t->table_comment();
  if(not table_comment.empty()) {
    sql << sql_comment(table_comment) << '\n';
  }
  sql << "  " << sql_field_name(sql_oid_column_name()) << " INTEGER PRIMARY KEY";
  if(!t->fields.empty())sql << ',';
  sql << '\n';
  for ( auto f : t->fields )
  {
    sql << "  " << sql_field_name(f->field_name) << ' ' << sql_type(f->field_d);
    if(f != t->fields.back())sql << ',';
    std::string field_comment = f->field_comment();
    if(not field_comment.empty()) {
      sql << sql_comment(field_comment,1,4,true);
    }
    sql << '\n';
  }
  sql << ")";
  return sql.str();
}

std::string SQLSerializer::sql_add_field_to_table(const SQLTableField* f)
{
  std::ostringstream sql;
  sql << "ALTER TABLE " << sql_table_name(f->table->table_name)
    << " ADD COLUMN " << sql_field_name(f->field_name) << ' ' << sql_type(f->field_d);
  std::string field_comment = f->field_comment();
  if(not field_comment.empty()) {
    sql << sql_comment(field_comment,1,4,true);
  }
  return sql.str();
}

calin::ix::io::sql_serializer::SQLTableAndFieldCollection*
SQLSerializer::sqltable_tree_as_proto(const SQLTable* t)
{
  auto* proto = new calin::ix::io::sql_serializer::SQLTableAndFieldCollection;
  t->iterate_over_tables([this,proto](const SQLTable* it) {
    proto->mutable_tables()->AddAllocated(this->table_as_proto(it));
    for(auto f : it->fields) {
      proto->mutable_fields()->AddAllocated(this->field_as_proto(f));
    }
  });
  return proto;
}

calin::ix::io::sql_serializer::SQLTable* SQLSerializer::table_as_proto(const SQLTable* t)
{
  auto* proto = new calin::ix::io::sql_serializer::SQLTable;
  auto* root_t = t;
  while(root_t->parent_table != nullptr)root_t = root_t->parent_table;
  proto->set_base_name(root_t->table_name);
  proto->set_table_name(t->table_name);
  proto->set_sql_table_name(sql_table_name(t->table_name));
  proto->set_description(t->table_desc);
  proto->set_units(t->table_units);
  if(t->table_d) {
    proto->set_proto_message_type(t->table_d->full_name());
  }
  return proto;
}

calin::ix::io::sql_serializer::SQLTableField* SQLSerializer::field_as_proto(const SQLTableField* f)
{
  auto* proto = new calin::ix::io::sql_serializer::SQLTableField;
  auto* root_t = f->table;
  while(root_t->parent_table != nullptr)root_t = root_t->parent_table;
  proto->set_base_name(root_t->table_name);
  proto->set_table_name(f->table->table_name);
  proto->set_field_name(f->field_name);
  proto->set_sql_table_name(sql_table_name(f->table->table_name));
  proto->set_sql_field_name(sql_table_name(f->field_name));
  proto->set_description(f->field_desc);
  proto->set_units(f->field_units);
  if(f->field_d)
  {
    proto->set_proto_message_type(f->field_d->containing_type()->full_name());
    proto->set_proto_field_name(f->field_d->name());
    proto->set_proto_field_number(f->field_d->number());
  }
  else if(f->oneof_d)
  {
    proto->set_proto_message_type(f->oneof_d->containing_type()->full_name());
    proto->set_proto_field_name(f->oneof_d->name());
  }
  return proto;
}

std::string SQLSerializer::sql_create_index(const SQLTable* t)
{
  std::vector<const SQLTableField*> keys;
  for ( auto f : t->fields )
    if(f->is_key())keys.push_back(f);
  if(keys.empty())return {};
  std::ostringstream sql;
  sql << "CREATE UNIQUE INDEX ";
  sql << sql_table_name(t->table_name + "$index") << "\n";
  sql << "  ON " << sql_table_name(t->table_name) << " (\n";
  for(auto f : keys)
  {
    sql << "    " << sql_field_name(f->field_name);
    if(f != keys.back())sql << ',';
    sql << '\n';
  }
  sql << "  )\n";
  return sql.str();
}

std::string SQLSerializer::sql_insert(const SQLTable* t)
{
  std::ostringstream sql;
  sql << "INSERT INTO " << sql_table_name(t->table_name);
  if(t->fields.empty()) {
    sql << " DEFAULT VALUES\n";
    return sql.str();
  }
  sql << " (\n";
  for ( auto f : t->fields )
  {
    sql << "  " << sql_field_name(f->field_name);
    if(f != t->fields.back())sql << ',';
    sql << '\n';
  }
  sql << ") VALUES (\n";
  for ( auto f : t->fields )
  {
    sql << "  " << sql_insert_field_spec(f);
    if(f != t->fields.back())sql << ',';
    sql << '\n';
  }
  sql << ')';
  return sql.str();
}
