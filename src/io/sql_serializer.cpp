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
  return do_create_or_extend_tables(t.get());
}

bool SQLSerializer::insert(const std::string& table_name, uint64_t& oid,
  const google::protobuf::Message* m, bool allow_mismatching_type)
{
  const google::protobuf::Descriptor* d = m->GetDescriptor();
  if(not allow_mismatching_type and db_tables_.count(table_name) and
      db_tables_.at(table_name)->proto_message_type() != d->full_name()) {
    throw std::runtime_error("SQLSerializer::insert: mistmatching types in table "
      + table_name + " (" + db_tables_.at(table_name)->proto_message_type()
      + " != " + d->full_name() + "). Retry with allow_mismatching_type=true to ignore.");
  }

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

    if(not do_create_or_extend_tables(t)) {
      return false;
    }

    test_sqltable_tree_db_presence(t);
    if(not t->db_all_tree_fields_present)
      throw std::runtime_error("SQLSerializer::insert: not all fields preset in database");
  }

  oid = 0;
  bool success = true;

  t->iterate_over_tables([this,&success](SQLTable* it) {
    if(success)
    {
      if(it->stmt_insert == nullptr) {
        it->stmt_insert = prepare_statement(sql_insert(it));
      }
      if(not it->stmt_insert->is_initialized())
      {
        LOG(ERROR) << "SQL error preparing INSERT: " << it->stmt_insert->error_message() << '\n'
          << "SQL: " << it->stmt_insert->sql();
        success = false;
      }
    }
  });

  if(not success)
  {
    return success;
  }

  begin_transaction();

  try {
    success = r_exec_insert(t, m, oid, 0, 0, false);
  } catch(...) {
    rollback_transaction();
    throw;
  }

  if(not success)
  {
    rollback_transaction();
    return success;
  }

  commit_transaction();
  return success;
}

bool SQLSerializer::retrieve_by_oid(const std::string& table_name, uint64_t oid,
  google::protobuf::Message* m)
{
  const google::protobuf::Descriptor* d = m->GetDescriptor();
  SQLTable* t = schema_[table_name][d];
  if(t == nullptr) {
    t = make_sqltable_tree(table_name, d);
    schema_[table_name][d] = t;
  }

  bool success = true;

  t->iterate_over_tables([this,&success](SQLTable* it) {
    if(success)
    {
      if(it->db_table_present) {
        if(it->stmt_select_oid == nullptr) {
          it->stmt_select_oid = prepare_statement(sql_select_where_oid_equals(it));
        }
        if(it->stmt_select_oid != nullptr and not it->stmt_select_oid->is_initialized())
        {
          LOG(ERROR) << "SQL error preparing SELECT: " << it->stmt_select_oid->error_message() << '\n'
            << "SQL: " << it->stmt_select_oid->sql();
          success = false;
        }
      }
    }
  });

  if(not success)
  {
    return success;
  }

  begin_transaction();

  try {
    t->iterate_over_tables([this,oid,m,&success](SQLTable* it) {
      if(success)
      {
        if(it->stmt_select_oid != nullptr) {
          success = exec_select_by_oid(it, oid, m, /* ignore_errors= */ false);
        }
      }
    });
  } catch(...) {
    rollback_transaction();
    throw;
  }

  if(not success)
  {
    rollback_transaction();
    return success;
  }

  commit_transaction();
  return success;
}

uint64_t SQLSerializer::count_entries_in_table(const std::string& table_name)
{
  std::unique_ptr<SQLStatement> stmt {
    prepare_statement(sql_count_entries(table_name)) };
  if(!stmt->is_initialized()) {
    LOG(ERROR) << "SQL error preparing COUNT : "
               << stmt->error_message() << '\n'
               << "SQL: " << stmt->sql();
    throw std::runtime_error("Could not prepare SQL count statement");
  }

  auto status = stmt->step();
  if(status == SQLStatement::ERROR) {
    LOG(ERROR) << "SQL error executing COUNT : "
               << stmt->error_message() << '\n'
               << "SQL: " << stmt->sql();
    throw std::runtime_error("SQL count returned error");
  } else if (status == SQLStatement::OK_NO_DATA) {
    LOG(ERROR) << "SQL error executing COUNT : "
               << stmt->error_message() << '\n'
               << "SQL: " << stmt->sql();
    throw std::logic_error("SQL count returned no data");
  }

  bool good = true;
  uint64_t entries = stmt->extract_uint64(0, &good);
  if(!good) {
    LOG(ERROR) << "Could not extract value from SQL count\n"
               << "SQL: " << stmt->sql();
    throw std::logic_error("SQL count could not be extracted");
  }

  if(stmt->step() != SQLStatement::OK_NO_DATA)
    throw std::runtime_error("SQLSerializer::count_entries_in_table: unexpected data returned");

  return entries;
}

std::vector<uint64_t> SQLSerializer::retrieve_all_oids(const std::string& table_name)
{
  std::unique_ptr<SQLStatement> stmt {
    prepare_statement(sql_select_oids(table_name)) };

  if(!stmt->is_initialized()) {
    LOG(ERROR) << "SQL error preparing SELECT OID : "
               << stmt->error_message() << '\n'
               << "SQL: " << stmt->sql();
    throw std::runtime_error("Could not prepare SQL select OID statement");
  }

  std::vector<uint64_t> oids;
  SQLStatement::StepStatus status = SQLStatement::ERROR;
  for(status = stmt->step(); status == SQLStatement::OK_HAS_DATA;
    status = stmt->step())
  {
    bool good = true;
    uint64_t oid = stmt->extract_uint64(0, &good);
    if(!good) {
      LOG(ERROR) << "Could not extract OID from SQL select\n"
                 << "SQL: " << stmt->sql();
      throw std::logic_error("SQL OID could not be extracted");
    }
    oids.emplace_back(oid);
  }

  if(status == SQLStatement::ERROR) {
    LOG(ERROR) << "SQL error executing SELECT OID : "
               << stmt->error_message() << '\n'
               << "SQL: " << stmt->sql();
    throw std::runtime_error("SQL SELECT OID returned error");
  }

  return oids;
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
  if(parent_table == nullptr) {
    t->root_table =              t;
  } else {
    t->root_table =              parent_table->root_table;
    t->root_field_d_path =       parent_table->root_field_d_path;
    if(parent_table->parent_field_d) {
      t->root_field_d_path.emplace_back(parent_table->parent_field_d);
    }
  }

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
    }
    else if(f->is_repeated())
    {
      sub_table = new SQLTable;
      sub_table->parent_table                = t;
      sub_table->table_name                  = sub_name(table_name,f->name());
      sub_table->parent_field_d              = f;
      sub_table->table_desc                  = field_desc;
      sub_table->table_units                 = field_units;

      sub_table->root_table                  = t->root_table;
      sub_table->root_field_d_path           = t->root_field_d_path;
      if(t->parent_field_d) {
        sub_table->root_field_d_path.emplace_back(t->parent_field_d);
      }

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

    if(f->is_repeated())
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
    f->field_name     = sub_name(t->parent_table->table_name, sql_oid_column_name());
    f->is_root_oid    = (t->parent_table == t->root_table);
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
    it->db_table_present = this->db_tables_.count(it->table_name);
    it->db_all_table_fields_present = it->db_table_present;
    it->db_all_tree_fields_present = it->db_table_present;
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

void SQLSerializer::force_sqltable_tree_db_presence(SQLTable* t)
{
  t->iterate_over_tables([](SQLTable* it) {
    it->db_table_present = true;
    it->db_all_table_fields_present = true;
    it->db_all_tree_fields_present = true;
    for(auto* f : it->fields) {
      f->db_field_present = true;
    }
  });
}

bool SQLSerializer::do_create_or_extend_tables(SQLTable* t,
  bool write_new_tables_and_fields_to_db)
{
  bool success = true;

  std::vector<SQLTable*> new_tables;
  std::vector<SQLTableField*> new_table_fields;

  begin_transaction();

  std::vector<calin::ix::io::sql_serializer::SQLTable*> new_table_protos;
  std::vector<calin::ix::io::sql_serializer::SQLTableField*> new_field_protos;

  try {
    t->iterate_over_tables([this, &success, &new_tables, &new_table_fields](SQLTable* itable) {
      if(success) {
        if(this->db_tables_.find(itable->table_name) == this->db_tables_.end()) {
          std::unique_ptr<SQLStatement> stmt { prepare_statement(sql_create_table(itable)) };
          if(not stmt->is_initialized()) {
            LOG(ERROR) << "SQL error preparing CREATE TABLE: " << stmt->error_message() << '\n'
                       << "SQL: " << stmt->sql();
            success = false;
          } else {
            success = execute_one_no_data_statement(stmt.get());
            if(success) {
              new_tables.push_back(itable);
            }
          }
        } else {
          for(auto ifield : itable->fields) {
            std::string fqdn = sub_name(itable->table_name, ifield->field_name);
            if(this->db_table_fields_.find(fqdn) == this->db_table_fields_.end()) {
              std::unique_ptr<SQLStatement> stmt { prepare_statement(sql_add_field_to_table(ifield))  };
              if(not stmt->is_initialized()) {
                LOG(ERROR) << "SQL error preparing EXTEND TABLE: " << stmt->error_message() << '\n'
                           << "SQL: " << stmt->sql();
                success = false;
              } else {
                success = execute_one_no_data_statement(stmt.get());
                if(success) {
                  new_table_fields.push_back(ifield);
                }
              }
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

    for(auto t : new_tables) {
      new_table_protos.push_back(table_as_proto(t));
      for(auto f : t->fields) {
        new_field_protos.push_back(field_as_proto(f));
      }
    }
    for(auto f : new_table_fields) {
      new_field_protos.push_back(field_as_proto(f));
    }

    if(write_new_tables_and_fields_to_db) {
      for(auto pt : new_table_protos) {
        if(not insert(internal_tables_tablename(), pt)) {
          rollback_transaction();
          return false;
        }
      }
      for(auto pf : new_field_protos) {
        if(not insert(internal_fields_tablename(), pf)) {
          rollback_transaction();
          return false;
        }
      }
    }
  } catch(...) {
    rollback_transaction();
    throw;
  }

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
  if(transaction_count_.fetch_add(1) == 0) {
    if(write_sql_to_log_) {
      LOG(INFO) << "BEGIN TRANSACTION";
    }
  }
  return true;
}

bool SQLSerializer::commit_transaction()
{
  if(transaction_count_.fetch_sub(1) == 1) {
    if(write_sql_to_log_) {
      LOG(INFO) << "COMMIT TRANSACTION";
    }
  }
  return true;
}

bool SQLSerializer::rollback_transaction()
{
  if(transaction_count_.fetch_sub(1) == 1) {
    if(write_sql_to_log_) {
      LOG(INFO) << "ROLLBACK TRANSACTION";
    }
  }
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

void SQLSerializer::set_const_data_pointers(SQLTable* t, const google::protobuf::Message* m,
  const uint64_t* parent_oid, const uint64_t* loop_id)
{
  for(auto f : t->fields)
  {
#if 0
    LOG(INFO) << f->field_name << ' ' << f->field_type << ' ' << f << ' '
              << f->field_origin << '\n';
#endif
    if(f->field_origin == f)
    {
      if(f->field_d != nullptr or f->oneof_d != nullptr)
      {
        const google::protobuf::Reflection* r = m->GetReflection();

        if(f->oneof_d or f->field_d->is_repeated() or r->HasField(*m, f->field_d)
          /* or f->field_d->message_type()==nullptr */) {
          f->set_data_const_message(m);
        } else {
          f->set_data_null();
        }
      }
      else
      {
        switch(f->field_type)
        {
        case SQLTableField::KEY_PARENT_OID:
          if(parent_oid != nullptr)
            f->set_data_const_uint64(parent_oid);
          break;
        case SQLTableField::KEY_LOOP_ID:
          if(loop_id != nullptr)
            f->set_data_const_uint64(loop_id);
          break;
        case SQLTableField::KEY_INHERITED:      // handled earlier in tree
        case SQLTableField::POD:                // handled in if clause above
          assert(0);
          break;
        }
      }
    }
  }
}

void SQLSerializer::bind_fields_from_data_pointers(const SQLTable* t, uint64_t loop_id,
  SQLStatement* stmt, bool bind_inherited_keys_only)
{
  unsigned ifield = 0;
  for(auto f : t->fields)
  {
    if(not f->db_field_present) {
      // if(f->is_key()) {
      throw std::runtime_error("SQLSerializer::bind_fields_from_data_pointers: "
        "column not present in DB : " + sub_name(t->table_name, f->field_name));
      // }
      continue;
    }

    //LOG(INFO) << t->table_name << " -> " << f->field_name;
    if(not bind_inherited_keys_only or f->is_inherited())
    {
      f = f->field_origin;

      if(f->is_data_null()) {
        stmt->bind_null(ifield);
      } else if(f->oneof_d != nullptr) {
        const google::protobuf::Message* m = f->data_const_message();
        const google::protobuf::Reflection* r = m->GetReflection();
        stmt->bind_uint64(ifield, r->HasOneof(*m, f->oneof_d) ?
          r->GetOneofFieldDescriptor(*m, f->oneof_d)->number() : 0);
      } else if(f->field_d == nullptr) {
        stmt->bind_uint64(ifield, *f->data_const_uint64());
      } else if(f->field_d->is_repeated()) {
        stmt->bind_repeated_field(ifield, loop_id, f->data_const_message(), f->field_d);
      } else {
        stmt->bind_field(ifield,  f->data_const_message(), f->field_d);
      }
      ifield++;
    }
  }
}

bool SQLSerializer::r_exec_insert(SQLTable* t, const google::protobuf::Message* m,
  uint64_t& oid, uint64_t parent_oid, uint64_t loop_id, bool ignore_errors)
{
  // This function is too long and complex. The basic idea is to
  // 1 iterate through all non-inherited fields (i.e. those whose values don't
  //   come from higher up in the tree) and place a pointer to their values in
  //   the f->data member. If the data is:
  //   - a protobuf field : f->data has a pointer to the appropriate message
  //     if the field is present in the mesage, null otherwise
  //   - a loop_id : f->data is a pointer to the uint64_t loop index
  //   - an oid : f->data is a pointer to the parent_oid for this table
  //   - a protobuf oneof index : f->data is the message, should it be present,
  //     null otherwise
  // 2 values for all fields (inherited and new) are bound to the sql statement
  //   - non repeated protobuf fields are bound using t->stmt->bind_field
  //   - repeated protobuf fields using t->stmt->bind_repeated_field
  //   - oneof by binding the probuf number of the selected field (bind_uint64)
  //     or null if no field is selected
  //   - non protobuf fields (oid, loop_id) using bind_uint64
  //   - fields without f->data are bound to NULL
  // 3 the statement is executed
  // 4 the oid for the insert is stored in "oid" to be passed back to the caller
  //   and on to any sub tables
  // 5 all sub tables are processed by recursive calls to this function
  //   - repeated sub tables are processed in a loop with loop_id passed in
  //   - simple sub tables are called directly if the appropriate field is
  //     set in the message

  set_const_data_pointers(t, m, &parent_oid, &loop_id);
  bind_fields_from_data_pointers(t, loop_id, t->stmt_insert);

  bool good = true;
  SQLStatement::StepStatus status = t->stmt_insert->step();
  if(write_sql_to_log_)LOG(INFO) << t->stmt_insert->bound_sql();
  switch(status)
  {
  case SQLStatement::ERROR:
    good = false;
    if(!ignore_errors)
    {
      LOG(ERROR) << "INSERT statement returned error: "
                 << t->stmt_insert->error_message();
      return false;
    }
    break;
  case SQLStatement::OK_NO_DATA:
    // this is what we expect
    break;
  case SQLStatement::OK_HAS_DATA:
    LOG(ERROR) << "INSERT statement returned data" << '\n'
               << t->stmt_insert->sql();
    return false;
  }
  oid = t->stmt_insert->get_oid();
  t->stmt_insert->reset();

  for(auto st : t->sub_tables)
  {
    const google::protobuf::Reflection* r = m->GetReflection();
    if(st->parent_field_d->is_repeated())
    {
      uint64_t nloop = r->FieldSize(*m, st->parent_field_d);
      for(uint64_t iloop = 0; iloop < nloop; iloop++)
      {
        const google::protobuf::Message* mi = m;
        if(st->parent_field_d->type() == FieldDescriptor::TYPE_MESSAGE) {
          mi = &r->GetRepeatedMessage(*m, st->parent_field_d, iloop); }
        uint64_t unused_sub_table_oid = 0;
        good &= r_exec_insert(st, mi, unused_sub_table_oid, oid, iloop,
                              ignore_errors);
        if(!ignore_errors and !good)return good;
      }
    }
    else
    {
      assert(st->parent_field_d->type() == FieldDescriptor::TYPE_MESSAGE);
      if(!r->HasField(*m, st->parent_field_d))goto next_sub_table;
      const auto* stm = &r->GetMessage(*m, st->parent_field_d);
      uint64_t unused_sub_table_oid = 0;
      good &= r_exec_insert(st, stm, unused_sub_table_oid, oid, 0,
                            ignore_errors);
      if(!ignore_errors and !good)return good;
    }
  next_sub_table:
    ;
  }

  return good;
}

namespace {
  bool field_selected(const google::protobuf::FieldDescriptor* f,
    const std::map<const google::protobuf::OneofDescriptor*,int>& oneof_map)
  {
    const google::protobuf::OneofDescriptor* c_oneof = f->containing_oneof();
    if(c_oneof != nullptr and
        (oneof_map.count(c_oneof)==0 or oneof_map.at(c_oneof) != f->number()))
      return false;
    return true;
  }
} // anonymous namespace

bool SQLSerializer::
exec_select_by_oid(SQLTable* t, uint64_t oid, google::protobuf::Message* m_root,
  bool ignore_errors)
{
  t->stmt_select_oid->bind_uint64(0, oid);

  google::protobuf::Message* m_parent = nullptr;
  std::vector<uint64_t> m_parent_loop_index;
  for(auto ifield : t->fields) {
    if(ifield->is_inherited() and not ifield->is_root_oid) {
      m_parent_loop_index.emplace_back(0);
    }
  }
  std::vector<uint64_t> step_loop_index = m_parent_loop_index;

  unsigned row_count = 0;
  SQLStatement::StepStatus status = t->stmt_select_oid->step();
  if(write_sql_to_log_)LOG(INFO) << t->stmt_select_oid->bound_sql();
  while(status == SQLStatement::OK_HAS_DATA)
  {
    bool good = true;
    auto ifield = t->fields.begin();
    unsigned icol = 0;
    uint64_t loop_id = 0;
    google::protobuf::Message* m = nullptr;
    std::map<const OneofDescriptor*, int> oneof_map;

    ++row_count;

    if((*ifield)->is_root_oid)++ifield;

    while(icol<step_loop_index.size()) {
      step_loop_index[icol] = t->stmt_select_oid->extract_uint64(icol, &good);
      if(!good)goto next_step;
      ++icol;
      ++ifield;
    }

    if(m_parent == nullptr or step_loop_index != m_parent_loop_index) {
      // Walk m_parent to correct place
      m_parent = m_root;
      m_parent_loop_index = step_loop_index;
      unsigned iloop_index = 0;
      for(auto* d : t->root_field_d_path) {
        auto* r = m_parent->GetReflection();
        if(d->is_repeated()) {
          assert(iloop_index < m_parent_loop_index.size());
          int loop_index = m_parent_loop_index[iloop_index];
          ++iloop_index;
          assert(loop_index < r->FieldSize(*m, d));
          // if(r->FieldSize(*m, d) <= loop_index)r->AddMessage(m_parent, d);
          m_parent = r->MutableRepeatedMessage(m_parent, d, loop_index);
        } else {
          m_parent = r->MutableMessage(m_parent, d);
        }
      }
      assert(iloop_index == m_parent_loop_index.size());
    }

    if(t->parent_table == nullptr) {
      m = m_root;
    } else if (not t->parent_field_d->is_repeated()) {
      assert(t->parent_field_d->type()==FieldDescriptor::TYPE_MESSAGE);
      m = m_parent->GetReflection()->MutableMessage(m_parent, t->parent_field_d);
    }

    if(ifield == t->fields.end()) goto next_step;

    if((*ifield)->field_type == SQLTableField::KEY_LOOP_ID) {
      loop_id = t->stmt_select_oid->extract_uint64(icol, &good);
      if(!good)goto next_step;
      ++icol;
      ++ifield;
    }

    if(m == nullptr) {
      // In this case we are talkling about a repeated field
      auto* r = m_parent->GetReflection();
      if(t->parent_field_d->type()==FieldDescriptor::TYPE_MESSAGE) {
        if(r->FieldSize(*m_parent, t->parent_field_d) <= loop_id)
          r->AddMessage(m_parent, t->parent_field_d);
        m = r->MutableRepeatedMessage(m_parent, t->parent_field_d, loop_id);
      } else {
        m = m_parent;
      }
    }

    while(ifield != t->fields.end() and (*ifield)->oneof_d) {
      if(not t->stmt_select_oid->column_is_null(icol)) {
        oneof_map[(*ifield)->oneof_d] = t->stmt_select_oid->extract_uint64(icol, &good);
      }
      if(!good)goto next_step;
      ++icol;
      ++ifield;
    }

    while(ifield != t->fields.end()) {
      if(not t->stmt_select_oid->column_is_null(icol) and
        field_selected((*ifield)->field_d, oneof_map))
      {
        if((*ifield)->field_d->is_repeated())
          good = t->stmt_select_oid->extract_repeated_field(icol, loop_id, m, (*ifield)->field_d);
        else
          good = t->stmt_select_oid->extract_field(icol, m, (*ifield)->field_d);
        if(!good)goto next_step;
      }
      ++icol;
      ++ifield;
    }

next_step:
    if(not good) {
      throw std::runtime_error(
        std::string("SQLSerializer::exec_select_by_oid: could not extract column ")
        + std::to_string(icol) + " [" + (*ifield)->field_name + "]");
    }

    status = t->stmt_select_oid->step();
  }

  if(write_sql_to_log_) {
    if(row_count!=1)LOG(INFO) << "Retrieved " << row_count << " rows";
    else LOG(INFO) << "Retreived " << row_count << " row";
  }

  if(status == SQLStatement::ERROR and not ignore_errors) {
    LOG(ERROR) << "SQL statement returned error: "
      << t->stmt_select_oid->error_message() << '\n'
      << "SQL: " << t->stmt_select_oid->sql();
  }

  t->stmt_select_oid->reset();
  return status != SQLStatement::ERROR;
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

std::string SQLSerializer::sql_select_field_spec(const SQLTableField* f)
{
  return sql_field_name(f->field_name);
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

std::string SQLSerializer::sql_select_where_oid_equals(const SQLTable* t)
{
  if(not t->db_table_present) {
    return {};
  }

  std::string oid_name;
  if(t == t->root_table) {
    oid_name = sql_oid_column_name();
  } else {
    oid_name = sub_name(t->root_table->table_name, sql_oid_column_name());
  }

  std::ostringstream sql;
  sql << "SELECT\n";
  bool has_field = false;
  for ( auto f : t->fields ) {
    if(f->db_field_present and f->field_name!=oid_name) {
      if(has_field) {
        sql << ",\n";
      }
      sql << "  " << sql_field_name(f->field_name);
      has_field = true;
    }
  }
  if(has_field) {
    sql << '\n';
  } else {
    sql << "NULL\n";
  }
  sql << "FROM " << sql_table_name(t->table_name) << " WHERE "
    << sql_field_name(oid_name) << " == ?";
  return sql.str();
}

std::string SQLSerializer::sql_count_entries(const std::string& table_name)
{
  std::ostringstream sql;
  sql << "SELECT COUNT(*) FROM " << sql_table_name(table_name);
  return sql.str();
}

std::string SQLSerializer::sql_select_oids(const std::string& table_name)
{
  std::ostringstream sql;
  sql << "SELECT " << sql_field_name(sql_oid_column_name()) << " FROM "
    << sql_table_name(table_name);
  return sql.str();
}

std::string SQLSerializer::internal_tables_tablename()
{
  return "calin.tables";
}

std::string SQLSerializer::internal_fields_tablename()
{
  return "calin.fields";
}

std::string SQLSerializer::internal_params_tablename()
{
  return "calin.params";
}

void SQLSerializer::retrieve_db_tables_and_fields()
{
  const google::protobuf::Descriptor* d;
  SQLTable* t;

  d = calin::ix::io::sql_serializer::SQLTable::descriptor();
  t = schema_[internal_tables_tablename()][d];
  if(t == nullptr) {
    t = make_sqltable_tree(internal_tables_tablename(), d, "calin database table list");
    force_sqltable_tree_db_presence(t);
    schema_[internal_tables_tablename()][d] = t;
  }

  d = calin::ix::io::sql_serializer::SQLTableField::descriptor();
  t = schema_[internal_fields_tablename()][d];
  if(t == nullptr) {
    t = make_sqltable_tree(internal_fields_tablename(), d, "calin database field list");
    force_sqltable_tree_db_presence(t);
    schema_[internal_fields_tablename()][d] = t;
  }

  std::vector<uint64_t> oids;

  begin_transaction();
  oids = retrieve_all_oids(internal_tables_tablename());
  for(auto oid : oids) {
    auto* proto_table = new calin::ix::io::sql_serializer::SQLTable;
    if(retrieve_by_oid(internal_tables_tablename(), oid, proto_table)) {
      db_tables_[proto_table->table_name()] = proto_table;
    }
  }

  oids = retrieve_all_oids(internal_fields_tablename());
  for(auto oid : oids) {
    auto* proto_field = new calin::ix::io::sql_serializer::SQLTableField;
    if(retrieve_by_oid(internal_fields_tablename(), oid, proto_field)) {
      db_table_fields_[sub_name(proto_field->table_name(), proto_field->field_name())] = proto_field;
    }
  }
  commit_transaction();
}

void SQLSerializer::create_db_tables_and_fields()
{
  std::unique_ptr<SQLTable> tp { make_sqltable_tree(internal_params_tablename(),
    calin::ix::io::sql_serializer::SQLSerializerParams::descriptor(), "calin sql serializer parameters") };
  std::unique_ptr<SQLTable> tt { make_sqltable_tree(internal_tables_tablename(),
    calin::ix::io::sql_serializer::SQLTable::descriptor(), "calin database table list") };
  std::unique_ptr<SQLTable> tf { make_sqltable_tree(internal_fields_tablename(),
    calin::ix::io::sql_serializer::SQLTableField::descriptor(), "calin database field list") };

  begin_transaction();

  if(not do_create_or_extend_tables(tp.get(), false)) {
    rollback_transaction();
    throw std::runtime_error("create_db_tables_and_fields: could not create table: " +
      internal_params_tablename());
  }

  if(not do_create_or_extend_tables(tt.get(), false)) {
    rollback_transaction();
    throw std::runtime_error("create_db_tables_and_fields: could not create table: " +
      internal_tables_tablename());
  }

  if(not do_create_or_extend_tables(tf.get(), false)) {
    rollback_transaction();
    throw std::runtime_error("create_db_tables_and_fields: could not create table: " +
      internal_fields_tablename());
  }

  calin::ix::io::sql_serializer::SQLSerializerParams params;
  params.set_version(native_db_version_);
  if(not insert(internal_params_tablename(), &params)) {
    rollback_transaction();
    throw std::runtime_error("create_db_tables_and_fields: could not insert into table: " +
      internal_params_tablename());
  }

  for(auto pt : db_tables_) {
    if(not insert(internal_tables_tablename(), pt.second)) {
      rollback_transaction();
      throw std::runtime_error("create_db_tables_and_fields: could not insert into table: " +
        internal_tables_tablename());
    }
  }

  for(auto pf : db_table_fields_) {
    if(not insert(internal_fields_tablename(), pf.second)) {
      rollback_transaction();
      throw std::runtime_error("create_db_tables_and_fields: could not insert into table: " +
        internal_tables_tablename());
    }
  }
  commit_transaction();
}
