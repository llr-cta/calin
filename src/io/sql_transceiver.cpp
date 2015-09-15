/* 

   calin/io/sql_transceiver.cpp -- Stephen Fegan -- 2015-09-08

   Provides reading and writing of protobuf structures to SQLite3 databases

*/

#include <deque>

#include "proto/calin.pb.h"
#include "io/sql_transceiver.hpp"

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
  SQLTable* t { make_sqltable_tree(table_name, d_data) };
  propagate_keys(t);
  iterate_over_tables(t, [this](const SQLTable* it){
      std::cout << sql_create_table(it) << '\n'; });
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
make_sqltable_tree(const std::string& table_name,
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
      sub_table = make_sqltable_tree(sub_name(table_name, f->name()),
                                     f->message_type(), t, f);

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
      tf->field_type        = SQLTableField::POD;
      tf->field_name        = f->name();
      tf->field_d           = f;

      sub_table->fields.push_back(tf);
    }
    else
    {
      SQLTableField* tf { new SQLTableField };
      tf->table             = t;
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
  iterate_over_tables(t, [](SQLTable* t) { if(t->fields.empty()) {
        for(auto st : t->sub_tables) { st->parent_table = t->parent_table; }
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

std::string SQLTransceiver::sql_type(const google::protobuf::FieldDescriptor* d)
{
  if(!d)return "int32";
  if(d->cpp_type_name())return d->cpp_type_name();
  return "unknown";
}

std::string SQLTransceiver::sql_create_table(const SQLTable* t) 
{
  std::ostringstream sql;
  std::vector<const SQLTableField*> keys;
  sql << "CREATE TABLE " << sql_table_name(t->table_name) << " ( \n";
  for ( auto f : t->fields )
  {
    if(f->is_key())keys.push_back(f);
    sql << "  " << sql_field_name(f->field_name) << ' ' << sql_type(f->field_d);
    if(f != t->fields.back())sql << ',';
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
  sql << ')';
  if(!keys.empty())
  {
    sql << " PRIMARY KEY (\n";
    for(auto f : keys)
    {
      sql << "  " << sql_field_name(f->field_name);
      if(f != keys.back())sql << ',';
      sql << '\n';
    }
    sql << ')';
  }
  return sql.str();
}
