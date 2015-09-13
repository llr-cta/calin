/* 

   calin/io/sql_transceiver.cpp -- Stephen Fegan -- 2015-09-08

   Provides reading and writing of protobuf structures to SQLite3 databases

*/

#include "proto/calin.pb.h"
#include "io/sql_transceiver.hpp"

using namespace calin::io::sql_transceiver;
using namespace google::protobuf;

std::vector<std::pair<std::string,std::string>>
SQLTransceiver::list_all_table_columns(const std::string& table_name,
                                       const google::protobuf::Descriptor* d)
{ 
  std::vector<std::pair<std::string,std::string>> fields;
  SQLTable* t = make_table(table_name, d);
  fields = extract_all_table_columns(t);
  return fields;
}

std::vector<std::pair<std::string,std::string> >
SQLTransceiver::extract_all_table_columns(SQLTable* t)
{
  std::vector<std::pair<std::string,std::string>> fields;
  for(auto f : t->fields)
    fields.push_back(std::make_pair(t->table_name, f->field_name));
  for(auto st : t->sub_tables)
  {
    std::vector<std::pair<std::string,std::string>> sub_fields =
        extract_all_table_columns(st);
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
make_table(const std::string& table_name,
           const google::protobuf::Descriptor* d,
           SQLTable* parent_table,
           const google::protobuf::FieldDescriptor* parent_field_d)
{
  SQLTable* t { new SQLTable };
  t->parent_table   = parent_table;
  t->table_name     = table_name;
  t->parent_field_d = parent_field_d;

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
      sub_table = make_table(sub_name(table_name, f->name()),
                             f->message_type(), t, f);

      if((parent_field_d and parent_field_d->is_map()) or
         (!f->is_repeated() and fopt->HasExtension(CFO) and
          fopt->GetExtension(CFO).sql().inline_message()))
      {
        // Message can be inlined - so move its fields into primary table

        for(auto ifield : sub_table->fields)
        {
          ifield->table = t;
          //if(!parent_field_d or !parent_field_d->is_map())
          ifield->field_name = sub_name(f->name(), ifield->field_name);
          ifield->table = t;
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
      tf->table        = sub_table;
      tf->field_type   = SQLTableField::POD;
      tf->field_name   = f->name();
      tf->field_d      = f;

      sub_table->fields.push_back(tf);
    }
    else
    {
      SQLTableField* tf { new SQLTableField };
      tf->table        = t;
      tf->field_type   = SQLTableField::POD;
      tf->field_name   = f->name();
      tf->field_d      = f;
      
      if(parent_table == nullptr and
         fopt->HasExtension(CFO) and fopt->GetExtension(CFO).sql().is_key())
        tf->field_type   = SQLTableField::KEY_PROTO_DEFINED;

      t->fields.push_back(tf);
      continue;
    }

    assert(sub_table);
    
    if(f->is_map())
    {
      sub_table->fields.front()->field_type   = SQLTableField::KEY_MAP_KEY;
    }
    else if(f->is_repeated())
    {
      SQLTableField* tf { new SQLTableField };
      tf->table        = sub_table;
      tf->field_type   = SQLTableField::KEY_LOOP_ID;
      tf->field_name   = sub_name(sub_table->table_name, "loop_id");
      tf->field_d      = nullptr;
      sub_table->fields.insert(sub_table->fields.begin(), tf);
    }

    t->sub_tables.insert(t->sub_tables.end(), sub_table);
  }
  return t;
}

std::string SQLTransceiver::
key_name(const std::string& name)
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
