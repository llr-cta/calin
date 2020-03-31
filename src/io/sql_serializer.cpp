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

SQLSerializer::~SQLSerializer()
{
  // nothing to see here
}

calin::io::sql_serializer::SQLTable* SQLSerializer::
make_sqltable_tree(const std::string& table_name,
                   const google::protobuf::Descriptor* d)
{
  return r_make_sqltable_tree(table_name, d,
    /* parent_table= */ nullptr, /* parent_field_d= */ nullptr);
}

calin::io::sql_serializer::SQLTable* SQLSerializer::
r_make_sqltable_tree(const std::string& table_name,
                     const google::protobuf::Descriptor* d,
                     SQLTable* parent_table,
                     const google::protobuf::FieldDescriptor* parent_field_d)
{
  SQLTable* t { new SQLTable };
  t->parent_table        = parent_table;
  t->table_name          = table_name;
  t->parent_field_d      = parent_field_d;

  for(int ioneof = 0; ioneof<d->oneof_decl_count(); ioneof++)
  {
    const OneofDescriptor* oo { d->oneof_decl(ioneof) };
    SQLTableField* tf { new SQLTableField };
    tf->table             = t;
    tf->field_origin      = tf;
    tf->field_type        = SQLTableField::POD;
    tf->field_name        = oo->name();
    tf->oneof_d           = oo;

    t->fields.push_back(tf);
  }

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
                                       f->message_type(), t, f);

      if(parent_field_d and parent_field_d->is_map())
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
