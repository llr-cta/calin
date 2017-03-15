/*

  calin/tools/protoc-gen-calin-cpp-extensions/calin_cpp_extensions_generator.cpp -- Stephen Fegan -- 2017-04-14

  Code generator for CPP extensions

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <iostream>
#include <sstream>
#include <numeric>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/descriptor.pb.h>
#include <google/protobuf/io/printer.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include <util/string.hpp>
#include <calin.pb.h>
#include "calin_cpp_extensions_generator.hpp"

using std::string;
using namespace calin::util::string;
using namespace calin::tools::calin_cpp_extensions_generator;
using google::protobuf::io::Printer;

CalinCppExtensionsGenerator::~CalinCppExtensionsGenerator()
{
  // nothing to see here
}

namespace {

string pb_to_cpp_filename(string pb_filename, string extension = ".pb.cc")
{
  return pb_filename.substr(0,pb_filename.find_last_of('.')) + extension;
}

string pb_to_hpp_filename(string pb_filename, string extension = ".pb.h")
{
  return pb_to_cpp_filename(pb_filename, extension);
}

bool is_numeric_type(google::protobuf::FieldDescriptor::Type type)
{
  switch(type)
  {
  case google::protobuf::FieldDescriptor::TYPE_DOUBLE:    // fall through
  case google::protobuf::FieldDescriptor::TYPE_FLOAT:     // fall through
  case google::protobuf::FieldDescriptor::TYPE_INT64:     // fall through
  case google::protobuf::FieldDescriptor::TYPE_UINT64:    // fall through
  case google::protobuf::FieldDescriptor::TYPE_INT32:     // fall through
  case google::protobuf::FieldDescriptor::TYPE_FIXED64:   // fall through
  case google::protobuf::FieldDescriptor::TYPE_FIXED32:   // fall through
  case google::protobuf::FieldDescriptor::TYPE_UINT32:    // fall through
  case google::protobuf::FieldDescriptor::TYPE_SFIXED32:  // fall through
  case google::protobuf::FieldDescriptor::TYPE_SFIXED64:  // fall through
  case google::protobuf::FieldDescriptor::TYPE_SINT32:    // fall through
  case google::protobuf::FieldDescriptor::TYPE_SINT64:    // fall through
    return true;
  default:
    return false;
  }
  assert(0);
  return false;
}

bool is_string_type(google::protobuf::FieldDescriptor::Type type)
{
  switch(type)
  {
  case google::protobuf::FieldDescriptor::TYPE_STRING:    // fall through
  case google::protobuf::FieldDescriptor::TYPE_BYTES:     // fall through
    return true;
  default:
    return false;
  }
  assert(0);
  return false;
}

bool is_non_numeric_scalar_type(google::protobuf::FieldDescriptor::Type type)
{
  switch(type)
  {
  case google::protobuf::FieldDescriptor::TYPE_BOOL:      // fall through
  case google::protobuf::FieldDescriptor::TYPE_ENUM:      // fall through
    return true;
  default:
    return false;
  }
  assert(0);
  return false;
}

bool integrate_field(const google::protobuf::FieldDescriptor* f,
  google::protobuf::io::Printer& printer, string * error)
{
  const google::protobuf::FieldOptions* fopt = &f->options();
  const calin::FieldOptions* cfo = &fopt->GetExtension(calin::CFO);
  if(is_numeric_type(f->type())) {
    if(f->is_repeated()) {
      switch(cfo->integration_algorithm()) {
      case calin::FieldOptions::IGNORE:
        break;
      case calin::FieldOptions::DEFAULT:
      case calin::FieldOptions::SUM:
        printer.Print(
          "for(int i=0; i<from.$name$_size(); i++) {\n"
          "  if(i<this->$name$_size())this->set_$name$(i, this->$name$(i)+from.$name$(i));\n"
          "  else this->add_$name$(from.$name$(i));\n"
          "}\n", "name", f->name());
        break;
      case calin::FieldOptions::REPLACE:
        printer.Print("this->clear_$name();\n", "name", f->name());
        // Fall through
      case calin::FieldOptions::APPEND:
        printer.Print(
          "for(int i = 0;i < from.$name$_size(); i++)\n"
          "  this->add_$name$(from.$name$(i));\n", "name", f->name());
        break;
      default:
        assert(false);
        break;
      }
    } else { // if(f->is_repeated()) {
      switch(cfo->integration_algorithm()) {
      case calin::FieldOptions::IGNORE:
        break;
      case calin::FieldOptions::DEFAULT:
      case calin::FieldOptions::SUM:
        printer.Print(
          "this->set_$name$(this->$name$()+from.$name$());\n", "name", f->name());
        break;
      case calin::FieldOptions::REPLACE:
        printer.Print(
          "this->set_$name$(from.$name$());\n", "name", f->name());
          break;
      case calin::FieldOptions::APPEND:
        *error = "Cannot select APPEND algorithm for field: "+f->name();
        return false;
      default:
        assert(false);
        break;
      }
    }
  } else if(is_string_type(f->type())) {
    if(f->is_repeated()) {
      switch(cfo->integration_algorithm()) {
      case calin::FieldOptions::IGNORE:
        break;
      case calin::FieldOptions::DEFAULT:
      case calin::FieldOptions::REPLACE:
        printer.Print("this->clear_$name();\n", "name", f->name());
        // Fall through
      case calin::FieldOptions::APPEND:
        printer.Print(
          "for(int i = 0;i < from.$name$_size(); i++)\n"
          "  this->add_$name$(from.$name$(i));\n", "name", f->name());
        break;
      case calin::FieldOptions::SUM:
        *error = "Cannot select SUM algorithm for field: "+f->name();
        return false;
      default:
        assert(false);
        break;
      }
    } else { // if(f->is_repeated()) {
      switch(cfo->integration_algorithm()) {
      case calin::FieldOptions::IGNORE:
        break;
      case calin::FieldOptions::DEFAULT:
      case calin::FieldOptions::REPLACE:
        printer.Print(
          "this->set_$name$(from.$name$());\n", "name", f->name());
          break;
      case calin::FieldOptions::APPEND:
        printer.Print(
          "this->set_$name$(this->$name$()+from.$name$());\n", "name", f->name());
        break;
      case calin::FieldOptions::SUM:
        *error = "Cannot select SUM algorithm for field: "+f->name();
        return false;
      default:
        assert(false);
        break;
      }
    }
  } else if(is_non_numeric_scalar_type(f->type())) {
    if(f->is_repeated()) {
      switch(cfo->integration_algorithm()) {
      case calin::FieldOptions::IGNORE:
        break;
      case calin::FieldOptions::DEFAULT:
      case calin::FieldOptions::REPLACE:
        printer.Print("this->clear_$name();\n", "name", f->name());
        // Fall through
      case calin::FieldOptions::APPEND:
        printer.Print(
          "for(int i = 0;i < from.$name$_size(); i++)\n"
          "  this->add_$name$(from.$name$(i));\n", "name", f->name());
        break;
      case calin::FieldOptions::SUM:
        *error = "Cannot select SUM algorithm for field: "+f->name();
        return false;
      default:
        assert(false);
        break;
      }
    } else { // if(f->is_repeated()) {
      switch(cfo->integration_algorithm()) {
      case calin::FieldOptions::IGNORE:
        break;
      case calin::FieldOptions::DEFAULT:
      case calin::FieldOptions::REPLACE:
        printer.Print(
          "this->set_$name$(from.$name$());\n", "name", f->name());
          break;
      case calin::FieldOptions::APPEND:
        *error = "Cannot select APPEND algorithm for field: "+f->name();
        return false;
      case calin::FieldOptions::SUM:
        *error = "Cannot select SUM algorithm for field: "+f->name();
        return false;
      default:
        assert(false);
        break;
      }
    }
  }
  return true;
}

}

bool CalinCppExtensionsGenerator::
Generate(const google::protobuf::FileDescriptor * file,
         const string & parameter,
         google::protobuf::compiler::GeneratorContext *context,
         string * error) const
{
  // Print IntegrateFrom functions for all messages that request them
  std::vector<const google::protobuf::Descriptor*> automatic_mif;
  for(int i=0;i<file->message_type_count();i++)
  {
    const google::protobuf::Descriptor* d = file->message_type(i);
    const google::protobuf::MessageOptions* mopt = &d->options();
    if(!mopt->HasExtension(calin::CMO))continue;
    const calin::MessageOptions* cmo = &mopt->GetExtension(calin::CMO);
    if(cmo->message_integration_function() != calin::MessageOptions::MIF_NONE)
    {
      google::protobuf::scoped_ptr<google::protobuf::io::ZeroCopyOutputStream> output(
          context->OpenForInsert(pb_to_hpp_filename(file->name()),
          std::string("class_scope:")+d->full_name()));
      google::protobuf::io::Printer printer(output.get(), '$');
      printer.Print("void IntegrateFrom(const $type$& from);\n", "type", d->name());
    }

    if(cmo->message_integration_function() == calin::MessageOptions::MIF_AUTOMATIC)
      automatic_mif.push_back(d);
  }

  if(not automatic_mif.empty())
  {
    google::protobuf::scoped_ptr<google::protobuf::io::ZeroCopyOutputStream> output(
        context->OpenForInsert(pb_to_cpp_filename(file->name()),
        std::string("namespace_scope")));
    google::protobuf::io::Printer printer(output.get(), '$');
    for(const google::protobuf::Descriptor* d : automatic_mif)
    {
      printer.Print("void $type$::IntegrateFrom(const $type$& from)\n"
                    "{\n", "type", d->name());
      printer.Indent();
      // Fields
      for(int i=0;i<d->field_count();i++)
      {
        const google::protobuf::FieldDescriptor* f = d->field(i);
        if(not integrate_field(f, printer, error))return false;
      }
      printer.Outdent();
      printer.Print("}\n\n");
    }
  }

  return true;
}
