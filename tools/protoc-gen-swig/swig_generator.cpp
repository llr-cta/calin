/* 

   calin/tools/swig_generator.cpp -- Stephen Fegan -- 2015-12-03

   Code generator for SWIG interface files

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include "swig_generator.hpp"

#include <google/protobuf/descriptor.h>
#include <google/protobuf/io/printer.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

using std::string;
using namespace calin::tools::protoc_gen_swig;
using google::protobuf::io::Printer;

SwigGenerator::~SwigGenerator()
{

}

namespace {

string pb_to_gen_filename(string pb_filename, string extension = ".pb.i")
{
  return pb_filename.substr(0,pb_filename.find_last_of('.')) + extension;
}

void print_includes(Printer* I, const google::protobuf::FileDescriptor * file,
                    string directive = "#include", string extension = ".pb.h")
{
  for(int i=0; i<file->dependency_count(); i++)
    I->Print("$directive$<$file$>\n","directive", directive, "file",
             pb_to_gen_filename(file->dependency(i)->name(),extension));
  for(int i=0; i<file->public_dependency_count(); i++)
    I->Print("$directive$<$file$>\n","directive", directive, "file",
             pb_to_gen_filename(file->public_dependency(i)->name(),extension));
  for(int i=0; i<file->weak_dependency_count(); i++) 
    I->Print("$directive$<$file$>\n","directive", directive, "file",
             pb_to_gen_filename(file->weak_dependency(i)->name(),extension));
}

void print_message(Printer* I, const google::protobuf::Descriptor* d)
{
  I->Print("class $name$\n"
           "{\n"
           " public:\n",
           "name", d->name());
  I->Indent();

  for(unsigned i=0;i<d->field_count();i++)
  {
    auto* f = d->field(i);
    if(i>0)I->Print("\n");
    std::map<string, string> vars;
    vars["name"]   = f->name();
    vars["type"]   = f->cpp_type_name();
    vars["index"]  = "";
    vars["get"]    = "";
    vars["set"]    = "set_";

    if(f->is_repeated())
    {
      vars["index"]  = "int index";
      vars["set"]    = "add_";
    }
    
    if(f->is_map())
    {
      
    }
    else if(f->type() == google::protobuf::FieldDescriptor::TYPE_MESSAGE)
    {
      vars["get"]    = "mutable_";
      vars["set"]    = "";
      vars["type"]   = f->message_type()->full_name();
    }
    else if(f->cpp_type_name() == std::string("string"))
    {
      vars["type"] = "const std::string&";
    }

    if(f->containing_oneof())
      I->Print(vars, "bool has_$name$() const;\n");

    if(f->is_repeated())
      I->Print(vars, "int $name$_size() const;\n");

    I->Print(vars, "$type$ $get$$name$($index$) const;\n");
    if(!vars["set"].empty())
      I->Print(vars, "void $set$$name$($type$ x);\n");
    I->Print(vars, "void clear_$name$();\n");
  }
  
  I->Outdent();
  I->Print("};\n");
}

} // anonymous namespace

bool SwigGenerator::
Generate(const google::protobuf::FileDescriptor * file,
         const string & parameter,
         google::protobuf::compiler::GeneratorContext *context,
         string * error) const
{
  auto I_stream = context->Open(pb_to_gen_filename(file->name(),".pb.i"));
  Printer* I = new Printer(I_stream,'$');

  I->Print("// Auto-generated from \"$file$\". "
           "Do not edit by hand.\n\n","file",file->name());
  
  std::map<string,string> vars;
  if(file->package().find('.') != string::npos)
  {
    auto index = file->package().find_last_of('.');
    I->Print("%module (package=\"$package$\") $module$\n",
             "package", file->package().substr(0,index),
             "module", file->package().substr(index+1));
  }
  else
  {
     I->Print("%module $module$\n",
              "module", file->package());
  }

  I->Print("\n%{\n");
  I->Indent();
  print_includes(I, file);
  I->Print("#define SWIG_FILE_WITH_INIT\n");
  I->Outdent();
  I->Print("%}\n\n");

  I->Print("%init %{\n"
           "  import_array();\n"
           "%}\n"
           "\n"
           "%include \"package_wide_definitions.i\"\n\n");

  print_includes(I, file, "%import", ".pb.i");

  for(int i=0;i<file->message_type_count();i++)
  {
    I->Print("\n");
    print_message(I, file->message_type(i));
  }
    
  delete I;
  delete I_stream;
  return true;
}
