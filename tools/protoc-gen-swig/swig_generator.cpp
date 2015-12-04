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

string pb_to_i_filename(string pb_filename)
{
  return pb_filename.substr(0,pb_filename.find_last_of('.')) + ".pb.i";
}

string pb_to_h_filename(string pb_filename)
{
  return pb_filename.substr(0,pb_filename.find_last_of('.')) + ".pb.h";
}

} // anonymous namespace

bool SwigGenerator::
Generate(const google::protobuf::FileDescriptor * file,
         const string & parameter,
         google::protobuf::compiler::GeneratorContext *context,
         string * error) const
{
  auto I_stream = context->Open(pb_to_i_filename(file->name()));
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
  for(int i=0; i<file->dependency_count(); i++)
    I->Print("#include<$file$>\n","file",
             pb_to_h_filename(file->dependency(i)->name()));
  for(int i=0; i<file->public_dependency_count(); i++)
    I->Print("#include<$file$>\n","file",
             pb_to_h_filename(file->public_dependency(i)->name()));
  for(int i=0; i<file->weak_dependency_count(); i++)
    I->Print("#include<$file$>\n","file",
             pb_to_h_filename(file->weak_dependency(i)->name()));
  I->Print("#define SWIG_FILE_WITH_INIT\n");
  I->Outdent();
  I->Print("%}\n\n");

  I->Print("%init %{\n"
           "  import_array();\n"
           "%}\n"
           "\n"
           "%include \"package_wide_definitions.i\"\n\n");

  for(int i=0; i<file->dependency_count(); i++)
    I->Print("#import \"$file$\"\n","file",
             pb_to_i_filename(file->dependency(i)->name()));
  for(int i=0; i<file->public_dependency_count(); i++)
    I->Print("#import \"$file$\"\n","file",
             pb_to_i_filename(file->public_dependency(i)->name()));
  for(int i=0; i<file->weak_dependency_count(); i++)
    I->Print("#import \"$file$\"\n","file",
             pb_to_i_filename(file->weak_dependency(i)->name()));
  

  delete I;
  delete I_stream;
  return true;
}
