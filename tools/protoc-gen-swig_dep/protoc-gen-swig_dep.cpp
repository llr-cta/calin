/*

   calin/tools/protoc-gen-swig_dep.cpp -- Stephen Fegan -- 2016-06-28

   Procobuf compiler plugin for generating SWIG dependencies interface files

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

#include <iostream>
#include <string>

#include <google/protobuf/compiler/code_generator.h>
#include <google/protobuf/compiler/plugin.h>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/descriptor.pb.h>
#include <google/protobuf/io/printer.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

using namespace std;
using google::protobuf::io::Printer;

namespace {

string pb_to_gen_filename(string pb_filename, string extension = ".pb.i")
{
  return pb_filename.substr(0,pb_filename.find_last_of('.')) + extension;
}

void write_dep(Printer* I, const string& dep, string& sep)
{
  if(dep == "calin.proto")return;
  if(dep == "google/protobuf/descriptor.proto")return;

  string cmake_sep = ";";
  I->Print("$sep$$file$", "sep", sep, "file", pb_to_gen_filename(dep));
  sep = cmake_sep;
}


} // anonymous namespace

class SwigDepGenerator: public google::protobuf::compiler::CodeGenerator
{
 public:
  virtual ~SwigDepGenerator();
  virtual bool Generate(const google::protobuf::FileDescriptor * file,
                        const std::string & parameter,
                        google::protobuf::compiler::GeneratorContext *context,
                        std::string * error) const;
};

SwigDepGenerator::~SwigDepGenerator()
{
  // nothing to see here
}

bool SwigDepGenerator::Generate(const google::protobuf::FileDescriptor * file,
  const std::string & parameter,
  google::protobuf::compiler::GeneratorContext *context,
  std::string * error) const
{
  auto I_stream = context->Open(pb_to_gen_filename(file->name(),".pb.swig_dep"));
  Printer* I = new Printer(I_stream,'$');

  string sep = "";

  for(int i=0; i<file->dependency_count(); i++)
    write_dep(I, file->dependency(i)->name(), sep);
  for(int i=0; i<file->public_dependency_count(); i++)
    write_dep(I, file->public_dependency(i)->name(), sep);
  for(int i=0; i<file->weak_dependency_count(); i++)
    write_dep(I, file->weak_dependency(i)->name(), sep);

  if(sep != "")I->Print("\n");

  delete I;
  delete I_stream;

  return true;
}

int main(int argc, char** argv) {
  SwigDepGenerator generator;
  return google::protobuf::compiler::PluginMain(argc, argv, &generator);
}
