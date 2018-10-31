/*

   calin/tools/protoc-gen-swig_dep.cpp -- Stephen Fegan -- 2016-06-28

   Procobuf compiler plugin for generating SWIG dependencies interface files

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

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
#include <unistd.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <pwd.h>
#include <cassert>
#include <cstdlib>

#include <google/protobuf/compiler/code_generator.h>
#include <google/protobuf/compiler/plugin.h>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/descriptor.pb.h>
#include <google/protobuf/io/printer.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

using namespace std;
using google::protobuf::io::Printer;

namespace {

string pb_basename(string pb_filename, string extension = ".pb.i")
{
  string fn = pb_filename.substr(0,pb_filename.find_last_of('.'));
  if(fn.find_last_of('/') != string::npos) {
    fn = fn.substr(fn.find_last_of('/')+1);
  }
  return fn;
}

string pb_to_gen_filename(string pb_filename, string extension = ".pb.i")
{
  return pb_filename.substr(0,pb_filename.find_last_of('.')) + extension;
}

void write_dep(Printer* I, const string& dep, string& sep,
  const std::string & src_dir, const std::string & bin_dir)
{
  string cmake_sep = ";";
  string swig_file = pb_to_gen_filename(dep);

  string test_file = src_dir + "/swig/" + swig_file;
  if(access(test_file.c_str(),R_OK)==0) {
    I->Print("$sep$$file$", "sep", sep, "file", test_file);
    sep = cmake_sep;
    return;
  }

  test_file = src_dir + "/proto/" + dep;
  if(access(test_file.c_str(),R_OK)==0) {
    test_file = "protoc-gen-swig-" + pb_basename(dep);
    // test_file = bin_dir + "/proto/" + swig_file;
    I->Print("$sep$$file$", "sep", sep, "file", test_file);
    sep = cmake_sep;
    return;
  }
}


} // anonymous namespace

class SwigDepGenerator: public google::protobuf::compiler::CodeGenerator
{
 public:
  virtual ~SwigDepGenerator();
  virtual bool Generate(const google::protobuf::FileDescriptor * file,
                        const std::string & base_dir,
                        google::protobuf::compiler::GeneratorContext *context,
                        std::string * error) const;
};

SwigDepGenerator::~SwigDepGenerator()
{
  // nothing to see here
}

bool SwigDepGenerator::Generate(const google::protobuf::FileDescriptor * file,
  const std::string & parameters,
  google::protobuf::compiler::GeneratorContext *context,
  std::string * error) const
{
  auto I_stream = context->Open(pb_to_gen_filename(file->name(),".pb.swig_dep"));
  Printer* I = new Printer(I_stream,'$');

  string sep = "";
  auto comma = parameters.find(',');
  string src_dir = parameters.substr(0, comma);
  string bin_dir = parameters.substr(comma+1);

  for(int i=0; i<file->dependency_count(); i++)
    write_dep(I, file->dependency(i)->name(), sep, src_dir, bin_dir);
  for(int i=0; i<file->public_dependency_count(); i++)
    write_dep(I, file->public_dependency(i)->name(), sep, src_dir, bin_dir);
  for(int i=0; i<file->weak_dependency_count(); i++)
    write_dep(I, file->weak_dependency(i)->name(), sep, src_dir, bin_dir);

  if(sep != "")I->Print("\n");

  delete I;
  delete I_stream;

  return true;
}

int main(int argc, char** argv) {
  SwigDepGenerator generator;
  return google::protobuf::compiler::PluginMain(argc, argv, &generator);
}
