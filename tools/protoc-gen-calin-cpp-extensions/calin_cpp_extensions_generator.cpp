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

}

bool CalinCppExtensionsGenerator::
Generate(const google::protobuf::FileDescriptor * file,
         const string & parameter,
         google::protobuf::compiler::GeneratorContext *context,
         string * error) const
{
  // Print IntegrateFrom functions for all messages that request them
  for(int i=0;i<file->message_type_count();i++)
  {
    const google::protobuf::Descriptor* d = file->message_type(i);
    const google::protobuf::MessageOptions* mopt = &d->options();
    if(!mopt->HasExtension(calin::CMO))continue;
    const calin::MessageOptions* cmo = &mopt->GetExtension(calin::CMO);
    if(cmo->message_integration_function() == calin::MessageOptions::MIF_NONE)
      continue;
    google::protobuf::scoped_ptr<google::protobuf::io::ZeroCopyOutputStream> output(
        context->OpenForInsert(pb_to_hpp_filename(file->name()),
        std::string("class_scope:")+d->full_name()));
    google::protobuf::io::Printer printer(output.get(), '$');
    printer.Print("void IntegrateFrom(const $type$& from);\n", "type", d->name());

    if(cmo->message_integration_function() == calin::MessageOptions::MIF_USER_GENERATED)
      continue;

  }

  return true;
}
