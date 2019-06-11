/*

   calin/tools/protoc-gen-calin-cpp-extensions/protoc-gen-calin-cpp-extensions.cpp -- Stephen Fegan -- 2017-04-14

   Procobuf compiler plugin for generating CPP extensions

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <google/protobuf/compiler/plugin.h>
#include "calin_cpp_extensions_generator.hpp"

int main(int argc, char** argv) {
  calin::tools::calin_cpp_extensions_generator::CalinCppExtensionsGenerator generator;
  return google::protobuf::compiler::PluginMain(argc, argv, &generator);
}
