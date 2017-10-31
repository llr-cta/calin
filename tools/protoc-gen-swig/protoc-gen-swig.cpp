/* 

   calin/tools/protoc-gen-swig.cpp -- Stephen Fegan -- 2015-12-03

   Procobuf compiler plugin for generating SWIG interface files

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

#include <google/protobuf/compiler/plugin.h>
#include "swig_generator.hpp"

int main(int argc, char** argv) {
  calin::tools::protoc_gen_swig::SwigGenerator generator;
  return google::protobuf::compiler::PluginMain(argc, argv, &generator);
}

