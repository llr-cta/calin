/*

   calin/tools/protoc_gen_swig/swig_generator.hpp -- Stephen Fegan -- 2015-12-03

   Code generator for SWIG interface files

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

#include <string>
#include <google/protobuf/compiler/code_generator.h>

namespace calin { namespace tools { namespace protoc_gen_swig {

class SwigGenerator: public google::protobuf::compiler::CodeGenerator
{
 public:
  virtual ~SwigGenerator();
  virtual bool Generate(const google::protobuf::FileDescriptor * file,
                        const std::string & parameter,
                        google::protobuf::compiler::GeneratorContext *context,
                        std::string * error) const;
};

} } } // namespace calin::tools::protoc_gen_swig
