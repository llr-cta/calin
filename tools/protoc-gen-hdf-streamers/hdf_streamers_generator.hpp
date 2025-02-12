/*

   calin/tools/protoc-gen-hdf-streamers/hdf_streamers_generator.hpp -- Stephen Fegan -- 2025-02-01

   Code generator for HDF streamers

   Copyright 2025, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <string>
#include <google/protobuf/compiler/code_generator.h>

namespace calin { namespace tools { namespace hdf_streamers_generator {

class HDFStreamersGenerator: public google::protobuf::compiler::CodeGenerator
{
 public:
  virtual ~HDFStreamersGenerator();
  virtual bool Generate(const google::protobuf::FileDescriptor * file,
                        const std::string & parameter,
                        google::protobuf::compiler::GeneratorContext *context,
                        std::string * error) const;
};

} } } // namespace calin::tools::hdf_streamers_generator
