/*

   calin/util/string_to_protobuf.hpp -- Stephen Fegan -- 2016-08-10

   String to and from protobuf

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

   Based on original, copyright 2006, Stephen Fegan, see notice below

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#pragma once

#include <string>

#include <google/protobuf/message.h>
#include <google/protobuf/descriptor.h>

namespace calin { namespace util { namespace string_to_protobuf {

bool string_to_protobuf_field(const std::string& s,
  google::protobuf::Message* m, const google::protobuf::FieldDescriptor* f);

} } } // namespace calin::util::string_to_protobuf
