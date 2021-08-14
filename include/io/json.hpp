/*

   calin/util/json.hpp -- Stephen Fegan -- 2021-08-08

   JSON IO functions

   Copyright 2021, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#pragma once

#include <string>
#include <google/protobuf/message.h>

namespace calin { namespace io { namespace json {

std::string escape_json(const std::string &s);

template<typename T> inline std::string json_value(T x) {
  return std::to_string(x);
}

inline std::string json_value(bool b) {
  return b ? "true" : "false";
}

inline std::string json_string_value(const std::string& s) {
  return "\"" + escape_json(s) + "\"";
}

inline std::string json_value(const std::string& s) {
  return json_string_value(s);
}

inline std::string json_value(const char* s) {
  return json_string_value(s);
}

template<typename T> inline std::string json_for_single_value(const std::string key, T val) {
  return "{\n \"" + key + "\": " + json_value(val) + "\n}\n";
}

std::string json_for_dictionary(const std::vector<std::pair<std::string,std::string> >& keyval);

std::string encode_protobuf_to_json_string(const google::protobuf::Message& message);
void save_protobuf_to_json_file(const std::string& filename,
  const google::protobuf::Message& message);
inline void save_protobuf_to_json_file(const google::protobuf::Message& message,
    const std::string& filename) {
  save_protobuf_to_json_file(filename, message);
}

void decode_protobuf_from_json_string(const std::string& buffer_str,
  google::protobuf::Message* message);
inline void decode_protobuf_from_json_string(google::protobuf::Message* message,
    const std::string& buffer_str){
  decode_protobuf_from_json_string(buffer_str, message);
}

void load_protobuf_from_json_file(const std::string& filename,
  google::protobuf::Message* message);
inline void load_protobuf_from_json_file(google::protobuf::Message* message,
    const std::string& filename) {
  load_protobuf_from_json_file(filename, message);
}


} } } // namespace calin::io::json
