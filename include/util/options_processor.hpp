/*

   calin/util/options_processor.hpp -- Stephen Fegan -- 2016-08-05

   Process command line options to protobuf

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

#include <vector>
#include <string>

#include <google/protobuf/message.h>

namespace calin { namespace util { namespace options_processor {

enum class OptionStatus { OK, NOT_FOUND, COULD_NOT_CONVERT_VALUE };

class OptionsProcessor
{
public:
  OptionsProcessor(google::protobuf::Message* message);
  void load_json_cfg(const std::string& json_file_name);
  void save_json_cfg(const std::string& json_file_name);
  void process_arguments(const std::vector<std::string>& args,
    const std::string& load_cfg_opt = {},
    bool first_arg_is_program_name = true);
  void process_arguments(int argc, char** argv,
    const std::string& load_cfg_opt = {},
    bool first_arg_is_program_name = true);
  bool has_unknown_options() { return not unknown_options_.empty(); }
  const std::string& program_name() { return program_name_; }
  const std::vector<std::string>& unknown_options() { return unknown_options_; }
  const std::vector<std::string>& problem_options() { return problem_options_; }
  const std::vector<std::string>& arguments() { return arguments_; }
  std::string usage();

protected:
  OptionStatus process_one_argument(google::protobuf::Message* m,
    const std::string& key, bool has_val, const std::string& val,
    std::vector<const google::protobuf::FieldDescriptor*> f_path = {});

  std::string program_name_;
  std::vector<std::string> unknown_options_;
  std::vector<std::string> problem_options_;
  std::vector<std::string> arguments_;
  google::protobuf::Message* default_message_ = nullptr;
  google::protobuf::Message* message_ = nullptr;
};

} } } // namespace calin::util::options_processor
