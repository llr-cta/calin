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
#include <map>
#include <string>

#include <google/protobuf/message.h>

namespace calin { namespace util { namespace options_processor {

enum class OptionHandlerResult {
  UNKNOWN_OPTION,
  EXCLUSIVE_OPTION_OK,
  SHARED_OPTION_OK,
  INVALID_OPTION_NAME,
  EXCLUSIVE_OPTION_INVALID_VALUE,
  SHARED_OPTION_INVALID_VALUE,
  TERMINATE_PROCESSING
};

struct OptionSpec
{
  std::string name;
  bool is_exclusive;
  bool takes_value;
  bool requires_value;
  std::string value_type;
  std::string value_default;
  std::string value_units;
  std::string description;
};

class OptionHandler
{
public:
  virtual ~OptionHandler();
  virtual OptionHandlerResult handle_option(const std::string& key,
    bool has_val, const std::string& val) = 0;
  virtual std::vector<OptionSpec> list_options() = 0;
};

#if 0
class SimpleOptionHandler: public OptionHandler
{
public:
  SimpleOptionHandler(): OptionHandler() { /* nothing to see here */ }
  SimpleOptionHandler(const std::string& key,
      OptionHandlerResult result = OptionHandlerResult::EXCLUSIVE_OPTION_OK):
    OptionHandler() { known_keys_[key] = result; }
  SimpleOptionHandler(std::map<std::string, OptionHandlerResult> known_keys):
    OptionHandler(), known_keys_(known_keys) { /* nothing to see here */ }
  virtual ~SimpleOptionHandler();
  OptionHandlerResult handle_option(const std::string& key,
    bool has_val, const std::string& val) override;
  std::vector<OptionSpec> list_options() override;
  bool was_option_handled() { return !keys_.empty(); }
  unsigned number_of_options_handled() { return keys_.size(); }
  const std::vector<std::string>& keys() const { return keys_; }
  const std::vector<bool> has_vals() const { return has_vals_; }
  const std::vector<std::string> vals() const { return vals_; }
protected:
  std::map<std::string, OptionHandlerResult> known_keys_;
  std::vector<std::string> keys_;
  std::vector<bool> has_vals_;
  std::vector<std::string> vals_;
};
#endif

class ProtobufOptionHandler: public OptionHandler
{
public:
  ProtobufOptionHandler(google::protobuf::Message* message);
  virtual ~ProtobufOptionHandler();
  OptionHandlerResult handle_option(const std::string& key,
    bool has_val, const std::string& val) override;
  std::vector<OptionSpec> list_options() override;
  void load_json_cfg(const std::string& json_file_name);
  void save_json_cfg(const std::string& json_file_name);
protected:
  std::vector<OptionSpec> r_list_options(const std::string& prefix,
    const google::protobuf::Message* m);
  google::protobuf::Message* default_message_ = nullptr;
  google::protobuf::Message* message_ = nullptr;
};

#if 0
class LoadCfgOptionHandler: public OptionHandler
{
public:
  LoadCfgOptionHandler(ProtobufOptionHandler* opt_proc,
      const std::string& key = "load_cfg"):
    OptionHandler(), opt_proc_(opt_proc), key_(key) {
      /* nothing to see here */ }
  virtual ~LoadCfgOptionHandler();
  OptionHandlerResult handle_option(const std::string& key,
    bool has_val, const std::string& val) override;
  std::vector<OptionSpec> list_options() override;
private:
  ProtobufOptionHandler* opt_proc_ = nullptr;
  std::string key_;
};
#endif

class OptionsProcessor
{
public:
  OptionsProcessor();
  OptionsProcessor(google::protobuf::Message* message);
  ~OptionsProcessor();
  void add_option_handler(OptionHandler* handler, bool adopt_handler = false);
  void process_arguments(const std::vector<std::string>& args,
    bool first_arg_is_program_name = true);
  void process_arguments(int argc, char** argv,
    bool first_arg_is_program_name = true);
  bool has_unknown_options() { return not unknown_options_.empty(); }
  const std::string& program_name() { return program_name_; }
  const std::vector<std::string>& unknown_options() { return unknown_options_; }
  const std::vector<std::string>& problem_options() { return problem_options_; }
  const std::vector<std::string>& arguments() { return arguments_; }
  std::vector<OptionSpec> list_options();
  std::string usage();

protected:
  ProtobufOptionHandler* priority_protobuf_handler_ = nullptr;
  std::vector<OptionHandler*> handlers_;
  std::vector<OptionHandler*> adopted_handlers_;

  std::string program_name_;
  std::vector<std::string> unknown_options_;
  std::vector<std::string> problem_options_;
  std::vector<std::string> arguments_;
  google::protobuf::Message* default_message_ = nullptr;
  google::protobuf::Message* message_ = nullptr;
};

} } } // namespace calin::util::options_processor
