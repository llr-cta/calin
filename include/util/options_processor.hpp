/*

   calin/util/options_processor.hpp -- Stephen Fegan -- 2016-08-05

   Process command line options to protobuf

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

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
#include <util/options_processor.pb.h>

namespace calin { namespace util { namespace options_processor {

enum class OptionHandlerResult {
  INVALID_OPTION_NAME,
  UNKNOWN_OPTION,
  OPTION_OK,
  KNOWN_OPTION_HAS_INVALID_VALUE,
  TERMINATE_PROCESSING
};

struct OptionSpec
{
  std::string name;
  bool takes_value;
  bool requires_value;
  std::string value_type;
  std::string value_default;
  std::string value_units;
  std::string description;
};

#ifdef SWIG
} } } // namespace calin::util::options_processor
%template (Vector_OptionsProcessor_OptionSpec) std::vector<calin::util::options_processor::OptionSpec>;
namespace calin { namespace util { namespace options_processor {
#endif

class OptionHandler
{
public:
  virtual ~OptionHandler();
  virtual OptionHandlerResult handle_option(const std::string& key,
    bool has_val, const std::string& val) = 0;
  virtual std::vector<OptionSpec> list_options() = 0;
};

class SimpleOptionHandler: public OptionHandler
{
public:
  SimpleOptionHandler(const std::string& key, const std::string& help_text):
    OptionHandler(), key_(key), description_(help_text)
  { /* nothing to see here */ }
  SimpleOptionHandler(const std::string& key, const std::string& alt_key,
      const std::string& help_text):
    OptionHandler(), key_(key), alt_key_(alt_key), description_(help_text)
  { /* nothing to see here */ }
  virtual ~SimpleOptionHandler();
  OptionHandlerResult handle_option(const std::string& key,
    bool has_val, const std::string& val) override;
  std::vector<OptionSpec> list_options() override;
  bool was_option_handled() { return option_handled_; }
protected:
  std::string key_;
  std::string alt_key_;
  std::string description_;
  bool option_handled_ = false;
};

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
  OptionsProcessor(google::protobuf::Message* message,
    bool add_help_option = false);
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
  std::string usage(unsigned width = 80);
  bool help_requested() { return help_handler_ != nullptr
    and help_handler_->was_option_handled(); }
  const calin::ix::util::options_processor::CommandLineArguments& command_line_arguments() { return cla_; }
protected:
  ProtobufOptionHandler* priority_protobuf_handler_ = nullptr;
  std::vector<OptionHandler*> handlers_;
  std::vector<OptionHandler*> adopted_handlers_;
  SimpleOptionHandler* help_handler_ = nullptr;

  std::string program_name_;
  std::vector<std::string> unknown_options_;
  std::vector<std::string> problem_options_;
  std::vector<std::string> arguments_;
  google::protobuf::Message* default_message_ = nullptr;
  google::protobuf::Message* message_ = nullptr;

  calin::ix::util::options_processor::CommandLineArguments cla_;
};

} } } // namespace calin::util::options_processor
