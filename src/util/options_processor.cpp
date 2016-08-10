/*

   calin/util/options_processor.cpp -- Stephen Fegan -- 2016-08-05

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

#include <google/protobuf/descriptor.h>

#include <util/options_processor.hpp>
#include <util/string_to_protobuf.hpp>

using namespace calin::util::options_processor;

OptionsProcessor::
OptionsProcessor(google::protobuf::Message* message):
  default_message_(message->New()), message_(message)
{
  default_message_->CopyFrom(*message);
}

void OptionsProcessor::process_arguments(const std::vector<std::string>& args)
{
  bool has_program_name = false;
  bool processing = true;
  for(const auto& iarg: args)
  {
    if(not has_program_name) {
      program_name_ = iarg;
      has_program_name = true;
      continue;
    }

    if(not processing) {
      arguments_.emplace_back(iarg);
      continue;
    }

    if(iarg.size() == 0) {
      // Purposefully empty arguments are passed to user
      arguments_.emplace_back(iarg);
      continue;
    }

    if(iarg[0] == '-')
    {
      if(iarg.size() == 1) {
        // A single dash argument is passed to user
        arguments_.emplace_back(iarg);
        continue;
      }

      // By here we know we have >=2 chars and first is a dash
      std::string::size_type ikey = 1;
      if(iarg[ikey] == '-') {
        if(iarg.size() == 2) {
          // The first double dash found terminates argument processing and
          // is itself absorbed
          processing = false;
          continue;
        } else {
          // Double dash with something following, advance past 2nd dash
          ikey++;
        }
      }

      if(not std::isalpha(iarg[ikey]) and iarg[ikey]!='_') {
        // Pass weird non-protobuf valid names to user to handle
        unknown_options_.emplace_back(iarg);
        continue;
      }

      std::string::size_type ifind = std::string::npos;
      if(ikey < iarg.size()-1)ifind = iarg.find_first_not_of(
        "abcdefghijklmnopqrstuvwxyz"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "1234567890_.", ikey+1);

      std::string key;
      std::string val;
      bool has_val = false;
      if(ifind == std::string::npos) {
        key = iarg.substr(ikey);
      } else if (iarg[ifind] == '=') {
        key = iarg.substr(ikey, ifind-ikey);
        if(ifind < iarg.size()-1) {
          has_val = true;
          val = iarg.substr(ifind+1);
        }
      } else {
        // Weird non-protobuf character found, pass argument to user to handle
        unknown_options_.emplace_back(iarg);
        continue;
      }

      switch(process_one_argument(message_, key, has_val, val))
      {
      case OptionStatus::OK:
        // Yay!
        break;
      case OptionStatus::NOT_FOUND:
        unknown_options_.emplace_back(iarg);
        break;
      case OptionStatus::COULD_NOT_CONVERT_VALUE:
        problem_options_.emplace_back(iarg);
        break;
      }
    }
    else
    {
      arguments_.emplace_back(iarg);
      continue;
    }
  }
}

void OptionsProcessor::process_arguments(int argc, char** argv)
{
  process_arguments(std::vector<std::string>(argv,argv+argc));
}

std::string OptionsProcessor::usage()
{

}

OptionStatus OptionsProcessor::
process_one_argument(google::protobuf::Message* m,
  const std::string& key, bool has_val, const std::string& val,
  std::vector<const google::protobuf::FieldDescriptor*> f_path)
{
  std::string::size_type ifind = key.find_first_of('.');
  std::string field_name;
  if(ifind == std::string::npos)
    field_name = key;
  else
    field_name = key.substr(0,ifind);

  if(field_name.empty())return OptionStatus::NOT_FOUND;

  const google::protobuf::Descriptor* d = nullptr;
  if(f_path.empty())d = m->GetDescriptor();
  else d = f_path.back()->message_type();
  assert(d);

  const google::protobuf::FieldDescriptor* f = d->FindFieldByName(field_name);
  if(f == nullptr)return OptionStatus::NOT_FOUND;

  if(f->message_type() != nullptr) {
    // Sub-message type
    if(f->is_repeated())
      return OptionStatus::NOT_FOUND; // Repeated messages not handled
    if(ifind == std::string::npos)
      return OptionStatus::NOT_FOUND; // Must have a sub-message field_name
    f_path.push_back(f);
    return process_one_argument(m, key.substr(ifind+1), has_val, val, f_path);
  } else {
    // Repeated or singular POD type
    if(string_to_protobuf::string_to_protobuf_field(val, m, f, f_path))
      return OptionStatus::OK;
    else
      return OptionStatus::COULD_NOT_CONVERT_VALUE;
  }
  assert(0);
  return OptionStatus::NOT_FOUND;
}
