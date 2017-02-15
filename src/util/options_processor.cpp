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

#include <sstream>
#include <google/protobuf/descriptor.pb.h>

#include <calin.pb.h>
#include <util/options_processor.hpp>
#include <util/string_to_protobuf.hpp>
#include <util/file.hpp>
#include <util/string.hpp>

using namespace calin::util::options_processor;

OptionHandler::~OptionHandler()
{
  // nothing to see here
}

// ============================================================================
// ============================================================================

// SimpleOptionHandler

// ============================================================================
// ============================================================================

SimpleOptionHandler::~SimpleOptionHandler()
{
  // nothing to see here
}

OptionHandlerResult SimpleOptionHandler::handle_option(const std::string& key,
  bool has_val, const std::string& val)
{
  if(key != key_ and key != alt_key_)
    return OptionHandlerResult::UNKNOWN_OPTION;
  if(has_val)
    return OptionHandlerResult::KNOWN_OPTION_HAS_INVALID_VALUE;
  option_handled_ = true;
  return OptionHandlerResult::OPTION_OK;
}

std::vector<OptionSpec> SimpleOptionHandler::list_options()
{
  std::vector<OptionSpec> options;
  OptionSpec o;
  o.name           = key_;
  o.takes_value    = false;
  o.requires_value = false;
  o.value_type     = "";
  o.value_default  = "";
  o.value_units    = "";
  o.description    = description_;
  options.emplace_back(o);
  if(not alt_key_.empty()) {
    o.name         = alt_key_;
    options.emplace_back(o);
  }
  return options;
}

// ============================================================================
// ============================================================================

// ProtobufOptionHandler

// ============================================================================
// ============================================================================

ProtobufOptionHandler::
ProtobufOptionHandler(google::protobuf::Message* message):
  OptionHandler(), default_message_(message->New()), message_(message)
{
  assert(message);
  default_message_->CopyFrom(*message);
}

ProtobufOptionHandler::~ProtobufOptionHandler()
{
  delete default_message_;
}

OptionHandlerResult ProtobufOptionHandler::
handle_option(const std::string& key, bool has_val, const std::string& val)
{
  std::vector<const google::protobuf::FieldDescriptor*> f_path;

  const google::protobuf::Descriptor* d = message_->GetDescriptor();
  assert(d);

  std::string::size_type istart = 0;
  while(istart < key.size())
  {
    std::string::size_type ifind = key.find_first_of('.', istart);

    std::string field_name;
    if(ifind == std::string::npos)
      field_name = key.substr(istart);
    else
      field_name = key.substr(istart,ifind);

    if(field_name.empty())return OptionHandlerResult::UNKNOWN_OPTION;

    const google::protobuf::FieldDescriptor* f = d->FindFieldByName(field_name);
    if(f == nullptr)return OptionHandlerResult::UNKNOWN_OPTION;

    if(f->message_type() != nullptr) {
      // Sub-message type
      if(f->is_repeated()) // Repeated messages not handled
        return OptionHandlerResult::UNKNOWN_OPTION;
      if(ifind == std::string::npos) // No sub-field specified
        return OptionHandlerResult::UNKNOWN_OPTION;
      f_path.push_back(f);
      d = f->message_type();
      istart = ifind+1;
    } else {
      if(ifind != std::string::npos) // More sub-fields specified - oops
        return OptionHandlerResult::UNKNOWN_OPTION;
      // Repeated or singular POD type
      if(not has_val) {
        if(not f->is_repeated() and
            f->type() == google::protobuf::FieldDescriptor::TYPE_BOOL) {
          // As a special case, singular boolean types don't need a value
          google::protobuf::Message* m = message_;
          const google::protobuf::Reflection* r = m->GetReflection();
          for(const auto* if_path : f_path) {
            m = r->MutableMessage(m, if_path);
            assert(m);
            r = m->GetReflection();
            assert(r);
          }
          r->SetBool(m, f, true);
          return OptionHandlerResult::OPTION_OK;
        } else {
          return OptionHandlerResult::KNOWN_OPTION_HAS_INVALID_VALUE;
        }
      }
      if(string_to_protobuf::string_to_protobuf_field(val, message_, f, f_path))
        return OptionHandlerResult::OPTION_OK;
      else
        return OptionHandlerResult::KNOWN_OPTION_HAS_INVALID_VALUE;
    }
  }
  assert(0);
  return OptionHandlerResult::UNKNOWN_OPTION;
}

std::vector<OptionSpec> ProtobufOptionHandler::list_options()
{
  return r_list_options("", default_message_);
}

void ProtobufOptionHandler::load_json_cfg(const std::string& json_file_name)
{
  calin::util::file::load_protobuf_from_json_file(json_file_name, message_);
}

void ProtobufOptionHandler::save_json_cfg(const std::string& json_file_name)
{
  calin::util::file::save_protobuf_to_json_file(json_file_name, message_);
}

std::vector<OptionSpec> ProtobufOptionHandler::
r_list_options(const std::string& prefix, const google::protobuf::Message* m)
{
  std::vector<OptionSpec> options;
  const google::protobuf::Descriptor* d = m->GetDescriptor();
  const google::protobuf::Reflection* r = m->GetReflection();
  int nfield = d->field_count();
  for(int ifield=0; ifield<nfield; ifield++)
  {
    const google::protobuf::FieldDescriptor* f = d->field(ifield);
    if(f->message_type() != nullptr)
    {
      // Sub-message type
      if(f->is_repeated())continue; // Repeated messages not handled
      std::vector<OptionSpec> sub_options =
        r_list_options(prefix+f->name()+".", &r->GetMessage(*m, f));
      options.insert(options.end(), sub_options.begin(), sub_options.end());
    }
    else
    {
      OptionSpec o;
      o.name             = prefix+f->name();
      o.takes_value      = true;
      switch(f->type())
      {
      case google::protobuf::FieldDescriptor::TYPE_DOUBLE:   // fallthrough
      case google::protobuf::FieldDescriptor::TYPE_FLOAT:    // fallthrough
      case google::protobuf::FieldDescriptor::TYPE_SFIXED64: // fallthrough
      case google::protobuf::FieldDescriptor::TYPE_SINT64:   // fallthrough
      case google::protobuf::FieldDescriptor::TYPE_INT64:    // fallthrough
      case google::protobuf::FieldDescriptor::TYPE_FIXED64:  // fallthrough
      case google::protobuf::FieldDescriptor::TYPE_UINT64:   // fallthrough
      case google::protobuf::FieldDescriptor::TYPE_SFIXED32: // fallthrough
      case google::protobuf::FieldDescriptor::TYPE_SINT32:   // fallthrough
      case google::protobuf::FieldDescriptor::TYPE_INT32:    // fallthrough
      case google::protobuf::FieldDescriptor::TYPE_FIXED32:  // fallthrough
      case google::protobuf::FieldDescriptor::TYPE_UINT32:   // fallthrough
      case google::protobuf::FieldDescriptor::TYPE_STRING:   // fallthrough
      case google::protobuf::FieldDescriptor::TYPE_BYTES:
        o.requires_value = true;
        o.value_type     = f->cpp_type_name();
        break;
      case google::protobuf::FieldDescriptor::TYPE_BOOL:
        o.requires_value = f->is_repeated(); // Singular bools are special
        o.value_type     = "{true,false}";
        break;
      case google::protobuf::FieldDescriptor::TYPE_ENUM:
        o.requires_value = true;
        {
          std::string vals = "{";
          const google::protobuf::EnumDescriptor* ed = f->enum_type();
          assert(ed);
          for(int ivalue=0; ivalue<ed->value_count(); ivalue++) {
            if(ivalue) vals += ",";
            vals += ed->value(ivalue)->name();
          }
          vals += "}";
          o.value_type    = vals;
        }
        break;
      case google::protobuf::FieldDescriptor::TYPE_MESSAGE:
      case google::protobuf::FieldDescriptor::TYPE_GROUP:
      default:
        assert(0);
        break;
      }
      if(f->is_repeated())
        o.value_type = std::string("vector<")+o.value_type+">";
      string_to_protobuf::protobuf_field_to_string(o.value_default, m, f);
      o.value_units      = "";
      o.description      = "";
      const google::protobuf::FieldOptions* fopt { &f->options() };
      if(fopt->HasExtension(CFO))
      {
        o.description    = fopt->GetExtension(CFO).desc();
        o.value_units    = fopt->GetExtension(CFO).units();
      }
      options.emplace_back(o);
    }
  }
  return options;
}

// ============================================================================
// ============================================================================

// OptionsProcessor

// ============================================================================
// ============================================================================

OptionsProcessor::OptionsProcessor()
{
  // nothing to see here
}

OptionsProcessor::OptionsProcessor(google::protobuf::Message* message,
  bool add_help_option):
    priority_protobuf_handler_(new ProtobufOptionHandler(message))
{
  if(add_help_option) {
    help_handler_ = new SimpleOptionHandler("h", "help", "Print help message.");
    add_option_handler(help_handler_, true);
  }
  add_option_handler(priority_protobuf_handler_, true);
}

OptionsProcessor::~OptionsProcessor()
{
  for(auto* ihandler : adopted_handlers_)delete ihandler;
}

void OptionsProcessor::
add_option_handler(OptionHandler* handler, bool adopt_handler)
{
  handlers_.emplace_back(handler);
  if(adopt_handler)adopted_handlers_.emplace_back(handler);
}

#if 0
void OptionsProcessor::load_json_cfg(const std::string& json_file_name)
{
  assert(priority_protobuf_handler_);
  priority_protobuf_handler_->load_json_cfg(json_file_name);
}

void OptionsProcessor::save_json_cfg(const std::string& json_file_name)
{
  assert(priority_protobuf_handler_);
  priority_protobuf_handler_->save_json_cfg(json_file_name);
}
#endif

struct OptionDetails
{
  std::string arg;
  std::string key;
  bool has_val;
  std::string val;
  OptionHandlerResult status;
};

void OptionsProcessor::process_arguments(const std::vector<std::string>& args,
  bool first_arg_is_program_name)
{
  std::vector<OptionDetails> options;

  bool arg_is_program_name = first_arg_is_program_name;
  bool processing = true;
  for(const auto& iarg: args)
  {
    if(arg_is_program_name) {
      program_name_ = iarg;
      arg_is_program_name = false;
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

      OptionDetails option {
        iarg, "", false, "", OptionHandlerResult::UNKNOWN_OPTION };

      if(not std::isalpha(iarg[ikey]) and iarg[ikey]!='_') {
        // Pass weird non-protobuf valid names to user to handle
        option.status = OptionHandlerResult::INVALID_OPTION_NAME;
        continue;
      }

      std::string::size_type ifind = std::string::npos;
      if(ikey < iarg.size()-1)ifind = iarg.find_first_not_of(
        "abcdefghijklmnopqrstuvwxyz"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "1234567890_.", ikey+1);

      if(ifind == std::string::npos) {
        option.key = iarg.substr(ikey);
      } else if (iarg[ifind] == '=') {
        option.key = iarg.substr(ikey, ifind-ikey);
        if(ifind < iarg.size()-1) {
          option.has_val = true;
          option.val = iarg.substr(ifind+1);
        }
      } else {
        // Weird non-protobuf character found, pass argument to user to handle
        option.status = OptionHandlerResult::INVALID_OPTION_NAME;
        continue;
      }

      options.emplace_back(std::move(option));
    } // if(iarg[0] == '-')
    else
    {
      arguments_.emplace_back(iarg);
      continue;
    }
  } // for(const auto& iarg: args)

  for(auto* ihandler : handlers_)
  {
    for(auto& ioption : options)
    {
      switch(ioption.status)
      {
      case OptionHandlerResult::UNKNOWN_OPTION:
        // OK to continue processing this option through this handler
        break;
      case OptionHandlerResult::TERMINATE_PROCESSING:
        assert(0);
        // Fall through to continue
      case OptionHandlerResult::OPTION_OK:
      case OptionHandlerResult::KNOWN_OPTION_HAS_INVALID_VALUE:
      case OptionHandlerResult::INVALID_OPTION_NAME:
        // Skip this option
        continue;
      }

      OptionHandlerResult result =
        ihandler->handle_option(ioption.key, ioption.has_val, ioption.val);

      switch(result)
      {
      case OptionHandlerResult::INVALID_OPTION_NAME:
      case OptionHandlerResult::UNKNOWN_OPTION:
        // Leave the status as it was
        break;
      case OptionHandlerResult::OPTION_OK:
      case OptionHandlerResult::KNOWN_OPTION_HAS_INVALID_VALUE:
        ioption.status = result;
        break;
      case OptionHandlerResult::TERMINATE_PROCESSING:
        goto terminate_processing;
      }
    }
  }

terminate_processing:
  for(auto& ioption : options)
  {
    switch(ioption.status)
    {
    case OptionHandlerResult::INVALID_OPTION_NAME:
    case OptionHandlerResult::UNKNOWN_OPTION:
      unknown_options_.emplace_back(ioption.arg);
      break;
    case OptionHandlerResult::KNOWN_OPTION_HAS_INVALID_VALUE:
      problem_options_.emplace_back(ioption.arg);
      break;
    case OptionHandlerResult::TERMINATE_PROCESSING:
      assert(0);
    case OptionHandlerResult::OPTION_OK:
      break;
    }
  }
}

void OptionsProcessor::
process_arguments(int argc, char** argv, bool first_arg_is_program_name)
{
  process_arguments(std::vector<std::string>(argv,argv+argc),
    first_arg_is_program_name);
}

std::vector<OptionSpec> OptionsProcessor::list_options()
{
  std::vector<OptionSpec> options;
  for(auto* ihandler : handlers_)
  {
    std::vector<OptionSpec> sub_options = ihandler->list_options();
    options.insert(options.end(), sub_options.begin(), sub_options.end());
  }
  return options;
}

std::string OptionsProcessor::usage(unsigned width)
{
  constexpr unsigned indent   = 4;
  constexpr unsigned l1maxopt = 42;
  constexpr unsigned l1indent = 45;
  std::string usage_string;
  std::vector<OptionSpec> options = list_options();
  for(auto ioption : options)
  {
    if(not usage_string.empty())usage_string += "\n\n";

    std::string s;
    s += '-';
    s += ioption.name;
    if(ioption.takes_value) {
      //if(not ioption.requires_value)s += '[';
      s += '=';
      if(ioption.value_type=="string" or ioption.value_type=="bytes")s += '"';
      s += ioption.value_default;
      if(ioption.value_type=="string" or ioption.value_type=="bytes")s += '"';
      //if(!ioption.requires_value)s += ']';
      if(not ioption.value_units.empty()) {
        s += " [";
        s += ioption.value_units;
        s += "]";
      }
      s += " (";
      s += ioption.value_type;
      s += ")";
    }

    if(not ioption.description.empty()) {
      if(s.size() <= l1maxopt and l1indent <= width) {
        std::string desc =
          calin::util::string::reflow(ioption.description,
            width, indent, width-l1indent, 0);
        if(desc.find('\n') == std::string::npos) {
          s += std::string(l1indent-s.size(), ' ');
          s += desc;
        } else {
          s += "\n";
          s += calin::util::string::reflow(ioption.description, width, indent);          
        }
      } else {
        s += "\n";
        s += calin::util::string::reflow(ioption.description, width, indent);
      }
    }

    usage_string += s;
  }
  return usage_string;
}
