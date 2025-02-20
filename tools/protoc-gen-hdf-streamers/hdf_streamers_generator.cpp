/*

  calin/tools/protoc-gen-hdf-streaners/hdf_streamers_generator.cpp -- Stephen Fegan -- 2025-02-04

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

#include <iostream>
#include <sstream>
#include <numeric>
#include <memory>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/descriptor.pb.h>
#include <google/protobuf/io/printer.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include <util/string.hpp>
#include <calin.pb.h>

#include "hdf_streamers_generator.hpp"

using std::string;
using namespace calin::util::string;
using namespace calin::tools::hdf_streamers_generator;
using google::protobuf::io::Printer;

HDFStreamersGenerator::~HDFStreamersGenerator()
{
  // nothing to see here
}

string pb_to_cpp_filename(string pb_filename, string extension = ".pb.cc")
{
  return pb_filename.substr(0,pb_filename.find_last_of('.')) + extension;
}

string pb_to_hpp_filename(string pb_filename, string extension = ".pb.h")
{
  return pb_to_cpp_filename(pb_filename, extension);
}

bool has_mutable(const google::protobuf::FieldDescriptor *f)
{
  return f->is_repeated() 
    or f->cpp_type()==google::protobuf::FieldDescriptor::CPPTYPE_STRING
    or f->cpp_type()==google::protobuf::FieldDescriptor::CPPTYPE_MESSAGE;
}

bool is_pod(google::protobuf::FieldDescriptor::Type type)
{
  switch(type)
  {
  case google::protobuf::FieldDescriptor::TYPE_DOUBLE:    // fall through
  case google::protobuf::FieldDescriptor::TYPE_FLOAT:     // fall through
  case google::protobuf::FieldDescriptor::TYPE_INT64:     // fall through
  case google::protobuf::FieldDescriptor::TYPE_UINT64:    // fall through
  case google::protobuf::FieldDescriptor::TYPE_INT32:     // fall through
  case google::protobuf::FieldDescriptor::TYPE_FIXED64:   // fall through
  case google::protobuf::FieldDescriptor::TYPE_FIXED32:   // fall through
  case google::protobuf::FieldDescriptor::TYPE_UINT32:    // fall through
  case google::protobuf::FieldDescriptor::TYPE_SFIXED32:  // fall through
  case google::protobuf::FieldDescriptor::TYPE_SFIXED64:  // fall through
  case google::protobuf::FieldDescriptor::TYPE_SINT32:    // fall through
  case google::protobuf::FieldDescriptor::TYPE_SINT64:    // fall through
  case google::protobuf::FieldDescriptor::TYPE_BOOL:      // fall through
  case google::protobuf::FieldDescriptor::TYPE_ENUM:      // fall through
  case google::protobuf::FieldDescriptor::TYPE_STRING:    // fall through
  case google::protobuf::FieldDescriptor::TYPE_BYTES:     // fall through
    return true;
  default:
    return false;
  }
  assert(0);
  return false;
}

std::string the_namespace(const google::protobuf::FileDescriptor* d) 
{
  auto cpts = split(d->package(),'.');
  return join(cpts, "::");
}

std::string message_name(const google::protobuf::Descriptor* d) 
{
  std::string name = d->name();
  while(d->containing_type()) { d=d->containing_type(); name = d->name() + "_" + name; }
  return name;
}

std::string message_fullname(const google::protobuf::Descriptor* d) 
{
  auto cpts = split(d->full_name(),'.');
  return join(cpts, "::");
}

std::string enum_fullname(const google::protobuf::EnumDescriptor* d) 
{
  auto cpts = split(d->full_name(),'.');
  return join(cpts, "::");
}

std::string stream_writer_name(const google::protobuf::Descriptor* d) 
{
  return message_name(d) + "_StreamWriter";
}

std::string stream_writer_fullname(const google::protobuf::Descriptor* d) 
{
  return the_namespace(d->file()) + "::" + stream_writer_name(d);
}

std::string hdf_stream_writer_name(const google::protobuf::Descriptor* d) 
{
  return stream_writer_name(d) + "_HDFImpl";
}

std::string hdf_stream_writer_fullname(const google::protobuf::Descriptor* d) 
{
  return the_namespace(d->file()) + "::" + hdf_stream_writer_name(d);
}

std::string stream_reader_name(const google::protobuf::Descriptor* d) 
{
  return message_name(d) + "_StreamReader";
}

std::string stream_reader_fullname(const google::protobuf::Descriptor* d) 
{
  return the_namespace(d->file()) + "::" + stream_reader_name(d);
}

std::string hdf_stream_reader_name(const google::protobuf::Descriptor* d) 
{
  return stream_reader_name(d) + "_HDFImpl";
}

std::string hdf_stream_reader_fullname(const google::protobuf::Descriptor* d) 
{
  return the_namespace(d->file()) + "::" + hdf_stream_reader_name(d);
}

std::string dsw_name(const google::protobuf::FieldDescriptor* f, const std::string& suffix="")
{
  return f->name() + "_" + std::to_string(f->number()) + suffix + "_";
}

std::string dsw_type(const google::protobuf::FieldDescriptor* f, bool is_array = false, const std::string& doer="Writer")
{
  if(f->is_map()) {
    const auto* dmap = f->message_type();
    const auto* kf = dmap->FindFieldByName("key");
    const auto* vf = dmap->FindFieldByName("value");
    return "Map"+doer+"<" + dsw_type(kf,/* is_array= */ true, doer) + "," + dsw_type(vf,/* is_array= */ true, doer) + " >";
  }

  const google::protobuf::FieldOptions* fopt = &f->options();
  const calin::FieldOptions* cfo = &fopt->GetExtension(calin::CFO);
  std::string i32t;
  switch(cfo->int32_type()) {
  case calin::FieldOptions::INT_32:
  default:
    i32t = "int32_t";
    break;
  case calin::FieldOptions::INT_16:
    i32t = "int16_t";
    break;
  case calin::FieldOptions::INT_8:
    i32t = "int8_t";
    break;
  }

  is_array |= f->is_repeated();
  switch(f->type())
  {
  case google::protobuf::FieldDescriptor::TYPE_DOUBLE:
    return is_array ?  "ArrayDataset"+doer+"<double>" : "Dataset"+doer+"<double>";
  case google::protobuf::FieldDescriptor::TYPE_FLOAT:
    return is_array ?  "ArrayDataset"+doer+"<float>" : "Dataset"+doer+"<float>";
  case google::protobuf::FieldDescriptor::TYPE_INT64:     // fall through
  case google::protobuf::FieldDescriptor::TYPE_FIXED64:   // fall through
  case google::protobuf::FieldDescriptor::TYPE_SFIXED64:  // fall through
  case google::protobuf::FieldDescriptor::TYPE_SINT64:
    return is_array ?  "ArrayDataset"+doer+"<int64_t>" : "Dataset"+doer+"<int64_t>";
  case google::protobuf::FieldDescriptor::TYPE_UINT64:    // fall through
    return is_array ?  "ArrayDataset"+doer+"<uint64_t>" : "Dataset"+doer+"<uint64_t>";
  case google::protobuf::FieldDescriptor::TYPE_INT32:     // fall through
  case google::protobuf::FieldDescriptor::TYPE_FIXED32:   // fall through
  case google::protobuf::FieldDescriptor::TYPE_SFIXED32:  // fall through
  case google::protobuf::FieldDescriptor::TYPE_SINT32:    // fall through
  case google::protobuf::FieldDescriptor::TYPE_ENUM:      
    return is_array ?  "ArrayDataset"+doer+"<"+i32t+">" : "Dataset"+doer+"<"+i32t+">";
  case google::protobuf::FieldDescriptor::TYPE_UINT32:    // fall through
  case google::protobuf::FieldDescriptor::TYPE_BOOL:
    return is_array ?  "ArrayDataset"+doer+"<u"+i32t+">" : "Dataset"+doer+"<u"+i32t+">";
  case google::protobuf::FieldDescriptor::TYPE_STRING:
    return is_array ?  "StringArrayDataset"+doer+"" : "StringDataset"+doer+"";
  case google::protobuf::FieldDescriptor::TYPE_BYTES:
    return is_array ?  "BytesArrayDataset"+doer+"" : "BytesDataset"+doer+"";
  case google::protobuf::FieldDescriptor::TYPE_MESSAGE:
    return "Message"+doer+"<"+message_fullname(f->message_type())+">";
  default:
    return "void";
  }
  assert(0);
  return "void";
}

std::string dsr_type(const google::protobuf::FieldDescriptor* f, bool is_array = false)
{
  return dsw_type(f, is_array, "Reader");
}

std::string oneof_name(const google::protobuf::OneofDescriptor* f, const std::string& suffix="")
{
  return f->name() + suffix + "_";
}

std::string oneof_dsw_type(const google::protobuf::OneofDescriptor* f)
{
  return "DatasetWriter<int>";
}

std::string oneof_dsr_type(const google::protobuf::OneofDescriptor* f)
{
  return "DatasetReader<int>";
}

void open_namespace(const google::protobuf::FileDescriptor* d,
  google::protobuf::io::Printer& printer) 
{
  auto cpts = split(d->package(),'.');
  for(auto cpt: cpts) {
    printer.Print("namespace $name$ { ", "name", cpt);
  }
  printer.Print("\n");
  printer.Indent();
}

void close_namespace(const google::protobuf::FileDescriptor* d,
  google::protobuf::io::Printer& printer) 
{
  auto cpts = split(d->package(),'.');
  std::string name = join(cpts, "::");
  printer.Outdent();
  for(auto cpt: cpts) {
    printer.Print("} ");
  }
  printer.Print("// namespace $name$\n", "name", name);
}


// 888    888                        888                  
// 888    888                        888                  
// 888    888                        888                  
// 8888888888  .d88b.   8888b.   .d88888  .d88b.  888d888 
// 888    888 d8P  Y8b     "88b d88" 888 d8P  Y8b 888P"   
// 888    888 88888888 .d888888 888  888 88888888 888     
// 888    888 Y8b.     888  888 Y88b 888 Y8b.     888     
// 888    888  "Y8888  "Y888888  "Y88888  "Y8888  888     

void generate_message_forward_declarations_headers(
  const google::protobuf::Descriptor* d,
  google::protobuf::io::Printer& printer)
{
  if(d->options().map_entry()) { return; }
  printer.Print(
    "class $stream_writer_name$;\n"
    "class $stream_reader_name$;\n",  
    "stream_writer_name", stream_writer_name(d),
    "stream_reader_name", stream_reader_name(d));
  for(int i=0;i<d->nested_type_count();++i) {
    generate_message_forward_declarations_headers(d->nested_type(i),printer);
  }
}

void generate_file_forward_declarations_headers(
  const google::protobuf::FileDescriptor * file,
  const string & parameter,
  google::protobuf::compiler::GeneratorContext *context,
  string * error)
{
  std::unique_ptr<google::protobuf::io::ZeroCopyOutputStream> output(
    context->OpenForInsert(pb_to_hpp_filename(file->name()), "includes"));
  google::protobuf::io::Printer printer(output.get(), '$');

  open_namespace(file, printer);
  printer.Print("// Forward declarations of stream writers / readers\n");
  for(int i=0;i<file->message_type_count();i++)
  {
    const google::protobuf::Descriptor* d = file->message_type(i);
    generate_message_forward_declarations_headers(d, printer);
  }
  close_namespace(file, printer);
}

void generate_message_stream_accessors(
  const google::protobuf::Descriptor* d,
  google::protobuf::compiler::GeneratorContext *context)
{
  if(d->options().map_entry()) { return; }
  if(true) {
    std::unique_ptr<google::protobuf::io::ZeroCopyOutputStream> output(
      context->OpenForInsert(pb_to_hpp_filename(d->file()->name()), "class_scope:"+d->full_name()));
    google::protobuf::io::Printer printer(output.get(), '$');

    printer.Print(
      "public:\n"
      "  using stream_writer = $stream_writer_name$;\n"
      "  using stream_reader = $stream_reader_name$;\n"
      "  static $stream_writer_name$* NewHDFStreamWriter(const std::string& filename, const std::string& groupname, bool truncate = false);\n" 
      "  static $stream_writer_name$* __NewHDFStreamWriter(const void* parent_ptr, const std::string& groupname);\n"
      "  static $stream_reader_name$* NewHDFStreamReader(const std::string& filename, const std::string& groupname);\n" 
      "  static $stream_reader_name$* __NewHDFStreamReader(const void* parent_ptr, const std::string& groupname);\n", 
      "stream_writer_name", stream_writer_name(d),
      "stream_reader_name", stream_reader_name(d));
  }
  for(int i=0;i<d->nested_type_count();++i) {
    generate_message_stream_accessors(d->nested_type(i), context);    
  }
}

void generate_file_stream_accessors(
  const google::protobuf::FileDescriptor * file,
  const string & parameter,
  google::protobuf::compiler::GeneratorContext *context,
  string * error)
{
  for(int i=0;i<file->message_type_count();i++)
  {
    const google::protobuf::Descriptor* d = file->message_type(i);
    generate_message_stream_accessors(d, context);
  }
}

void generate_message_stream_headers(
  const google::protobuf::Descriptor* d,
  google::protobuf::io::Printer& printer)
{
  if(d->options().map_entry()) { return; }
  printer.Print(
    "class $stream_writer_name$\n"
    "{\n"
    "public:\n"
    "  virtual ~$stream_writer_name$();\n"
    "  virtual uint64_t nrow() = 0;\n"
    "  virtual void write(const $message_name$& m) = 0;\n"
    "  virtual void flush() = 0;\n"
    "};\n"
    "class $stream_reader_name$\n"
    "{\n"
    "public:\n"
    "  virtual ~$stream_reader_name$();\n"
    "  virtual uint64_t nrow() = 0;\n"
    "  virtual bool read(uint64_t irow, $message_name$* m) = 0;\n"
    "  virtual bool preload(uint64_t start, uint64_t count) = 0;\n"
    "};\n",
    "message_name", message_fullname(d),
    "stream_writer_name", stream_writer_name(d),
    "stream_reader_name", stream_reader_name(d));

  for(int i=0;i<d->nested_type_count();++i) {
    generate_message_stream_headers(d->nested_type(i), printer);    
  }
}

void generate_file_stream_headers(
  const google::protobuf::FileDescriptor * file,
  const string & parameter,
  google::protobuf::compiler::GeneratorContext *context,
  string * error)
{
  std::unique_ptr<google::protobuf::io::ZeroCopyOutputStream> output(
    context->OpenForInsert(pb_to_hpp_filename(file->name()), "global_scope"));
  google::protobuf::io::Printer printer(output.get(), '$');

  open_namespace(file, printer);
  for(int i=0;i<file->message_type_count();i++)
  {
    const google::protobuf::Descriptor* d = file->message_type(i);
    generate_message_stream_headers(d, printer);
  }
  close_namespace(file, printer);    
}


// 888       888         d8b 888                    
// 888   o   888         Y8P 888                    
// 888  d8b  888             888                    
// 888 d888b 888 888d888 888 888888 .d88b.  888d888 
// 888d88888b888 888P"   888 888   d8P  Y8b 888P"   
// 88888P Y88888 888     888 888   88888888 888     
// 8888P   Y8888 888     888 Y88b. Y8b.     888     
// 888P     Y888 888     888  "Y888 "Y8888  888     
                                            
void generate_message_stream_writers_impl(
  const google::protobuf::Descriptor* d,
  google::protobuf::io::Printer& printer)
{
  if(d->options().map_entry()) { return; }

  // ===================================================================================================
  // ===================================================================================================
  // CLASS DEFINITION
  // ===================================================================================================
  // ===================================================================================================

  printer.Print(
    "class $hdf_stream_writer_name$:\n"
    "  public $stream_writer_name$,\n"
    "  public HDFStreamWriterBase\n"
    "{\n"
    "public:\n"
    "  $hdf_stream_writer_name$(const std::string& filename, const std::string& groupname, bool truncate):\n"
    "    $stream_writer_name$(),\n"
    "    HDFStreamWriterBase(filename, groupname, truncate, \"$message_name$\")\n"
    "  {\n"
    "    open_datasets();\n"
    "  }\n"
    "  $hdf_stream_writer_name$(const calin::protobuf_extensions::hdf_streamer::HDFStreamWriterBase* parent_ptr, const std::string& groupname):\n"
    "    $stream_writer_name$(),\n"
    "    HDFStreamWriterBase(parent_ptr, groupname, \"$message_name$\")\n"
    "  {\n"
    "    open_datasets();\n"
    "  }\n"
    "  virtual ~$hdf_stream_writer_name$();\n"
    "  uint64_t nrow() final;\n"
    "  void write(const $message_name$& m) final;\n"
    "  void flush() final;\n"
    "private:\n"
    "  void open_datasets();\n\n",
    "message_name", message_fullname(d),
    "stream_writer_name", stream_writer_fullname(d),
    "hdf_stream_writer_name", hdf_stream_writer_name(d));
  printer.Indent();
  for(int ifield=0; ifield<d->real_oneof_decl_count(); ++ifield) {
    const google::protobuf::OneofDescriptor* f = d->oneof_decl(ifield);
    printer.Print(
      "std::unique_ptr<$oneof_type$ > $oneof_name$;\n",
      "oneof_type", oneof_dsw_type(f),
      "oneof_name", oneof_name(f));
  }
  for(int ifield=0; ifield<d->field_count(); ++ifield) {
    const google::protobuf::FieldDescriptor* f = d->field(ifield);
    const google::protobuf::FieldOptions* fopt = &f->options();
    const calin::FieldOptions* cfo = &fopt->GetExtension(calin::CFO);
    if(cfo->dont_store())continue;
    printer.Print(
      "std::unique_ptr<$dsw_type$ > $dsw_name$;\n",
      "dsw_type", dsw_type(f),
      "dsw_name", dsw_name(f));
  }
  printer.Outdent();
  printer.Print("};\n\n");

  // ===================================================================================================
  // ===================================================================================================
  // NROW
  // ===================================================================================================
  // ===================================================================================================

  printer.Print(
    "uint64_t $hdf_stream_writer_name$::nrow() {\n"
    "  return HDFStreamWriterBase::nrow();\n"
    "}\n\n",
    "hdf_stream_writer_name", hdf_stream_writer_name(d));

  // ===================================================================================================
  // ===================================================================================================
  // DESTRUCTOR
  // ===================================================================================================
  // ===================================================================================================

  printer.Print(
    "$hdf_stream_writer_name$::~$hdf_stream_writer_name$() {\n"
    "  if(h5f_>=0) {\n"
    "    flush();\n"
    "  }\n",
    "hdf_stream_writer_name", hdf_stream_writer_name(d));
  printer.Print("}\n\n");

  // ===================================================================================================
  // ===================================================================================================
  // OPEN DATASETS
  // ===================================================================================================
  // ===================================================================================================

  printer.Print(
    "void $hdf_stream_writer_name$::open_datasets() {\n",
    "hdf_stream_writer_name", hdf_stream_writer_name(d));
  printer.Indent();
  for(int ifield=0; ifield<d->real_oneof_decl_count(); ++ifield) {
    const google::protobuf::OneofDescriptor* f = d->oneof_decl(ifield);
    printer.Print(
      "$oneof_name$ = std::make_unique<$oneof_type$ >(this, \"$name$::case\");\n",
      "oneof_type", oneof_dsw_type(f),
      "oneof_name", oneof_name(f),
      "name", f->name());
  }
  for(int ifield=0; ifield<d->field_count(); ++ifield) {
    const google::protobuf::FieldDescriptor* f = d->field(ifield);
    const google::protobuf::FieldOptions* fopt = &f->options();
    const calin::FieldOptions* cfo = &fopt->GetExtension(calin::CFO);
    if(cfo->dont_store())continue;
    auto name_f = f->name();
    printer.Print(
      "$dsw_name$ = std::make_unique<$dsw_type$ >(this, \"$name$\");\n",
      "dsw_type", dsw_type(f),
      "dsw_name", dsw_name(f),
      "name", f->name());
  }
  printer.Outdent();
  printer.Print("}\n\n");

  // ===================================================================================================
  // ===================================================================================================
  // WRITE MESSAGE
  // ===================================================================================================
  // ===================================================================================================

  printer.Print(
    "void $hdf_stream_writer_name$::write(const $message_name$& m) {\n",
    "message_name", message_fullname(d),
    "hdf_stream_writer_name", hdf_stream_writer_name(d));
  printer.Indent();
  for(int ifield=0; ifield<d->real_oneof_decl_count(); ++ifield) {
    const google::protobuf::OneofDescriptor* f = d->oneof_decl(ifield);
    printer.Print(
      "$oneof_name$->write(m.$name$_case());\n",
      "oneof_name", oneof_name(f),
      "name", f->name());
  }
  for(int ifield=0; ifield<d->field_count(); ++ifield) {
    const google::protobuf::FieldDescriptor* f = d->field(ifield);
    const google::protobuf::FieldOptions* fopt = &f->options();
    const calin::FieldOptions* cfo = &fopt->GetExtension(calin::CFO);
    if(cfo->dont_store())continue;
    auto dsw_name_f = dsw_name(f);
    if(is_pod(f->type()) or f->is_repeated()) {
      printer.Print(
        "$dsw_name$->write(m.$name$());\n",
        "dsw_name", dsw_name(f),
        "name", f->name());
    } else {
      printer.Print(
        "$dsw_name$->write(m.has_$name$() ? &m.$name$() : nullptr);\n",
        "dsw_name", dsw_name(f),
        "name", f->name());
    }
  }
  printer.Print(
    "HDFStreamWriterBase::insert_row();\n"
    "if(HDFStreamWriterBase::cache_full()) {\n"
    "  flush();\n"
    "}\n");
  printer.Outdent();
  printer.Print("}\n\n");

  // ===================================================================================================
  // ===================================================================================================
  // FLUSH CACHE
  // ===================================================================================================
  // ===================================================================================================

  printer.Print(
    "void $hdf_stream_writer_name$::flush() {\n",
    "hdf_stream_writer_name", hdf_stream_writer_name(d));
  printer.Indent();
  for(int ifield=0; ifield<d->real_oneof_decl_count(); ++ifield) {
    const google::protobuf::OneofDescriptor* f = d->oneof_decl(ifield);
    printer.Print(
      "$oneof_name$->flush();\n",
      "oneof_name", oneof_name(f));
  }
  for(int ifield=0; ifield<d->field_count(); ++ifield) {
    const google::protobuf::FieldDescriptor* f = d->field(ifield);
    const google::protobuf::FieldOptions* fopt = &f->options();
    const calin::FieldOptions* cfo = &fopt->GetExtension(calin::CFO);
    if(cfo->dont_store())continue;
    printer.Print(
      "$dsw_name$->flush();\n",
      "dsw_name", dsw_name(f));
  }
  printer.Print("HDFStreamWriterBase::flush();\n");
  printer.Outdent();
  printer.Print("}\n\n");

  for(int i=0;i<d->nested_type_count();++i) {
    generate_message_stream_writers_impl(d->nested_type(i), printer);    
  }
}

void generate_message_stream_writer_accessors_impl(
  const google::protobuf::Descriptor* d,
  google::protobuf::io::Printer& printer)
{
  if(d->options().map_entry()) { return; }

  printer.Print(
    "$stream_writer_name$::~$stream_writer_name$() { /* nothing to see here */ }\n"
    "$stream_writer_name$* $name$::NewHDFStreamWriter(const std::string& filename, const std::string& groupname, bool truncate) {\n"
    "  return new $hdf_stream_writer_name$(filename,groupname,truncate);\n"
    "}\n"
    "$stream_writer_name$* $name$::__NewHDFStreamWriter(const void* parent_ptr, const std::string& groupname) {\n"
    "  return new $hdf_stream_writer_name$(reinterpret_cast<const calin::protobuf_extensions::hdf_streamer::HDFStreamWriterBase*>(parent_ptr),groupname);\n"
    "}\n\n", 
    "name", message_name(d),
    "stream_writer_name", stream_writer_name(d),
    "hdf_stream_writer_name", hdf_stream_writer_name(d));

  for(int i=0;i<d->nested_type_count();++i) {
    generate_message_stream_writer_accessors_impl(d->nested_type(i), printer);    
  }
}
  
void generate_file_stream_writers_impl(
  const google::protobuf::FileDescriptor * file,
  const string & parameter,
  google::protobuf::compiler::GeneratorContext *context,
  string * error)
{
  std::unique_ptr<google::protobuf::io::ZeroCopyOutputStream> output(
    context->OpenForInsert(pb_to_cpp_filename(file->name()), "global_scope"));
  google::protobuf::io::Printer printer(output.get(), '$');

  printer.Print(
    "#include <memory>\n"
    "#include <protobuf_extensions/hdf_stream_writer.hpp>\n"
    "using namespace calin::protobuf_extensions::hdf_streamer;\n");
  printer.Print("namespace {\n");
  printer.Indent();

  for(int i=0;i<file->message_type_count();i++)
  {
    const google::protobuf::Descriptor* d = file->message_type(i);
    generate_message_stream_writers_impl(d, printer);
  }

  printer.Outdent();
  printer.Print("} // namespace anonymous\n\n");

  open_namespace(file, printer);
  for(int i=0;i<file->message_type_count();i++)
  {
    const google::protobuf::Descriptor* d = file->message_type(i);
    generate_message_stream_writer_accessors_impl(d, printer);    
  }
  close_namespace(file, printer);    
}


// 8888888b.                        888                  
// 888   Y88b                       888                  
// 888    888                       888                  
// 888   d88P .d88b.   8888b.   .d88888  .d88b.  888d888 
// 8888888P" d8P  Y8b     "88b d88" 888 d8P  Y8b 888P"   
// 888 T88b  88888888 .d888888 888  888 88888888 888     
// 888  T88b Y8b.     888  888 Y88b 888 Y8b.     888     
// 888   T88b "Y8888  "Y888888  "Y88888  "Y8888  888     

void generate_message_stream_readers_impl(
  const google::protobuf::Descriptor* d,
  google::protobuf::io::Printer& printer)
{
  if(d->options().map_entry()) { return; }

  // ===================================================================================================
  // ===================================================================================================
  // CLASS DEFINITION
  // ===================================================================================================
  // ===================================================================================================

  printer.Print(
    "class $hdf_stream_reader_name$:\n"
    "  public $stream_reader_name$,\n"
    "  public HDFStreamReaderBase\n"
    "{\n"
    "public:\n"
    "  $hdf_stream_reader_name$(const std::string& filename, const std::string& groupname):\n"
    "    $stream_reader_name$(),\n"
    "    HDFStreamReaderBase(filename, groupname)\n"
    "  {\n"
    "    open_datasets();\n"
    "  }\n"
    "  $hdf_stream_reader_name$(const calin::protobuf_extensions::hdf_streamer::HDFStreamReaderBase* parent_ptr, const std::string& groupname):\n"
    "    $stream_reader_name$(),\n"
    "    HDFStreamReaderBase(parent_ptr, groupname)\n"
    "  {\n"
    "    open_datasets();\n"
    "  }\n"
    "  virtual ~$hdf_stream_reader_name$();\n"
    "  uint64_t nrow() final;\n"
    "  bool read(uint64_t irow, $message_name$* m) final;\n"
    "  bool preload(uint64_t start, uint64_t count) final;\n"
    "private:\n"
    "  void open_datasets();\n\n",
    "message_name", message_fullname(d),
    "stream_reader_name", stream_reader_fullname(d),
    "hdf_stream_reader_name", hdf_stream_reader_name(d));
  printer.Indent();
  for(int ifield=0; ifield<d->real_oneof_decl_count(); ++ifield) {
    const google::protobuf::OneofDescriptor* f = d->oneof_decl(ifield);
    printer.Print(
      "std::unique_ptr<$oneof_type$ > $oneof_name$;\n",
      "oneof_type", oneof_dsr_type(f),
      "oneof_name", oneof_name(f));
  }
  for(int ifield=0; ifield<d->field_count(); ++ifield) {
    const google::protobuf::FieldDescriptor* f = d->field(ifield);
    const google::protobuf::FieldOptions* fopt = &f->options();
    const calin::FieldOptions* cfo = &fopt->GetExtension(calin::CFO);
    if(cfo->dont_store())continue;
    printer.Print(
      "std::unique_ptr<$dsr_type$ > $dsr_name$;\n",
      "dsr_type", dsr_type(f),
      "dsr_name", dsw_name(f));
  }
  printer.Outdent();
  printer.Print("};\n\n");

  // ===================================================================================================
  // ===================================================================================================
  // NROW
  // ===================================================================================================
  // ===================================================================================================

  printer.Print(
    "uint64_t $hdf_stream_reader_name$::nrow() {\n"
    "  return HDFStreamReaderBase::nrow();\n"
    "}\n\n",
    "hdf_stream_reader_name", hdf_stream_reader_name(d));

  // ===================================================================================================
  // ===================================================================================================
  // DESTRUCTOR
  // ===================================================================================================
  // ===================================================================================================

  printer.Print(
    "$hdf_stream_reader_name$::~$hdf_stream_reader_name$() {\n",
    "hdf_stream_reader_name", hdf_stream_reader_name(d));
  printer.Print("}\n\n");

  // ===================================================================================================
  // ===================================================================================================
  // OPEN DATASETS
  // ===================================================================================================
  // ===================================================================================================

  printer.Print(
    "void $hdf_stream_reader_name$::open_datasets() {\n"
    "  if(h5g_ < 0)return;\n",
    "hdf_stream_reader_name", hdf_stream_reader_name(d));
  printer.Indent();
  for(int ifield=0; ifield<d->real_oneof_decl_count(); ++ifield) {
    const google::protobuf::OneofDescriptor* f = d->oneof_decl(ifield);
    printer.Print(
      "$oneof_name$ = std::make_unique<$oneof_type$ >(this, \"$name$::case\");\n",
      "oneof_type", oneof_dsr_type(f),
      "oneof_name", oneof_name(f),
      "name", f->name());
  }
  for(int ifield=0; ifield<d->field_count(); ++ifield) {
    const google::protobuf::FieldDescriptor* f = d->field(ifield);
    const google::protobuf::FieldOptions* fopt = &f->options();
    const calin::FieldOptions* cfo = &fopt->GetExtension(calin::CFO);
    if(cfo->dont_store())continue;
    auto name_f = f->name();
    printer.Print(
      "$dsr_name$ = std::make_unique<$dsr_type$ >(this, \"$name$\");\n",
      "dsr_type", dsr_type(f),
      "dsr_name", dsw_name(f),
      "name", f->name());
  }
  printer.Outdent();
  printer.Print("}\n\n");

  // ===================================================================================================
  // ===================================================================================================
  // READ MESSAGE
  // ===================================================================================================
  // ===================================================================================================

  printer.Print(
    "bool $hdf_stream_reader_name$::read(uint64_t irow, $message_name$* m) {\n"
    "  if(h5g_ < 0)return false;\n"
    "  bool good = false;\n",
    "message_name", message_fullname(d),
    "hdf_stream_reader_name", hdf_stream_reader_name(d));
  printer.Indent();
  printer.Print(
    "if(HDFStreamReaderBase::cache_preload_required(irow)) {\n"
    "  preload(irow, cache_size_);\n"
    "}\n");
  for(int ifield=0; ifield<d->real_oneof_decl_count(); ++ifield) {
    const google::protobuf::OneofDescriptor* f = d->oneof_decl(ifield);
    printer.Print(
      "int $name$_case = 0;\n"
      "good |= $oneof_name$->read(irow,&$name$_case);\n",
      "oneof_name", oneof_name(f),
      "name", f->name());
  }
  for(int ifield=0; ifield<d->field_count(); ++ifield) {
    const google::protobuf::FieldDescriptor* f = d->field(ifield);
    const google::protobuf::FieldOptions* fopt = &f->options();
    const calin::FieldOptions* cfo = &fopt->GetExtension(calin::CFO);
    if(cfo->dont_store())continue;
    if(f->containing_oneof() != nullptr) {
      const google::protobuf::OneofDescriptor* oof = f->containing_oneof();
      printer.Print(
        "if($oneof_name$_case == $field_id$) {\n",
        "oneof_name", oof->name(),
        "field_id", std::to_string(f->number()));
      printer.Indent();
    }
    if(f->type() == google::protobuf::FieldDescriptor::TYPE_ENUM and not f->is_repeated()) {
      printer.Print(
        "m->set_$name$(static_cast<$enum_type$>($dsw_name$->read(irow, good)));\n",
        "dsw_name", dsw_name(f),
        "name", f->name(),
        "enum_type", enum_fullname(f->enum_type()));
    } else if(not has_mutable(f)) {
      printer.Print(
        "m->set_$name$($dsw_name$->read(irow, good));\n",
        "dsw_name", dsw_name(f),
        "name", f->name());
    } else if(is_pod(f->type()) or f->is_repeated()) {
      printer.Print(
        "good |= $dsw_name$->read(irow, m->mutable_$name$());\n",
        "dsw_name", dsw_name(f),
        "name", f->name());
    } else {
      printer.Print(
        "if($dsw_name$->count(irow)>0) {\n"
        "  good |= $dsw_name$->read(irow, m->mutable_$name$());\n"
        "} else {\n"
        "  m->clear_$name$();\n"
        "}\n",
        "dsw_name", dsw_name(f),
        "name", f->name());
    }
    if(f->containing_oneof() != nullptr) {
      printer.Outdent();
      printer.Print(
        "} else {\n"
        "  m->clear_$name$();\n"
        "}\n",
        "name", f->name());
    }
  }
  printer.Outdent();
  printer.Print(
    "  return good;\n"
    "}\n\n");

  // ===================================================================================================
  // ===================================================================================================
  // PRELOAD CACHE
  // ===================================================================================================
  // ===================================================================================================

  printer.Print(
    "bool $hdf_stream_reader_name$::preload(uint64_t start, uint64_t count) {\n"
    "  if(h5g_ < 0)return false;\n"
    "  bool good = false;\n",
    "hdf_stream_reader_name", hdf_stream_reader_name(d));
  printer.Indent();
  for(int ifield=0; ifield<d->real_oneof_decl_count(); ++ifield) {
    const google::protobuf::OneofDescriptor* f = d->oneof_decl(ifield);
    printer.Print(
      "good |= $oneof_name$->preload(start,count);\n",
      "oneof_name", oneof_name(f));
  }
  for(int ifield=0; ifield<d->field_count(); ++ifield) {
    const google::protobuf::FieldDescriptor* f = d->field(ifield);
    const google::protobuf::FieldOptions* fopt = &f->options();
    const calin::FieldOptions* cfo = &fopt->GetExtension(calin::CFO);
    if(cfo->dont_store())continue;
    printer.Print(
      "good |= $dsw_name$->preload(start,count);\n",
      "dsw_name", dsw_name(f));
  }
  printer.Print(
    "HDFStreamReaderBase::set_cache(start,count);\n"
    "return good;\n");
  printer.Outdent();
  printer.Print("}\n\n");

  for(int i=0;i<d->nested_type_count();++i) {
    generate_message_stream_readers_impl(d->nested_type(i), printer);    
  }
}

void generate_message_stream_reader_accessors_impl(
  const google::protobuf::Descriptor* d,
  google::protobuf::io::Printer& printer)
{
  if(d->options().map_entry()) { return; }

  printer.Print(
    "$stream_reader_name$::~$stream_reader_name$() { /* nothing to see here */ }\n"
    "$stream_reader_name$* $name$::NewHDFStreamReader(const std::string& filename, const std::string& groupname) {\n"
    "  return new $hdf_stream_reader_name$(filename,groupname);\n"
    "}\n"
    "$stream_reader_name$* $name$::__NewHDFStreamReader(const void* parent_ptr, const std::string& groupname) {\n"
    "  return new $hdf_stream_reader_name$(reinterpret_cast<const calin::protobuf_extensions::hdf_streamer::HDFStreamReaderBase*>(parent_ptr),groupname);\n"
    "}\n\n", 
    "name", message_name(d),
    "stream_reader_name", stream_reader_name(d),
    "hdf_stream_reader_name", hdf_stream_reader_name(d));

  for(int i=0;i<d->nested_type_count();++i) {
    generate_message_stream_reader_accessors_impl(d->nested_type(i), printer);    
  }
}

void generate_file_stream_readers_impl(
  const google::protobuf::FileDescriptor * file,
  const string & parameter,
  google::protobuf::compiler::GeneratorContext *context,
  string * error)
{
  std::unique_ptr<google::protobuf::io::ZeroCopyOutputStream> output(
    context->OpenForInsert(pb_to_cpp_filename(file->name()), "global_scope"));
  google::protobuf::io::Printer printer(output.get(), '$');

  printer.Print(
    "#include <protobuf_extensions/hdf_stream_reader.hpp>\n"
    "using namespace calin::protobuf_extensions::hdf_streamer;\n"
    "namespace {\n");
  printer.Indent();

  for(int i=0;i<file->message_type_count();i++)
  {
    const google::protobuf::Descriptor* d = file->message_type(i);
    generate_message_stream_readers_impl(d, printer);
  }

  printer.Outdent();
  printer.Print("} // namespace anonymous\n\n");

  open_namespace(file, printer);
  for(int i=0;i<file->message_type_count();i++)
  {
    const google::protobuf::Descriptor* d = file->message_type(i);
    generate_message_stream_reader_accessors_impl(d, printer);
  }
  close_namespace(file, printer);    
}

bool HDFStreamersGenerator::
Generate(const google::protobuf::FileDescriptor * file,
        const string & parameter,
        google::protobuf::compiler::GeneratorContext *context,
        string * error) const
{
  generate_file_forward_declarations_headers(file, parameter, context, error);
  generate_file_stream_accessors(file, parameter, context, error);
  generate_file_stream_headers(file, parameter, context, error);
  generate_file_stream_writers_impl(file, parameter, context, error);
  generate_file_stream_readers_impl(file, parameter, context, error);
  return true;
}
