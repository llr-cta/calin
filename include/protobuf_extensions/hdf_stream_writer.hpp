/*

  calin/protobuf_extensions/protohdf_stream_writer.hpp -- Stephen Fegan -- 2025-02-04

  Base classes and utility function for HDF stream writers.

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

#pragma once 

#include <vector>
#include <memory>
#include <iostream>
#include <cstring>

#include <hdf5.h>

namespace calin { namespace protobuf_extensions { namespace hdf_streamer {

template<typename T> inline hid_t h5_datatype_selector(hid_t h5t_in = -1) { return H5Tcopy(h5t_in>0?h5t_in:-1); }
template<> inline hid_t h5_datatype_selector<double>(hid_t h5t_in) { return H5Tcopy(h5t_in>0?h5t_in:H5T_NATIVE_DOUBLE); }
template<> inline hid_t h5_datatype_selector<float>(hid_t h5t_in) { return H5Tcopy(h5t_in>0?h5t_in:H5T_NATIVE_FLOAT); }
template<> inline hid_t h5_datatype_selector<uint64_t>(hid_t h5t_in) { return H5Tcopy(h5t_in>0?h5t_in:H5T_NATIVE_UINT64); }
template<> inline hid_t h5_datatype_selector<int64_t>(hid_t h5t_in) { return H5Tcopy(h5t_in>0?h5t_in:H5T_NATIVE_INT64); }
template<> inline hid_t h5_datatype_selector<uint32_t>(hid_t h5t_in) { return H5Tcopy(h5t_in>0?h5t_in:H5T_NATIVE_UINT32); }
template<> inline hid_t h5_datatype_selector<int32_t>(hid_t h5t_in) { return H5Tcopy(h5t_in>0?h5t_in:H5T_NATIVE_INT32); }
template<> inline hid_t h5_datatype_selector<uint16_t>(hid_t h5t_in) { return H5Tcopy(h5t_in>0?h5t_in:H5T_NATIVE_UINT16); }
template<> inline hid_t h5_datatype_selector<int16_t>(hid_t h5t_in) { return H5Tcopy(h5t_in>0?h5t_in:H5T_NATIVE_INT16); }
template<> inline hid_t h5_datatype_selector<uint8_t>(hid_t h5t_in) { return H5Tcopy(h5t_in>0?h5t_in:H5T_NATIVE_UINT8); }
template<> inline hid_t h5_datatype_selector<int8_t>(hid_t h5t_in) { return H5Tcopy(h5t_in>0?h5t_in:H5T_NATIVE_INT8); }
template<> inline hid_t h5_datatype_selector<std::string>(hid_t h5t_in) { 
  hid_t string_datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_datatype, H5T_VARIABLE);
  return string_datatype;
}
template<> inline hid_t h5_datatype_selector<char*>(hid_t h5t_in) { 
  hid_t string_datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_datatype, H5T_VARIABLE);
  return string_datatype;
}
template<> inline hid_t h5_datatype_selector<const char*>(hid_t h5t_in) { 
  hid_t string_datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_datatype, H5T_VARIABLE);
  return string_datatype;
}

// 888888b.                              
// 888  "88b                             
// 888  .88P                             
// 8888888K.   8888b.  .d8888b   .d88b.  
// 888  "Y88b     "88b 88K      d8P  Y8b 
// 888    888 .d888888 "Y8888b. 88888888 
// 888   d88P 888  888      X88 Y8b.     
// 8888888P"  "Y888888  88888P'  "Y8888  

class HDFStreamWriterBase
{
public:
  HDFStreamWriterBase(const std::string& filename, const std::string& groupname, bool truncate, const std::string& messagetype, uint64_t cache_size = 1024);
  HDFStreamWriterBase(const HDFStreamWriterBase* parent, const std::string& groupname, const std::string& messagetype, uint64_t cache_size = 1024);
  ~HDFStreamWriterBase();
  
  hid_t gid() const { return h5g_; }
  void flush();
  uint64_t nrow() const { return nrow_; }
  uint64_t ncached() const { return ncached_; }
  bool cache_full() const { return h5f_!=-1 and ncached_>=cache_size_; }

  template<typename T> void write_attribute(const std::string& name, T value) {
    hid_t attribute_id = H5Aopen(h5g_, name.c_str(), H5P_DEFAULT);
    hid_t type_id = h5_datatype_selector<T>();
    if(attribute_id < 0) {
      hid_t space_id = H5Screate(H5S_SCALAR);
      attribute_id = H5Acreate(h5g_, name.c_str(), type_id, space_id, H5P_DEFAULT, H5P_DEFAULT);
      H5Sclose(space_id);
      if (attribute_id < 0) {
        H5Tclose(type_id);
        throw std::runtime_error("Failed to create attribute: " + name);
      }
    }
    if(H5Awrite(attribute_id, type_id, &value) < 0) {
      H5Tclose(type_id);
      H5Aclose(attribute_id);
      throw std::runtime_error("Failed to write attribute: " + name);
    }
    H5Tclose(type_id);
    H5Aclose(attribute_id);
  }

  void write_attribute(const std::string& name, const std::string& value) {
    write_attribute(name, value.c_str());
  }

  template<typename T> bool read_attribute(const std::string& name, T* value) {
    hid_t attribute_id = H5Aopen(h5g_, name.c_str(), H5P_DEFAULT);
    if(attribute_id < 0) {
      *value = T();
      return false;
    }
    hid_t type_id = h5_datatype_selector<T>();
    if(H5Aread(attribute_id, type_id, value) < 0) {
      H5Tclose(type_id);
      H5Aclose(attribute_id);
      *value = T();
      return false;
    }
    H5Tclose(type_id);
    H5Aclose(attribute_id);
    return true;
  }

  bool read_attribute(const std::string& name, std::string* value) {
    char* cstr = nullptr;
    if(read_attribute(name, &cstr)) {
      value->assign(cstr);
      ::free(cstr);
      return true;
    }
    return false;
  }

protected:
  void open_group(hid_t file_id, const std::string& groupname, const std::string& messagetype);
  void insert_row() { nrow_+=1; ncached_+=1; }

  const HDFStreamWriterBase* parent_ = nullptr;
  uint64_t cache_size_ = 0;
  hid_t h5f_ = -1;
  hid_t h5g_ = -1;
  uint64_t nrow_ = 0;
  uint64_t ncached_ = 0;
  hid_t h5d_nrow_ = -1;
};


// 8888888b.   .d88888b.  8888888b.  
// 888   Y88b d88P" "Y88b 888  "Y88b 
// 888    888 888     888 888    888 
// 888   d88P 888     888 888    888 
// 8888888P"  888     888 888    888 
// 888        888     888 888    888 
// 888        Y88b. .d88P 888  .d88P 
// 888         "Y88888P"  8888888P"  

template<typename T> class DatasetWriter{
public:
  DatasetWriter(const HDFStreamWriterBase* base_ptr, const std::string dataset_name, T fill_value = T()):
    dataset_name_(dataset_name), datatype_(h5_datatype_selector<T>())
  {
    dataset_id_ = H5Dopen(base_ptr->gid(), dataset_name.c_str(), H5P_DEFAULT);
    if (dataset_id_ < 0) {
      hsize_t initial_size = 0;               // Start with an empty dataset
      hsize_t max_size = H5S_UNLIMITED;       // Unlimited extension

      // Create unlimited dataspace
      hid_t dataspace_id = H5Screate_simple(1, &initial_size, &max_size);

      // Set chunking properties (required for extensible datasets)
      hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
      hsize_t chunk_dims = 8192;
      H5Pset_chunk(plist_id, 1, &chunk_dims);
      H5Pset_shuffle(plist_id);
      H5Pset_deflate(plist_id, 2);

      // Create the dataset with chunking enabled
      dataset_id_ = H5Dcreate2(base_ptr->gid(), dataset_name.c_str(), datatype_, 
          dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);          
      H5Pclose(plist_id);
      H5Sclose(dataspace_id);

      if(dataset_id_ < 0) {
        throw std::runtime_error("Cannot create dataset: " + dataset_name);
      }
    }

    hid_t dataspace_id = H5Dget_space(dataset_id_);
    H5Sget_simple_extent_dims(dataspace_id, &nrow_, NULL);
    H5Sclose(dataspace_id);

    auto nfill = base_ptr->nrow();
    while(nrow_ < nfill) {
      write(fill_value);
      if(cache_.size() == 1024) flush(); 
    }
    if(not cache_.empty()) {
      flush();
    }

    if(nfill>0 and nrow_>nfill) {
      throw std::runtime_error("Too many rows in dataset \"" + dataset_name + "\" : " 
          + std::to_string(nrow_) + " > " + std::to_string(nfill));
    }
  }

  ~DatasetWriter() {
    flush();
    H5Dclose(dataset_id_);
    H5Tclose(datatype_);
  }

  uint64_t nrow() { return  nrow_; }
  uint64_t ncached() { return cache_.size(); }

  void write(const T& x) {
    cache_.push_back(x);
    nrow_ += 1;
  }

  void flush() {
    hsize_t start;
    hid_t space_id = H5Dget_space(dataset_id_);
    H5Sget_simple_extent_dims(space_id, &start, NULL);
    H5Sclose(space_id);

    hsize_t count = cache_.size();
    hsize_t after = start + count;
    if(after != nrow_) {
      throw std::runtime_error("Unexpected size in dataset \"" + dataset_name_ + "\" : " + std::to_string(after) + " != " + std::to_string(nrow_));
    }
    H5Dset_extent(dataset_id_, &after);

    hid_t file_space_id = H5Dget_space(dataset_id_);
    H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, &start, NULL, &count, NULL);
    hid_t mem_space_id = H5Screate_simple(1, &count, NULL);
    H5Dwrite(dataset_id_, datatype_, mem_space_id, file_space_id, H5P_DEFAULT, cache_.data());
    
    H5Sclose(mem_space_id);
    H5Sclose(file_space_id);

    cache_.clear();
  }

private:
  std::string dataset_name_;
  hid_t dataset_id_ = -1;
  hid_t datatype_ = -1;
  hsize_t nrow_ = 0;
  std::vector<T> cache_;
};


// 8888888b.   .d88888b.  8888888b.        d8888                                  
// 888   Y88b d88P" "Y88b 888  "Y88b      d88888                                  
// 888    888 888     888 888    888     d88P888                                  
// 888   d88P 888     888 888    888    d88P 888 888d888 888d888 8888b.  888  888 
// 8888888P"  888     888 888    888   d88P  888 888P"   888P"      "88b 888  888 
// 888        888     888 888    888  d88P   888 888     888    .d888888 888  888 
// 888        Y88b. .d88P 888  .d88P d8888888888 888     888    888  888 Y88b 888 
// 888         "Y88888P"  8888888P" d88P     888 888     888    "Y888888  "Y88888 
//                                                                            888 
//                                                                       Y8b d88P 
//                                                                        "Y88P"  

template<typename T> class ArrayDatasetWriter{
public:
  ArrayDatasetWriter(const HDFStreamWriterBase* base_ptr, const std::string dataset_name):
    dataset_name_(dataset_name)
  {
    // Create variable-length datatype
    hid_t datatype = h5_datatype_selector<T>();
    array_datatype_ = H5Tvlen_create(datatype);
    H5Tclose(datatype);

    dataset_id_ = H5Dopen(base_ptr->gid(), dataset_name.c_str(), H5P_DEFAULT);
    if (dataset_id_ < 0) {
      hsize_t initial_size = 0;               // Start with an empty dataset
      hsize_t max_size = H5S_UNLIMITED;       // Unlimited extension

      // Create unlimited dataspace
      hid_t dataspace_id = H5Screate_simple(1, &initial_size, &max_size);

      // Set chunking properties (required for extensible datasets)
      hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
      hsize_t chunk_dims = 8192;
      H5Pset_chunk(plist_id, 1, &chunk_dims);
      H5Pset_shuffle(plist_id);
      H5Pset_deflate(plist_id, 2);

      // Create the dataset with chunking enabled
      dataset_id_ = H5Dcreate2(base_ptr->gid(), dataset_name.c_str(), array_datatype_, 
          dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
                              
      H5Pclose(plist_id);
      H5Sclose(dataspace_id);

      if(dataset_id_ < 0) {
        throw std::runtime_error("Cannot create dataset: " + dataset_name);
      }
    }

    hid_t dataspace_id = H5Dget_space(dataset_id_);
    H5Sget_simple_extent_dims(dataspace_id, &nrow_, NULL);
    H5Sclose(dataspace_id);

    auto nfill = base_ptr->nrow();
    std::vector<T> fill_value;
    while(nrow_ < nfill) {
      write(fill_value);
      if(cache_.size() == 1024) flush(); 
    }
    if(not cache_.empty()) {
      flush();
    }

    if(nfill>0 and nrow_>nfill) {
      throw std::runtime_error("Too many rows in dataset \"" + dataset_name + "\" : " 
          + std::to_string(nrow_) + " > " + std::to_string(nfill));
    }
  }

  ~ArrayDatasetWriter() {
    flush();
    H5Dclose(dataset_id_);
    H5Tclose(array_datatype_);
  }

  uint64_t nrow() { return  nrow_; }
  uint64_t ncached() { return cache_.size(); }

  template<typename Container, typename Extractor> void write(const Container& c, const Extractor& f) {
    std::vector<T> v;
    v.reserve(c.size());
    for(const auto& x : c) {
      v.emplace_back(f(x));
    }
    cache_.emplace_back(std::move(v));
    nrow_ += 1;
  }

  template<typename Container> void write(const Container& c) {
    std::vector<T> v(c.begin(), c.end());
    cache_.emplace_back(std::move(v));
    nrow_ += 1;
  }

  void flush() {
    std::vector<hvl_t> hvl;
    hvl.reserve(cache_.size());
    for(auto& x : cache_) { 
      hvl.emplace_back(hvl_t{x.size(), x.data()}); 
    }

    hsize_t start;
    hid_t space_id = H5Dget_space(dataset_id_);
    H5Sget_simple_extent_dims(space_id, &start, NULL);
    H5Sclose(space_id);

    hsize_t count = hvl.size();
    hsize_t after = start + count;
    if(after != nrow_) {
      throw std::runtime_error("Unexpected size in dataset \"" + dataset_name_ + "\" : " + std::to_string(after) + " != " + std::to_string(nrow_));
    }
    H5Dset_extent(dataset_id_, &after);

    hid_t file_space_id = H5Dget_space(dataset_id_);
    H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, &start, NULL, &count, NULL);
    hid_t mem_space_id = H5Screate_simple(1, &count, NULL);
    H5Dwrite(dataset_id_, array_datatype_, mem_space_id, file_space_id, H5P_DEFAULT, hvl.data());
    
    H5Sclose(mem_space_id);
    H5Sclose(file_space_id);

    cache_.clear();
  }

private:
  std::string dataset_name_;
  hid_t dataset_id_ = -1;
  hid_t array_datatype_ = -1;
  hsize_t nrow_ = 0;
  std::vector<std::vector<T> > cache_;
};


//  .d8888b.  888            d8b                   
// d88P  Y88b 888            Y8P                   
// Y88b.      888                                  
//  "Y888b.   888888 888d888 888 88888b.   .d88b.  
//     "Y88b. 888    888P"   888 888 "88b d88P"88b 
//       "888 888    888     888 888  888 888  888 
// Y88b  d88P Y88b.  888     888 888  888 Y88b 888 
//  "Y8888P"   "Y888 888     888 888  888  "Y88888 
//                                             888 
//                                        Y8b d88P 
//                                         "Y88P"  

class StringDatasetWriter{
public:
  StringDatasetWriter(const HDFStreamWriterBase* base_ptr, const std::string dataset_name);
  ~StringDatasetWriter();

  uint64_t nrow() { return  nrow_; }
  uint64_t ncached() { return cache_.size(); }

  void write(const std::string& x) {
    write(x.c_str());
  }

  void write(const char* x) {
    cache_.emplace_back(::strdup(x));
    nrow_ += 1;
  }

  void flush();

private:
  std::string dataset_name_;
  hid_t dataset_id_ = -1;
  hid_t string_datatype_ = -1;
  hsize_t nrow_ = 0;
  std::vector<char*> cache_;
};

class StringArrayDatasetWriter{
public:
  StringArrayDatasetWriter(const HDFStreamWriterBase* base_ptr, const std::string dataset_name);
  ~StringArrayDatasetWriter();

  uint64_t nrow() { return  nrow_; }
  uint64_t ncached() { return cache_.size(); }

  template<typename Container, typename Extractor> void write(const Container& c, const Extractor& f) {
    std::vector<char*> v;
    v.reserve(c.size());
    for(const auto& x : c) {
      v.emplace_back(::strdup(f(x).c_str()));
    }
    cache_.emplace_back(std::move(v));
    nrow_ += 1;
  }

  template<typename Container> void write(const Container& c) {
    std::vector<char*> v;
    for(const auto& x : c) {
      v.emplace_back(::strdup(x.c_str()));
    }
    cache_.emplace_back(std::move(v));
    nrow_ += 1;
  }

  void flush();

private:
  std::string dataset_name_;
  hid_t dataset_id_ = -1;
  hid_t array_datatype_ = -1;
  hsize_t nrow_ = 0;
  std::vector<std::vector<char*> > cache_;
};


// 888888b.            888                     
// 888  "88b           888                     
// 888  .88P           888                     
// 8888888K.  888  888 888888 .d88b.  .d8888b  
// 888  "Y88b 888  888 888   d8P  Y8b 88K      
// 888    888 888  888 888   88888888 "Y8888b. 
// 888   d88P Y88b 888 Y88b. Y8b.          X88 
// 8888888P"   "Y88888  "Y888 "Y8888   88888P' 
//                 888                         
//            Y8b d88P                         
//             "Y88P"                          

class BytesDatasetWriter{
public:
  BytesDatasetWriter(const HDFStreamWriterBase* base_ptr, const std::string dataset_name);
  ~BytesDatasetWriter();

  uint64_t nrow() { return  nrow_; }
  uint64_t ncached() { return cache_.size(); }

  void write(const std::string& x) {
    cache_.emplace_back(x);
    nrow_ += 1;
  }

  void flush();

private:
  std::string dataset_name_;
  hid_t dataset_id_ = -1;
  hid_t array_datatype_ = -1;
  hsize_t nrow_ = 0;
  std::vector<std::string> cache_;
};

class BytesArrayDatasetWriter{
public:
  BytesArrayDatasetWriter(const HDFStreamWriterBase* base_ptr, const std::string dataset_name);
  ~BytesArrayDatasetWriter();

  uint64_t nrow() { return  nrow_; }
  uint64_t ncached() { return cache_.size(); }

  template<typename Container, typename Extractor> void write(const Container& c, const Extractor& f) {
    std::vector<std::string> v;
    v.reserve(c.size());
    for(const auto& x : c) {
      v.emplace_back(f(x));
    }
    cache_.emplace_back(std::move(v));
    nrow_ += 1;
  }


  template<typename Container> void write(const Container& c) {
    std::vector<std::string> v;
    for(const auto& x : c) {
      v.emplace_back(x);
    }
    cache_.emplace_back(std::move(v));
    nrow_ += 1;
  }

  void flush();

private:
  std::string dataset_name_;
  hid_t dataset_id_ = -1;
  hid_t array_datatype_ = -1;
  hsize_t nrow_ = 0;
  std::vector<std::vector<std::string> > cache_;
};


// 888b     d888                   
// 8888b   d8888                   
// 88888b.d88888                   
// 888Y88888P888  8888b.  88888b.  
// 888 Y888P 888     "88b 888 "88b 
// 888  Y8P  888 .d888888 888  888 
// 888   "   888 888  888 888 d88P 
// 888       888 "Y888888 88888P"  
//                        888      
//                        888      
//                        888      

template<typename KeyWriter, typename ValueWriter> class MapWriter {
public:
  MapWriter(const HDFStreamWriterBase* base_ptr, const std::string field_name) {
    key_writer_ = std::make_unique<KeyWriter>(base_ptr, field_name+"::key");
    value_writer_ = std::make_unique<ValueWriter>(base_ptr, field_name);
  }

  ~MapWriter() {
    flush();
  }

  template<typename Container> void write(const Container& c) {
    key_writer_->write(c, [](const auto& x) { return x.first; });
    value_writer_->write(c, [](const auto& x) { return x.second; });
  }

  void flush() {
    key_writer_->flush();
    value_writer_->flush();
  }

private:
  std::unique_ptr<KeyWriter> key_writer_;
  std::unique_ptr<ValueWriter> value_writer_;
};


// 888b     d888                                                       
// 8888b   d8888                                                       
// 88888b.d88888                                                        
// 888Y88888P888  .d88b.  .d8888b  .d8888b   8888b.   .d88b.   .d88b.  
// 888 Y888P 888 d8P  Y8b 88K      88K          "88b d88P"88b d8P  Y8b 
// 888  Y8P  888 88888888 "Y8888b. "Y8888b. .d888888 888  888 88888888 
// 888   "   888 Y8b.          X88      X88 888  888 Y88b 888 Y8b.     
// 888       888  "Y8888   88888P'  88888P' "Y888888  "Y88888  "Y8888  
//                                                        888          
//                                                   Y8b d88P          
//                                                    "Y88P"           

template<typename Message> class MessageWriter {
public:
  MessageWriter(const HDFStreamWriterBase* base_ptr, const std::string field_name):
    base_ptr_(base_ptr), field_name_(field_name)
  {
    // We defer opening sub-messages until they are used to allow nested messages to
    // be handled without getting stuck in an infinite loop.
  }

  ~MessageWriter() {
    flush();
  }

  void write(const Message* m) {
    if(not message_writer_) {
      if(m) deferred_open_datasets();
      else return;
    }

    start_writer_->write(message_writer_->nrow());
    if(m) {
      count_writer_->write(1);
      message_writer_->write(*m);
    } else {
      count_writer_->write(0);
    }
  }

  void write(const Message& m) {
    write(&m);
  }

  template<typename Container> void write(const Container& c) {
    if(not message_writer_) {
      if(not c.empty()) deferred_open_datasets();
      else return;
    }
    start_writer_->write(message_writer_->nrow());
    count_writer_->write(c.size());
    for(const auto& m : c) {
      message_writer_->write(m);
    }
  }

  template<typename Container, typename Extractor> void write(const Container& c, const Extractor& f) {
    if(not message_writer_) {
      if(not c.empty()) deferred_open_datasets();
      else return;
    }
    start_writer_->write(message_writer_->nrow());
    count_writer_->write(c.size());
    for(const auto& m : c) {
      message_writer_->write(f(m));
    }
  }

  void flush() {
    if(message_writer_) {
      start_writer_->flush();
      count_writer_->flush();
      message_writer_->flush();
    }
  }

private:
  void deferred_open_datasets() {
    message_writer_.reset(Message::__NewHDFStreamWriter(base_ptr_, field_name_));
    start_writer_ = std::make_unique<DatasetWriter<uint64_t> >(base_ptr_, field_name_+"::start", message_writer_->nrow());
    count_writer_ = std::make_unique<DatasetWriter<uint64_t> >(base_ptr_, field_name_+"::count", 0);
  }

  const HDFStreamWriterBase* base_ptr_ = nullptr;
  const std::string field_name_;

  std::unique_ptr<typename Message::stream_writer> message_writer_;
  std::unique_ptr<DatasetWriter<uint64_t> > start_writer_;
  std::unique_ptr<DatasetWriter<uint64_t> > count_writer_ ;
};

} } } // namespace calin::protobuf_extensions::hdf_streamer
