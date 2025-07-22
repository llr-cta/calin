/*

  calin/protobuf_extensions/hdf_stream_reader.hpp -- Stephen Fegan -- 2025-02-07

  Base classes and utility function for HDF stream readers.

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
#include <hdf5.h>
#include <iostream>

#include <protobuf_extensions/hdf_stream_writer.hpp>

namespace calin { namespace protobuf_extensions { namespace hdf_streamer {

// 888888b.                              
// 888  "88b                             
// 888  .88P                             
// 8888888K.   8888b.  .d8888b   .d88b.  
// 888  "Y88b     "88b 88K      d8P  Y8b 
// 888    888 .d888888 "Y8888b. 88888888 
// 888   d88P 888  888      X88 Y8b.     
// 8888888P"  "Y888888  88888P'  "Y8888  

class HDFStreamReaderBase
{
public:
  HDFStreamReaderBase(const std::string& filename, const std::string& groupname, uint64_t cache_size = 1024);
  HDFStreamReaderBase(const HDFStreamReaderBase* parent, const std::string& groupname, uint64_t cache_size = 1024);
  ~HDFStreamReaderBase();

  hid_t gid() const { return h5g_; }
  uint64_t nrow() const { return nrow_; }

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
  void open_group(hid_t file_id, const std::string& groupname);
  bool cache_preload_required(uint64_t irow) const { 
    return h5f_>=0 and (irow<cache_start_ or irow>=cache_end_); }
  void set_cache(uint64_t start, uint64_t count) { 
    cache_start_=std::min(start,nrow_); cache_end_=std::min(cache_start_+count,nrow_); }

  const HDFStreamReaderBase* parent_ = nullptr;
  uint64_t cache_size_ = 0;
  hid_t h5f_ = -1;
  hid_t h5g_ = -1;
  uint64_t nrow_ = 0;
  uint64_t cache_start_ = 0;
  uint64_t cache_end_ = 0;
};


// 8888888b.   .d88888b.  8888888b.  
// 888   Y88b d88P" "Y88b 888  "Y88b 
// 888    888 888     888 888    888 
// 888   d88P 888     888 888    888 
// 8888888P"  888     888 888    888 
// 888        888     888 888    888 
// 888        Y88b. .d88P 888  .d88P 
// 888         "Y88888P"  8888888P"  

template<typename T> class DatasetReader{
public:
  DatasetReader(const HDFStreamReaderBase* base_ptr, const std::string dataset_name):
    dataset_name_(dataset_name), datatype_(h5_datatype_selector<T>())
  {
    dataset_id_ = H5Dopen(base_ptr->gid(), dataset_name.c_str(), H5P_DEFAULT);
    if(dataset_id_ > 0) {
      hid_t dataspace_id = H5Dget_space(dataset_id_);
      H5Sget_simple_extent_dims(dataspace_id, &nrow_, NULL);
      H5Sclose(dataspace_id);
    }
  }

  ~DatasetReader() {
    if(dataset_id_ > 0) {
      H5Dclose(dataset_id_);
    }
    H5Tclose(datatype_);
  }

  uint64_t nrow() { return  nrow_; }

  bool read(uint64_t irow, T* x) {
    if(irow<cache_start_ or irow>=cache_end_){ *x = T(); return false; }
    *x = cache_[irow-cache_start_];
    return true;
  }

  T read(uint64_t irow, bool& good) {
    if(irow<cache_start_ or irow>=cache_end_) { /* note good is not set false */ return T(); }
    good = true;
    return cache_[irow-cache_start_];
  }

  T read(uint64_t irow) {
    bool good;
    return read(irow, good);
  }

  bool preload(uint64_t start, uint64_t count) {
    if(dataset_id_ < 0) {
      return false;
    }

    cache_start_ = start;
    cache_count_ = count; 

    cache_end_   = std::min(cache_start_ + cache_count_, nrow_);
    cache_start_ = std::min(cache_start_, nrow_);
    cache_count_ = cache_end_ - cache_start_;

    cache_.resize(cache_count_);

    if(cache_count_==0) {
      return true;
    }

    hid_t file_space_id = H5Dget_space(dataset_id_);
    H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, &cache_start_, NULL, &cache_count_, NULL);

    hid_t mem_space_id = H5Screate_simple(1, &cache_count_, NULL);
    herr_t res = H5Dread(dataset_id_, datatype_, mem_space_id, file_space_id, H5P_DEFAULT, cache_.data());

    H5Sclose(mem_space_id);
    H5Sclose(file_space_id);

    return res>=0;
  }

private:
  std::string dataset_name_;
  hid_t dataset_id_ = -1;
  hid_t datatype_ = -1;
  hsize_t nrow_ = 0;
  hsize_t cache_start_ = 0;
  hsize_t cache_count_ = 0;
  hsize_t cache_end_ = 0;
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

template<typename T> class ArrayDatasetReader{
public:
  ArrayDatasetReader(const HDFStreamReaderBase* base_ptr, const std::string dataset_name):
    dataset_name_(dataset_name)
  {
    // Create variable-length datatype
    hid_t datatype = h5_datatype_selector<T>();
    array_datatype_ = H5Tvlen_create(datatype);
    H5Tclose(datatype);

    dataset_id_ = H5Dopen(base_ptr->gid(), dataset_name.c_str(), H5P_DEFAULT);
    if(dataset_id_ > 0) {
      hid_t dataspace_id = H5Dget_space(dataset_id_);
      H5Sget_simple_extent_dims(dataspace_id, &nrow_, NULL);
      H5Sclose(dataspace_id);
    }
  }

  ~ArrayDatasetReader() {
    if(dataset_id_ > 0) {
      H5Dclose(dataset_id_);
    }
    H5Tclose(array_datatype_);
    free_cache();
  }

  uint64_t nrow() { return  nrow_; }

  bool count(uint64_t irow, uint64_t& count) {
    if(irow<cache_start_ or irow>=cache_end_) { count = 0; return false; }
    count = cache_[irow-cache_start_].len;
    return true;
  }

  uint64_t count(uint64_t irow) {
    uint64_t _count = 0;
    count(irow, _count);
    return _count;
  }

  bool read_one(uint64_t irow, uint64_t ielement, T* t) {
    if(irow<cache_start_ or irow>=cache_end_) { *t = T(); return false; }
    hvl_t x = cache_[irow-cache_start_];
    if(ielement >= x.len) { *t = T(); return false; }
    *t = reinterpret_cast<T*>(x.p)[ielement];
    return true;
  }

  template<typename Container> bool read(uint64_t irow, Container* c) {
    if(irow<cache_start_ or irow>=cache_end_) { c->Clear(); return false; }
    hvl_t x = cache_[irow-cache_start_];
    auto count = x.len;
    auto* p = reinterpret_cast<T*>(x.p);
    c->Resize(count, T());
    std::copy(p, p+count, c->begin());
    return true;
  }

  bool preload(uint64_t start, uint64_t count) {
    if(dataset_id_ < 0) {
      return false;
    }

    cache_start_ = start;
    cache_count_ = count; 

    cache_end_   = std::min(cache_start_ + cache_count_, nrow_);
    cache_start_ = std::min(cache_start_, nrow_);
    cache_count_ = cache_end_ - cache_start_;

    free_cache();
    cache_.resize(cache_count_);

    if(cache_count_==0) {
      return true;
    }

    hid_t file_space_id = H5Dget_space(dataset_id_);
    H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, &cache_start_, NULL, &cache_count_, NULL);

    hid_t mem_space_id = H5Screate_simple(1, &cache_count_, NULL);    
    herr_t res = H5Dread(dataset_id_, array_datatype_, mem_space_id, file_space_id, H5P_DEFAULT, cache_.data());

    H5Sclose(mem_space_id);
    H5Sclose(file_space_id);

    return res>=0;
  }

private:
  void free_cache() { for(auto& x: cache_) { ::free(x.p); } }

  std::string dataset_name_;
  hid_t dataset_id_ = -1;
  hid_t array_datatype_ = -1;
  hsize_t nrow_ = 0;
  hsize_t cache_start_ = 0;
  hsize_t cache_count_ = 0;
  hsize_t cache_end_ = 0;
  std::vector<hvl_t> cache_;
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

class StringDatasetReader{
public:
  StringDatasetReader(const HDFStreamReaderBase* base_ptr, const std::string dataset_name);
  ~StringDatasetReader();

  uint64_t nrow() { return  nrow_; }

  bool read(uint64_t irow, std::string* s);
  bool preload(uint64_t start, uint64_t count);

private:
  void free_cache();

  std::string dataset_name_;
  hid_t dataset_id_ = -1;
  hid_t string_datatype_ = -1;
  hsize_t nrow_ = 0;
  hsize_t cache_start_ = 0;
  hsize_t cache_count_ = 0;
  hsize_t cache_end_ = 0;
  std::vector<char*> cache_;
};

class StringArrayDatasetReader{
public:
  StringArrayDatasetReader(const HDFStreamReaderBase* base_ptr, const std::string dataset_name);
  ~StringArrayDatasetReader();

  uint64_t nrow() { return  nrow_; }

  bool count(uint64_t irow, uint64_t& count) {
    if(irow<cache_start_ or irow>=cache_end_) { count = 0; return false; }
    count = cache_[irow-cache_start_].len;
    return true;
  }

  uint64_t count(uint64_t irow) {
    uint64_t _count = 0;
    count(irow, _count);
    return _count;
  }

  bool read_one(uint64_t irow, uint64_t ielement, std::string* s) {
    if(irow<cache_start_ or irow>=cache_end_) { *s = std::string(); return false; }
    hvl_t x = cache_[irow-cache_start_];
    if(ielement >= x.len) { *s = std::string(); return false; }
    *s = reinterpret_cast<char**>(x.p)[ielement];
    return true;
  }

  template<typename Container> bool read(uint64_t irow, Container* c) {
    if(irow<cache_start_ or irow>=cache_end_) { c->Clear(); return false; }
    hvl_t x = cache_[irow-cache_start_];
    auto count = x.len;
    char** p = reinterpret_cast<char**>(x.p);
    while(c->size() > count) {
      c->RemoveLast();
    }
    for(int i=0;i<c->size();++i) {
      *(c->Mutable(i)) = p[i];
    }
    for(uint64_t i=c->size();i<count;++i) {
      *(c->Add()) = p[i];
    }
    return true;
  }

  bool preload(uint64_t start, uint64_t count);

private:
  void free_cache();

  std::string dataset_name_;
  hid_t dataset_id_ = -1;
  hid_t array_datatype_ = -1;
  hsize_t nrow_ = 0;
  hsize_t cache_start_ = 0;
  hsize_t cache_count_ = 0;
  hsize_t cache_end_ = 0;
  std::vector<hvl_t> cache_;
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

class BytesDatasetReader{
public:
  BytesDatasetReader(const HDFStreamReaderBase* base_ptr, const std::string dataset_name);
  ~BytesDatasetReader();

  uint64_t nrow() { return  nrow_; }

  bool read(uint64_t irow, std::string* x);
  bool preload(uint64_t start, uint64_t count);

private:
  void free_cache();
  
  std::string dataset_name_;
  hid_t dataset_id_ = -1;
  hid_t byte_string_datatype_ = -1;
  hsize_t nrow_ = 0;
  hsize_t cache_start_ = 0;
  hsize_t cache_count_ = 0;
  hsize_t cache_end_ = 0;
  std::vector<hvl_t> cache_;
};

class BytesArrayDatasetReader{
public:
  BytesArrayDatasetReader(const HDFStreamReaderBase* base_ptr, const std::string dataset_name);
  ~BytesArrayDatasetReader();

  uint64_t nrow() { return  nrow_; }

  bool count(uint64_t irow, uint64_t& count) {
    if(irow<cache_start_ or irow>=cache_end_) { count = 0; return false; }
    count = cache_[irow-cache_start_].len;
    return true;
  }

  uint64_t count(uint64_t irow) {
    uint64_t _count = 0;
    count(irow, _count);
    return _count;
  }

  bool read_one(uint64_t irow, uint64_t ielement, std::string* s) {
    if(irow<cache_start_ or irow>=cache_end_) { *s = std::string(); return false; }
    hvl_t x = cache_[irow-cache_start_];
    if(ielement >= x.len) { *s = std::string(); return false; }
    hvl_t y = reinterpret_cast<hvl_t*>(x.p)[ielement];
    s->assign(reinterpret_cast<const char*>(y.p), y.len);
    return true;
  }

  template<typename Container> bool read(uint64_t irow, Container* c) {
    if(irow<cache_start_ or irow>=cache_end_) { c->Clear(); return false; }
    hvl_t x = cache_[irow-cache_start_];
    uint64_t count = x.len;
    hvl_t* p = reinterpret_cast<hvl_t*>(x.p);

    while(c->size() > count) {
      c->RemoveLast();
    }
    for(int i=0;i<c->size();++i) {
      hvl_t y = p[i];
      c->Mutable(i)->assign(reinterpret_cast<const char*>(y.p), y.len);
    }
    for(uint64_t i=c->size();i<count;++i) {
      hvl_t y = p[i];
      c->Add()->assign(reinterpret_cast<const char*>(y.p), y.len);
    }
    return true;
  }

  bool preload(uint64_t start, uint64_t count);

private:
  void free_cache();
  
  std::string dataset_name_;
  hid_t dataset_id_ = -1;
  hid_t array_datatype_ = -1;
  hsize_t nrow_ = 0;
  hsize_t cache_start_ = 0;
  hsize_t cache_count_ = 0;
  hsize_t cache_end_ = 0;
  std::vector<hvl_t> cache_;
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

template<typename KeyReader, typename ValueReader> class MapReader {
public:
  MapReader(const HDFStreamReaderBase* base_ptr, const std::string field_name) {
    key_reader_ = std::make_unique<KeyReader>(base_ptr, field_name+"::key");
    value_reader_ = std::make_unique<ValueReader>(base_ptr, field_name);
  }

  ~MapReader() {
    // nothing to see here
  }

  bool count(uint64_t irow, uint64_t& count) {
    return key_reader_->count(irow, count);
  }

  uint64_t count(uint64_t irow) {
    uint64_t _count = 0;
    count(irow, _count);
    return _count;
  }

  template<typename Container> bool read(uint64_t irow, Container* c) {
    uint64_t n;
    if(not key_reader_->count(irow, n)) {
      c->clear();
      return false;
    }
    for(unsigned i=0; i<n; ++i) {
      typename Container::key_type key;
      if(key_reader_->read_one(irow, i, &key)) {
        value_reader_->read_one(irow, i, &(*c)[key]);
      }
    }
    return true;
  }

  bool preload(uint64_t start, uint64_t count) {
    return key_reader_->preload(start,count) and value_reader_->preload(start,count);
  }

private:
  std::unique_ptr<KeyReader> key_reader_;
  std::unique_ptr<ValueReader> value_reader_;
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

template<typename Message> class MessageReader {
public:
  MessageReader(const HDFStreamReaderBase* base_ptr, const std::string field_name) {
    message_reader_.reset(Message::__NewHDFStreamReader(base_ptr, field_name));
    start_reader_ = std::make_unique<DatasetReader<uint64_t> >(base_ptr, field_name+"::start");
    count_reader_ = std::make_unique<DatasetReader<uint64_t> >(base_ptr, field_name+"::count");
  }

  ~MessageReader() {
    // nothing to see here
  }

  bool count(uint64_t irow, uint64_t& count) {
    return count_reader_->read(irow, &count);
  }

  uint64_t count(uint64_t irow) {
    uint64_t _count = 0;
    count(irow, _count);
    return _count;
  }

  bool read_one(uint64_t irow, uint64_t ielement, Message* m) {
    uint64_t mstart = 0;
    uint64_t mcount = 0;
    bool good = start_reader_->read(irow, &mstart) and count_reader_->read(irow, &mcount);
    if(not good or ielement>=mcount) { m->Clear(); return false; }
    good = message_reader_->read(mstart+ielement, m);
    return good;
  }

  bool read(uint64_t irow, Message* m) {
    return read_one(irow, 0, m);
  }

  template<typename Container> bool read(uint64_t irow, Container* c) {
    uint64_t mstart = 0;
    uint64_t mcount = 0;
    bool good = start_reader_->read(irow, &mstart) and count_reader_->read(irow, &mcount);
    if(good) {
      while(c->size() > mcount) {
        c->RemoveLast();
      }
      for(int i=0;i<c->size();++i) {
        good &= message_reader_->read(mstart+i, c->Mutable(i));
      }
      for(uint64_t i=c->size();i<mcount;++i) {
        good &= message_reader_->read(mstart+i, c->Add());
      }
    } else {
      c->Clear();
    }
    return good;
  }

  bool preload(uint64_t start, uint64_t count) {
    bool good = true;
    good &= start_reader_->preload(start, count);
    good &= count_reader_->preload(start, count);
    uint64_t mstart = 0;
    uint64_t mcount = 0;
    start_reader_->read(start, &mstart);
    for(uint64_t irow=start,icount=0; count_reader_->read(irow,&icount); ++irow) {
      mcount += icount;
    }
    good &= message_reader_->preload(mstart, mcount);
    return good;
  }

private:
std::unique_ptr<typename Message::stream_reader> message_reader_;
  std::unique_ptr<DatasetReader<uint64_t> > start_reader_;
  std::unique_ptr<DatasetReader<uint64_t> > count_reader_;
};

} } } // namespace calin::protobuf_extensions::hdf_streamer
