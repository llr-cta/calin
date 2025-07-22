/*

  calin/src/protobuf_extensions/hdf_stream_reader.cpp -- Stephen Fegan -- 2025-02-97

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

#include <string>

#include <protobuf_extensions/hdf_stream_reader.hpp>

using namespace calin::protobuf_extensions::hdf_streamer;

// 888888b.                              
// 888  "88b                             
// 888  .88P                             
// 8888888K.   8888b.  .d8888b   .d88b.  
// 888  "Y88b     "88b 88K      d8P  Y8b 
// 888    888 .d888888 "Y8888b. 88888888 
// 888   d88P 888  888      X88 Y8b.     
// 8888888P"  "Y888888  88888P'  "Y8888  

HDFStreamReaderBase::
HDFStreamReaderBase(const std::string& filename, const std::string& groupname, uint64_t cache_size):
  cache_size_(cache_size)
{
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);
  
  h5f_ = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (h5f_ < 0)
    throw std::runtime_error("Failed to open file: " + filename);

  open_group(h5f_, groupname);
  if(h5g_ < 0) {
    throw std::runtime_error("Failed to open group: " + groupname);
  }
}

HDFStreamReaderBase::
HDFStreamReaderBase(const HDFStreamReaderBase* parent, const std::string& groupname, uint64_t cache_size):
  parent_(parent), cache_size_(cache_size)
{
  open_group(parent_->h5g_, groupname);
  // We don't check the error, this just means the sub-message isn't present
  // which allows us to read older versions of the file.
}

HDFStreamReaderBase::
~HDFStreamReaderBase()
{
  if (h5g_ >= 0) {
    H5Gclose(h5g_);
  }
  if(h5f_ >= 0) {
    H5Fclose(h5f_);
  }
}

void HDFStreamReaderBase::
open_group(hid_t file_id, const std::string& groupname)
{
  h5g_ = H5Gopen(file_id, groupname.c_str(), H5P_DEFAULT);
  if (h5g_ >= 0) {
    hid_t h5d_nrow = H5Aopen(h5g_, "nrow", H5P_DEFAULT);
    if (h5d_nrow < 0) {
      throw std::runtime_error("Failed to open attribute: nrow");
    }
    if(H5Aread(h5d_nrow, H5T_NATIVE_UINT64, &nrow_) < 0) {
      H5Aclose(h5d_nrow);
      throw std::runtime_error("Failed to read attribute: nrow");
    }
    H5Aclose(h5d_nrow);
  }
}

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

StringDatasetReader::
StringDatasetReader(const HDFStreamReaderBase* base_ptr, const std::string dataset_name):
  dataset_name_(dataset_name)
{
  // Create variable-length datatype
  string_datatype_ = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_datatype_, H5T_VARIABLE); // Variable-length strings

  dataset_id_ = H5Dopen(base_ptr->gid(), dataset_name.c_str(), H5P_DEFAULT);
  if(dataset_id_ >= 0) {
    hid_t dataspace_id = H5Dget_space(dataset_id_);
    H5Sget_simple_extent_dims(dataspace_id, &nrow_, NULL);
    H5Sclose(dataspace_id);
  }
}

StringDatasetReader::~StringDatasetReader()
{
  if(dataset_id_ > 0) {
    H5Dclose(dataset_id_);
  }
  H5Tclose(string_datatype_);
  free_cache();
}

bool StringDatasetReader::read(uint64_t irow, std::string* s)
{
  if(irow<cache_start_ or irow>=cache_end_){ *s = std::string(); return false; }
  *s = cache_[irow-cache_start_];
  return true;
}

bool StringDatasetReader::preload(uint64_t start, uint64_t count)
{
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
  herr_t res = H5Dread(dataset_id_, string_datatype_, mem_space_id, file_space_id, H5P_DEFAULT, cache_.data());

  H5Sclose(mem_space_id);
  H5Sclose(file_space_id);

  return res>=0;
}

void StringDatasetReader::free_cache()
{
  for(auto& x: cache_) { ::free(x); }
}

//  .d8888b.  888            d8b                          d8888                                  
// d88P  Y88b 888            Y8P                         d88888                                  
// Y88b.      888                                       d88P888                                  
//  "Y888b.   888888 888d888 888 88888b.   .d88b.      d88P 888 888d888 888d888 8888b.  888  888 
//     "Y88b. 888    888P"   888 888 "88b d88P"88b    d88P  888 888P"   888P"      "88b 888  888 
//       "888 888    888     888 888  888 888  888   d88P   888 888     888    .d888888 888  888 
// Y88b  d88P Y88b.  888     888 888  888 Y88b 888  d8888888888 888     888    888  888 Y88b 888 
//  "Y8888P"   "Y888 888     888 888  888  "Y88888 d88P     888 888     888    "Y888888  "Y88888 
//                                             888                                           888 
//                                        Y8b d88P                                      Y8b d88P 
//                                         "Y88P"                                        "Y88P"  

StringArrayDatasetReader::
StringArrayDatasetReader(const HDFStreamReaderBase* base_ptr, const std::string dataset_name):
  dataset_name_(dataset_name)
{
  // Create string datatype
  hid_t string_datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_datatype, H5T_VARIABLE); // Variable-length strings

  // Create variable-length datatype
  array_datatype_ = H5Tvlen_create(string_datatype);
  H5Tclose(string_datatype);

  dataset_id_ = H5Dopen(base_ptr->gid(), dataset_name.c_str(), H5P_DEFAULT);
  if (dataset_id_ >= 0) {
    hid_t dataspace_id = H5Dget_space(dataset_id_);
    H5Sget_simple_extent_dims(dataspace_id, &nrow_, NULL);
    H5Sclose(dataspace_id);
  }
}

StringArrayDatasetReader::~StringArrayDatasetReader() {
  if(dataset_id_ > 0) {
    H5Dclose(dataset_id_);
  }
  H5Tclose(array_datatype_);
  free_cache();
}

bool StringArrayDatasetReader::preload(uint64_t start, uint64_t count)
{
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

void StringArrayDatasetReader::free_cache()
{  
  for(auto& x: cache_) { 
    auto count = x.len;
    char** s = reinterpret_cast<char**>(x.p);
    for(unsigned i=0;i<count;++i,++s) {
      ::free(*s); 
    }
    ::free(x.p);
  }
}

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

BytesDatasetReader::
BytesDatasetReader(const HDFStreamReaderBase* base_ptr, const std::string dataset_name):
  dataset_name_(dataset_name)
{
  // Create opaque 1-byte datatype
  hid_t byte_datatype = H5Tcreate(H5T_OPAQUE, 1);

  // Create variable-length datatype
  byte_string_datatype_ = H5Tvlen_create(byte_datatype);
  H5Tclose(byte_datatype);

  dataset_id_ = H5Dopen(base_ptr->gid(), dataset_name.c_str(), H5P_DEFAULT);
  if (dataset_id_ >= 0) {
    hid_t dataspace_id = H5Dget_space(dataset_id_);
    H5Sget_simple_extent_dims(dataspace_id, &nrow_, NULL);
    H5Sclose(dataspace_id);
  }
}

BytesDatasetReader::~BytesDatasetReader() {
  if(dataset_id_ > 0) {
    H5Dclose(dataset_id_);
  }
  H5Tclose(byte_string_datatype_);
  free_cache();
}

bool BytesDatasetReader::read(uint64_t irow, std::string* s)
{
  if(irow<cache_start_ or irow>=cache_end_) { *s = std::string(); return false; }
  hvl_t x = cache_[irow-cache_start_];
  s->assign(reinterpret_cast<const char*>(x.p), x.len);
  return true;
}

bool BytesDatasetReader::preload(uint64_t start, uint64_t count)
{
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
  herr_t res = H5Dread(dataset_id_, byte_string_datatype_, mem_space_id, file_space_id, H5P_DEFAULT, cache_.data());

  H5Sclose(mem_space_id);
  H5Sclose(file_space_id);

  return res>=0;
}

void BytesDatasetReader::free_cache()
{
  for(auto& x: cache_) { ::free(x.p); }
}


// 888888b.            888                            d8888                                  
// 888  "88b           888                           d88888                                  
// 888  .88P           888                          d88P888                                  
// 8888888K.  888  888 888888 .d88b.  .d8888b      d88P 888 888d888 888d888 8888b.  888  888 
// 888  "Y88b 888  888 888   d8P  Y8b 88K         d88P  888 888P"   888P"      "88b 888  888 
// 888    888 888  888 888   88888888 "Y8888b.   d88P   888 888     888    .d888888 888  888 
// 888   d88P Y88b 888 Y88b. Y8b.          X88  d8888888888 888     888    888  888 Y88b 888 
// 8888888P"   "Y88888  "Y888 "Y8888   88888P' d88P     888 888     888    "Y888888  "Y88888 
//                 888                                                                   888 
//            Y8b d88P                                                              Y8b d88P 
//             "Y88P"                                                                "Y88P"  

BytesArrayDatasetReader::
BytesArrayDatasetReader(const HDFStreamReaderBase* base_ptr, const std::string dataset_name):
  dataset_name_(dataset_name)
{
  // Create opaque 1-byte datatype
  hid_t byte_datatype = H5Tcreate(H5T_OPAQUE, 1);

  // Create variable-length datatype (bytes string)
  hid_t byte_string_datatype = H5Tvlen_create(byte_datatype);
  H5Tclose(byte_datatype);

  // Create variable-length datatype (array of bytes strings)
  array_datatype_ = H5Tvlen_create(byte_string_datatype);
  H5Tclose(byte_string_datatype);

  dataset_id_ = H5Dopen(base_ptr->gid(), dataset_name.c_str(), H5P_DEFAULT);
  if (dataset_id_ >= 0) {
    hid_t dataspace_id = H5Dget_space(dataset_id_);
    H5Sget_simple_extent_dims(dataspace_id, &nrow_, NULL);
    H5Sclose(dataspace_id);
  }
}

BytesArrayDatasetReader::~BytesArrayDatasetReader() {
  if(dataset_id_ > 0) {
    H5Dclose(dataset_id_);
  }
  H5Tclose(array_datatype_);
  free_cache();
}

bool BytesArrayDatasetReader::preload(uint64_t start, uint64_t count)
{
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

void BytesArrayDatasetReader::free_cache()
{  
  for(auto& x: cache_) { 
    auto count = x.len;
    hvl_t* y = reinterpret_cast<hvl_t*>(x.p);
    for(unsigned i=0;i<count;++i,++y) {
      ::free(y->p); 
    }
    ::free(x.p);
  }
}
