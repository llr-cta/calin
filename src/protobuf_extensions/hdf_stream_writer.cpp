/*

  calin/src/protobuf_extensions/hdf_stream_writer.cpp -- Stephen Fegan -- 2025-02-04

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

#include <string>

#include <protobuf_extensions/hdf_stream_writer.hpp>

using namespace calin::protobuf_extensions::hdf_streamer;

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

StringDatasetWriter::
StringDatasetWriter(hid_t gid, const std::string dataset_name, hsize_t nfill):
  dataset_name_(dataset_name)
{
  // Create variable-length datatype
  string_datatype_ = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_datatype_, H5T_VARIABLE); // Variable-length strings

  dataset_id_ = H5Dopen(gid, dataset_name.c_str(), H5P_DEFAULT);
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
    dataset_id_ = H5Dcreate2(gid, dataset_name.c_str(), string_datatype_, 
        dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
                            
    H5Pclose(plist_id);
    H5Sclose(dataspace_id);
  }

  hid_t dataspace_id = H5Dget_space(dataset_id_);
  H5Sget_simple_extent_dims(dataspace_id, &nrow_, NULL);
  H5Sclose(dataspace_id);

  while(nrow_ < nfill) {
    write("");
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

StringDatasetWriter::~StringDatasetWriter() 
{
  flush();
  H5Dclose(dataset_id_);
  H5Tclose(string_datatype_);
}

void StringDatasetWriter::flush() 
{
  hsize_t start;
  hid_t space_id = H5Dget_space(dataset_id_);
  H5Sget_simple_extent_dims(space_id, &start, NULL);
  H5Sclose(space_id);

  hsize_t count = cache_.size();
  hsize_t after = start + count;
  if(after != nrow_) {
    throw std::runtime_error("Unexpected dataset size : " + std::to_string(after) + " != " + std::to_string(nrow_));
  }
  H5Dset_extent(dataset_id_, &after);

  hid_t file_space_id = H5Dget_space(dataset_id_);
  H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, &start, NULL, &count, NULL);
  hid_t mem_space_id = H5Screate_simple(1, &count, NULL);
  H5Dwrite(dataset_id_, string_datatype_, mem_space_id, file_space_id, H5P_DEFAULT, cache_.data());
  
  H5Sclose(mem_space_id);
  H5Sclose(file_space_id);

  for(auto& x : cache_)free(x);
  cache_.clear();
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
StringArrayDatasetWriter::
StringArrayDatasetWriter(hid_t gid, const std::string dataset_name, hsize_t nfill):
  dataset_name_(dataset_name)
{
  // Create string datatype
  hid_t string_datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_datatype, H5T_VARIABLE); // Variable-length strings

  // Create variable-length datatype
  array_datatype_ = H5Tvlen_create(string_datatype);
  H5Tclose(string_datatype);

  dataset_id_ = H5Dopen(gid, dataset_name.c_str(), H5P_DEFAULT);
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
    dataset_id_ = H5Dcreate2(gid, dataset_name.c_str(), array_datatype_, 
        dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
                            
    H5Pclose(plist_id);
    H5Sclose(dataspace_id);
  }

  hid_t dataspace_id = H5Dget_space(dataset_id_);
  H5Sget_simple_extent_dims(dataspace_id, &nrow_, NULL);
  H5Sclose(dataspace_id);

  std::vector<std::string> fill_value;
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

StringArrayDatasetWriter::~StringArrayDatasetWriter() {
  flush();
  H5Dclose(dataset_id_);
  H5Tclose(array_datatype_);
}

void StringArrayDatasetWriter::flush() 
{
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
    throw std::runtime_error("Unexpected dataset size : " + std::to_string(after) + " != " + std::to_string(nrow_));
  }
  H5Dset_extent(dataset_id_, &after);

  hid_t file_space_id = H5Dget_space(dataset_id_);
  H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, &start, NULL, &count, NULL);
  hid_t mem_space_id = H5Screate_simple(1, &count, NULL);
  H5Dwrite(dataset_id_, array_datatype_, mem_space_id, file_space_id, H5P_DEFAULT, hvl.data());
  
  H5Sclose(mem_space_id);
  H5Sclose(file_space_id);

  for(auto& v : cache_) {
    for(auto& x : v) {
      free(x);
    }
  }
  cache_.clear();
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

BytesDatasetWriter::
BytesDatasetWriter(hid_t gid, const std::string dataset_name, hsize_t nfill):
  dataset_name_(dataset_name)
{
  // Create opaque 1-byte datatype
  hid_t byte_datatype = H5Tcreate(H5T_OPAQUE, 1);

  // Create variable-length datatype
  array_datatype_ = H5Tvlen_create(byte_datatype);
  H5Tclose(byte_datatype);

  dataset_id_ = H5Dopen(gid, dataset_name.c_str(), H5P_DEFAULT);
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
    dataset_id_ = H5Dcreate2(gid, dataset_name.c_str(), array_datatype_, 
        dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
                            
    H5Pclose(plist_id);
    H5Sclose(dataspace_id);
  }

  hid_t dataspace_id = H5Dget_space(dataset_id_);
  H5Sget_simple_extent_dims(dataspace_id, &nrow_, NULL);
  H5Sclose(dataspace_id);

  std::string fill_value;
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

BytesDatasetWriter::~BytesDatasetWriter() {
  flush();
  H5Dclose(dataset_id_);
  H5Tclose(array_datatype_);
}

void BytesDatasetWriter::flush() 
{
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
    throw std::runtime_error("Unexpected dataset size : " + std::to_string(after) + " != " + std::to_string(nrow_));
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

BytesArrayDatasetWriter::
BytesArrayDatasetWriter(hid_t gid, const std::string dataset_name, hsize_t nfill):
  dataset_name_(dataset_name)
{
  // Create opaque 1-byte datatype
  hid_t byte_datatype = H5Tcreate(H5T_OPAQUE, 1);

  // Create variable-length datatype (bytes string)
  hid_t bytes_string_datatype = H5Tvlen_create(byte_datatype);
  H5Tclose(byte_datatype);

  // Create variable-length datatype (array of bytes strings)
  array_datatype_ = H5Tvlen_create(bytes_string_datatype);
  H5Tclose(bytes_string_datatype);

  dataset_id_ = H5Dopen(gid, dataset_name.c_str(), H5P_DEFAULT);
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

    H5Tclose(bytes_string_datatype);
    H5Tclose(byte_datatype);

    // Create the dataset with chunking enabled
    dataset_id_ = H5Dcreate2(gid, dataset_name.c_str(), array_datatype_, 
        dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
                            
    H5Pclose(plist_id);
    H5Sclose(dataspace_id);
  }

  hid_t dataspace_id = H5Dget_space(dataset_id_);
  H5Sget_simple_extent_dims(dataspace_id, &nrow_, NULL);
  H5Sclose(dataspace_id);

  std::vector<std::string> fill_value;
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

BytesArrayDatasetWriter::~BytesArrayDatasetWriter() {
  flush();
  H5Dclose(dataset_id_);
  H5Tclose(array_datatype_);
}

void BytesArrayDatasetWriter::flush() 
{
  std::vector<std::vector<hvl_t> > hvl_array;
  hvl_array.reserve(cache_.size());
  for(auto& v : cache_) { 
    std::vector<hvl_t> hvl_v;
    for(auto& x : v) { 
      hvl_v.emplace_back(hvl_t{x.size(), x.data()});
    }
    hvl_array.emplace_back(std::move(hvl_v));
  }

  std::vector<hvl_t> hvl;
  for(auto& x : hvl_array) { 
    hvl.emplace_back(hvl_t{x.size(), x.data()});
  }

  hsize_t start;
  hid_t space_id = H5Dget_space(dataset_id_);
  H5Sget_simple_extent_dims(space_id, &start, NULL);
  H5Sclose(space_id);

  hsize_t count = hvl.size();
  hsize_t after = start + count;
  if(after != nrow_) {
    throw std::runtime_error("Unexpected dataset size : " + std::to_string(after) + " != " + std::to_string(nrow_));
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

// 888888b.                              
// 888  "88b                             
// 888  .88P                             
// 8888888K.   8888b.  .d8888b   .d88b.  
// 888  "Y88b     "88b 88K      d8P  Y8b 
// 888    888 .d888888 "Y8888b. 88888888 
// 888   d88P 888  888      X88 Y8b.     
// 8888888P"  "Y888888  88888P'  "Y8888  

HDFStreamWriterBase::
HDFStreamWriterBase(const std::string& filename, const std::string& groupname, bool truncate, uint64_t cache_size):
  cache_size_(cache_size)
{
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);
  
  if(!truncate)
    h5f_ = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (h5f_ < 0)
    h5f_ = H5Fcreate(filename.c_str(), truncate?H5F_ACC_TRUNC:H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
  if (h5f_ < 0)
    throw std::runtime_error("Failed to create file: " + filename);

  open_group(h5f_, groupname);
}

HDFStreamWriterBase::
HDFStreamWriterBase(hid_t gid, const std::string& groupname, uint64_t cache_size):
  cache_size_(cache_size)
{
  open_group(gid, groupname);
}

HDFStreamWriterBase::
~HDFStreamWriterBase()
{
  flush();
  H5Dclose(h5d_nrow_);
  H5Gclose(h5g_);
  if(h5f_ > 0) {
    H5Fclose(h5f_);
  }
}

void HDFStreamWriterBase::flush()
{
  H5Awrite(h5d_nrow_, H5T_NATIVE_UINT64, &nrow_);
  if(h5f_>0) {
    H5Fflush(h5f_, H5F_SCOPE_LOCAL);
  }
  ncached_ = 0;
}

void HDFStreamWriterBase::
open_group(hid_t file_id, const std::string& groupname)
{
  h5g_ = H5Gopen(file_id, groupname.c_str(), H5P_DEFAULT);
  if (h5g_ < 0) {
    hid_t lcpl_id = H5Pcreate(H5P_LINK_CREATE);
    if (lcpl_id < 0) {
      throw std::runtime_error("Failed to creare link-creation property");
    }
    H5Pset_create_intermediate_group(lcpl_id, 1);
    h5g_ = H5Gcreate(file_id, groupname.c_str(), lcpl_id, H5P_DEFAULT, H5P_DEFAULT);
  }
  if (h5g_ < 0) {
    throw std::runtime_error("Failed to create group: " + groupname);
  }

  h5d_nrow_ = H5Aopen(h5g_, "nrow", H5P_DEFAULT);
  if (h5d_nrow_ < 0) {
    hid_t space_id = H5Screate(H5S_SCALAR);
    h5d_nrow_ = H5Acreate(h5g_, "nrow", H5T_INTEL_U64, space_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space_id);
    if (h5d_nrow_ < 0) {
      throw std::runtime_error("Failed to create attribute: nrow");
    }
    if(H5Awrite(h5d_nrow_, H5T_NATIVE_UINT64, &nrow_) < 0) {
      throw std::runtime_error("Failed to write attribute: nrow");
    }
  } else {
    if(H5Aread(h5d_nrow_, H5T_NATIVE_UINT64, &nrow_) < 0) {
      throw std::runtime_error("Failed to read attribute: nrow");
    }
  }
}
