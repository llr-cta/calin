/*

   calin/util/memory.hpp -- Stephen Fegan -- 2018-01-18

   Some memory utility functions

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>
#include <memory>
#include <cstdlib>
#include <type_traits>

#include <provenance/system_info.hpp>

namespace calin { namespace util { namespace memory {

#if defined(__AVX512F__)
constexpr size_t default_align = 512/8;
#elif defined(__AVX__)
constexpr size_t default_align = 256/8;
#else
constexpr size_t default_align = 128/8;
#endif

template<typename T> inline T* aligned_calloc(const size_t num, size_t align = default_align)
{
  align = std::max(std::max(sizeof(T),sizeof(void*)), align);
  void* ptr = nullptr;
  const size_t alloc_size = num*sizeof(T)+align;
  if(::posix_memalign(&ptr, align, alloc_size)==0)
    return reinterpret_cast<T*>(ptr);
  return nullptr;
}

template<typename T> inline bool aligned_calloc(T*& ptr, const size_t num, const size_t align = default_align)
{
  ptr = aligned_calloc<T>(num, align);
  return ptr != nullptr;
}

template<typename T> inline bool aligned_calloc(T*__restrict__& ptr, const size_t num, const size_t align = default_align)
{
  ptr = aligned_calloc<T>(num, align);
  return ptr != nullptr;
}

template<typename T> inline bool aligned_recalloc(T*& ptr, const size_t num, const size_t align = default_align)
{
  ::free(ptr);
  ptr = aligned_calloc<T>(num, align);
  return ptr != nullptr;
}

template<typename T> inline bool aligned_recalloc(T*__restrict__& ptr, const size_t num, const size_t align = default_align)
{
  ::free(ptr);
  ptr = aligned_calloc<T>(num, align);
  return ptr != nullptr;
}

template<typename T> inline T* safe_aligned_calloc(const size_t num, const size_t align = default_align)
{
  T* ptr = aligned_calloc<T>(num, align);
  if(ptr != nullptr)return ptr;
  throw std::runtime_error("Could not allocate array of " + std::to_string(num) +
    " elements of size " + std::to_string(sizeof(T)) + " bytes");
}

template<typename T> inline void safe_aligned_calloc(T*& ptr, const size_t num, const size_t align = default_align)
{
  ptr = safe_aligned_calloc<T>(num, align);
}

template<typename T> inline void safe_aligned_calloc(T*__restrict__& ptr, const size_t num, const size_t align = default_align)
{
  ptr = safe_aligned_calloc<T>(num, align);
}

template<typename T> inline void safe_aligned_recalloc(T*& ptr, const size_t num, const size_t align = default_align)
{
  free(ptr);
  ptr = safe_aligned_calloc<T>(num, align);
}

template<typename T> inline void safe_aligned_recalloc(T*__restrict__& ptr, const size_t num, const size_t align = default_align)
{
  free(ptr);
  ptr = safe_aligned_calloc<T>(num, align);
}

template<typename T> inline T* safe_aligned_calloc_and_fill(const size_t num, const T& fill_val = T(), const size_t align = default_align)
{
  T* ptr = aligned_calloc<T>(num, align);
  if(ptr != nullptr) {
    std::fill(ptr, ptr+num, fill_val);
    return ptr;
  }
  throw std::runtime_error("Could not allocate array of " + std::to_string(num) +
    " elements of size " + std::to_string(sizeof(T)) + " bytes");
}

template<typename T> inline void safe_aligned_calloc_and_fill(T*& ptr, const size_t num, const T& fill_val = T(), const size_t align = default_align)
{
  ptr = safe_aligned_calloc_and_fill<T>(num, fill_val, align);
}

template<typename T> inline void safe_aligned_calloc_and_fill(T*__restrict__& ptr, const size_t num, const T& fill_val = T(), const size_t align = default_align)
{
  ptr = safe_aligned_calloc_and_fill<T>(num, fill_val, align);
}

template<typename T> inline void safe_aligned_recalloc_and_fill(T*& ptr, const size_t num, const T& fill_val = T(), const size_t align = default_align)
{
  free(ptr);
  ptr = safe_aligned_calloc_and_fill<T>(num, fill_val, align);
}

template<typename T> inline void safe_aligned_recalloc_and_fill(T*__restrict__& ptr, const size_t num, const T& fill_val = T(), const size_t align = default_align)
{
  free(ptr);
  ptr = safe_aligned_calloc_and_fill<T>(num, fill_val, align);
}


} } } // namespace calin::util::memory
