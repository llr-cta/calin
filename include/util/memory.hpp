/*

   calin/util/memory.hpp -- Stephen Fegan -- 2018-01-18

   Some memory utility functions

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

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

namespace calin { namespace util { namespace memory {

template<typename T> inline T* aligned_calloc(const size_t num, const size_t log2_align = 3)
{
  const size_t align = 1<<log2_align;
  assert(align >= sizeof(void*));
  void* ptr = nullptr;
  const size_t alloc_size = ((num*sizeof(T)+align-1)>>log2_align)<<log2_align;
//  std::cout << num << ' ' << sizeof(T) << ' ' << alloc_size << ' ' << align << ' ' << log2_align << '\n';
  if(::posix_memalign(&ptr, align, alloc_size)==0)
    return reinterpret_cast<T*>(ptr);
  return nullptr;
}

template<typename T> inline bool aligned_calloc(T*& ptr, const size_t num, const size_t log2_align = 3)
{
  ptr = aligned_calloc<T>(num, log2_align);
  return ptr != nullptr;
}

template<typename T> inline bool aligned_calloc(T*__restrict__& ptr, const size_t num, const size_t log2_align = 3)
{
  ptr = aligned_calloc<T>(num, log2_align);
  return ptr != nullptr;
}

template<typename T> inline bool aligned_recalloc(T*& ptr, const size_t num, const size_t log2_align = 3)
{
  ::free(ptr);
  ptr = aligned_calloc<T>(num, log2_align);
  return ptr != nullptr;
}

template<typename T> inline bool aligned_recalloc(T*__restrict__& ptr, const size_t num, const size_t log2_align = 3)
{
  ::free(ptr);
  ptr = aligned_calloc<T>(num, log2_align);
  return ptr != nullptr;
}

template<typename T> inline T* safe_aligned_calloc(const size_t num, const size_t log2_align = 3)
{
  T* ptr = aligned_calloc<T>(num, log2_align);
  if(ptr != nullptr)return ptr;
  throw std::runtime_error("Could not allocate array of " + std::to_string(num) +
    " elements of size " + std::to_string(sizeof(T)) + " bytes");
}

template<typename T> inline void safe_aligned_calloc(T*& ptr, const size_t num, const size_t log2_align)
{
  ptr = safe_aligned_calloc<T>(num, log2_align);
}

template<typename T> inline void safe_aligned_calloc(T*__restrict__& ptr, const size_t num, const size_t log2_align)
{
  ptr = safe_aligned_calloc<T>(num, log2_align);
}

template<typename T> inline void safe_aligned_recalloc(T*& ptr, const size_t num, const size_t log2_align)
{
  free(ptr);
  ptr = safe_aligned_calloc<T>(num, log2_align);
}

template<typename T> inline void safe_aligned_recalloc(T*__restrict__& ptr, const size_t num, const size_t log2_align)
{
  free(ptr);
  ptr = safe_aligned_calloc<T>(num, log2_align);
}

} } } // namespace calin::util::memory
