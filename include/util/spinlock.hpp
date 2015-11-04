/* 

   calin/util/spinlock.hpp -- Stephen Fegan -- 2015-05-15

   Spinlock implementation cobbled together from various places on the web

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

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

namespace calin { namespace util { namespace parallel {

#if defined(__i386__) || defined(__x86_64__) || \
    defined(__amd64__) || defined(__ia64__)
#define SPINWAIT()  do { __asm __volatile("pause") } while (0)
#else
#define SPINWAIT()  do { /* nothing */ } while (0)
#endif

class spinlock
{
 public:
  inline void lock() {
    while (locked.test_and_set(std::memory_order_acquire)) { SPINWAIT(); }
  }
  
  inline void unlock() {
    locked.clear(std::memory_order_release);
  }

 private:
  std::atomic_flag locked = ATOMIC_FLAG_INIT ;
};

#undef SPINWAIT()

} } } // namespace calin::util::thread
