/* 

   calin/util/spinlock.hpp -- Stephen Fegan -- 2015-05-15

   Spinlock implementation cobbled together from various places on the web

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
