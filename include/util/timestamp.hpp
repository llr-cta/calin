/*

   calin/util/timestamp.hpp -- Stephen Fegan -- 2017-10-29

   Simple (non-precision) timestamp class

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <sys/time.h>

#include <util/timestamp.pb.h>

namespace calin { namespace util { namespace timestamp {

class Timestamp
{
public:
  Timestamp(): nsec_() { }
  Timestamp(int64_t tval, int64_t tmultiplier = 1LL): nsec_(tval * tmultiplier) { }
  Timestamp(const timespec &tp): nsec_(tp.tv_sec * 1000000000LL + tp.tv_nsec) { }
  Timestamp(const timeval &tp): nsec_(tp.tv_sec * 1000000000LL + tp.tv_usec*1000LL) { }
  uint64_t unix_sec() const { return nsec_/1000000000LL; }
  uint32_t unix_nsec() const { return nsec_%1000000000LL; }
  std::string as_string(bool utc = false) const;
#ifndef SWIG
  calin::ix::util::timestamp::Timestamp* as_proto(calin::ix::util::timestamp::Timestamp* x = nullptr, bool utc = false) const;
  calin::ix::util::timestamp::Timestamp* as_proto(bool utc, calin::ix::util::timestamp::Timestamp* x = nullptr) const;
#else
  calin::ix::util::timestamp::Timestamp* as_proto(bool utc = false) const;
  void as_proto(calin::ix::util::timestamp::Timestamp* x, bool utc = false) const;
  void as_proto(bool utc, calin::ix::util::timestamp::Timestamp* x) const;
#endif
  static Timestamp now();
  double seconds_since(const Timestamp& then) const;
  int64_t nanoseconds_since(const Timestamp& then) const;
private:
  int64_t nsec_ = 0;
};

} } } // namespace calin::util::timestamp
