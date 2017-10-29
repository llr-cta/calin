/*

   calin/util/timestamp.hpp -- Stephen Fegan -- 2017-10-29

   Simple (non-precision) timestamp class

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <util/timestamp.pb.h>

namespace calin { namespace util { namespace timestamp {

class Timestamp
{
public:
  Timestamp(uint64_t sec = 0, uint32_t usec = 0): sec_(sec), usec_(usec) { }
  uint64_t sec() const { return sec_; }
  uint32_t usec() const { return usec_; }
  std::string as_string() const;
#ifndef SWIG
  calin::ix::util::timestamp::Timestamp* as_proto(calin::ix::util::timestamp::Timestamp* x = nullptr) const;
#else
  calin::ix::util::timestamp::Timestamp* as_proto() const;
  void as_proto(calin::ix::util::timestamp::Timestamp* x = nullptr) const;
#endif
  static Timestamp now();
  double seconds_since(const Timestamp& then) const;
private:
  uint64_t sec_ = 0;
  uint32_t usec_ = 0;
};

} } } // namespace calin::util::timestamp
