/*

   calin/diagnostics/range.hpp -- Stephen Fegan -- 2018-10-28

   Various helper functions for ranges

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

#include <diagnostics/range.pb.h>

namespace calin { namespace diagnostics { namespace range {

template<typename RunLengthEncodingType, typename ValueType>
inline void encode_value(RunLengthEncodingType* rle, ValueType value)
{
  assert(rle->value_size() == rle->count_size());
  if(rle->value_size()==0 or rle->value(rle->value_size()-1) != value) {
    rle->add_count(1);
    rle->add_value(value);
  } else {
    unsigned icount = rle->count_size()-1;
    rle->set_count(icount, rle->count(icount)+1);
  }
}

} } } // namespace calin::diagnostics::range
