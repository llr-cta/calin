/*

   calin/util/vcl.cpp -- Stephen Fegan -- 2018-08-14

   Vector class library helpers

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

#include <util/vcl.hpp>

using namespace calin::util::vcl;

const char VCL128Architecture::architecture_name[];
const unsigned VCL128Architecture::vec_bits;
const unsigned VCL128Architecture::vec_bytes;
const unsigned VCL128Architecture::num_int8;
const unsigned VCL128Architecture::num_uint8;
const unsigned VCL128Architecture::num_int16;
const unsigned VCL128Architecture::num_uint16;
const unsigned VCL128Architecture::num_int32;
const unsigned VCL128Architecture::num_uint32;
const unsigned VCL128Architecture::num_int64;
const unsigned VCL128Architecture::num_uint64;
const unsigned VCL128Architecture::num_float;
const unsigned VCL128Architecture::num_double;

const char VCL256Architecture::architecture_name[];
const unsigned VCL256Architecture::vec_bits;
const unsigned VCL256Architecture::vec_bytes;
const unsigned VCL256Architecture::num_int8;
const unsigned VCL256Architecture::num_uint8;
const unsigned VCL256Architecture::num_int16;
const unsigned VCL256Architecture::num_uint16;
const unsigned VCL256Architecture::num_int32;
const unsigned VCL256Architecture::num_uint32;
const unsigned VCL256Architecture::num_int64;
const unsigned VCL256Architecture::num_uint64;
const unsigned VCL256Architecture::num_float;
const unsigned VCL256Architecture::num_double;

const char VCL512Architecture::architecture_name[];
const unsigned VCL512Architecture::vec_bits;
const unsigned VCL512Architecture::vec_bytes;
const unsigned VCL512Architecture::num_int8;
const unsigned VCL512Architecture::num_uint8;
const unsigned VCL512Architecture::num_int16;
const unsigned VCL512Architecture::num_uint16;
const unsigned VCL512Architecture::num_int32;
const unsigned VCL512Architecture::num_uint32;
const unsigned VCL512Architecture::num_int64;
const unsigned VCL512Architecture::num_uint64;
const unsigned VCL512Architecture::num_float;
const unsigned VCL512Architecture::num_double;
