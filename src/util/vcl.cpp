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

constexpr char VCL128Architecture::architecture_name[];
constexpr unsigned VCL128Architecture::vec_bits;
constexpr unsigned VCL128Architecture::vec_bytes;
constexpr unsigned VCL128Architecture::num_int8;
constexpr unsigned VCL128Architecture::num_uint8;
constexpr unsigned VCL128Architecture::num_int16;
constexpr unsigned VCL128Architecture::num_uint16;
constexpr unsigned VCL128Architecture::num_int32;
constexpr unsigned VCL128Architecture::num_uint32;
constexpr unsigned VCL128Architecture::num_int64;
constexpr unsigned VCL128Architecture::num_uint64;
constexpr unsigned VCL128Architecture::num_float;
constexpr unsigned VCL128Architecture::num_double;

constexpr char VCL256Architecture::architecture_name[];
constexpr unsigned VCL256Architecture::vec_bits;
constexpr unsigned VCL256Architecture::vec_bytes;
constexpr unsigned VCL256Architecture::num_int8;
constexpr unsigned VCL256Architecture::num_uint8;
constexpr unsigned VCL256Architecture::num_int16;
constexpr unsigned VCL256Architecture::num_uint16;
constexpr unsigned VCL256Architecture::num_int32;
constexpr unsigned VCL256Architecture::num_uint32;
constexpr unsigned VCL256Architecture::num_int64;
constexpr unsigned VCL256Architecture::num_uint64;
constexpr unsigned VCL256Architecture::num_float;
constexpr unsigned VCL256Architecture::num_double;

constexpr char VCL512Architecture::architecture_name[];
constexpr unsigned VCL512Architecture::vec_bits;
constexpr unsigned VCL512Architecture::vec_bytes;
constexpr unsigned VCL512Architecture::num_int8;
constexpr unsigned VCL512Architecture::num_uint8;
constexpr unsigned VCL512Architecture::num_int16;
constexpr unsigned VCL512Architecture::num_uint16;
constexpr unsigned VCL512Architecture::num_int32;
constexpr unsigned VCL512Architecture::num_uint32;
constexpr unsigned VCL512Architecture::num_int64;
constexpr unsigned VCL512Architecture::num_uint64;
constexpr unsigned VCL512Architecture::num_float;
constexpr unsigned VCL512Architecture::num_double;
