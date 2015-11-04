/* 

   calin/math/accumulator.cpp -- Stephen Fegan -- 2015-02-10

   Classes to sum up double precision values in a reliable way

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

#include "math/accumulator.hpp"

#include <iostream>

using namespace calin::math::accumulator;

#if 0
void CascadingIntegerAccumulator::accumulate(const value_type x)
{
  uint64_t uix { *reinterpret_cast<const uint64_t*>(&x) };
  uint64_t mantissa { uix&mantissa_mask };
  uint64_t exponent { (uix>>mantissa_digits) & exponent_mask };

  //acc_[exponent] += x;
  std::cout << std::dec << mantissa_digits << ' ' << std::hex
            << mantissa_mask << ' ' << exponent_mask << ' '
            << uix << ' '
            << mantissa << ' ' << exponent << '\n';
}

auto CascadingIntegerAccumulator::total() const -> value_type
{
  return 0;
}

#endif
