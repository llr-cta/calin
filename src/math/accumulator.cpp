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
