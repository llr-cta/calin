/*

   calin/math/accumulator.hpp -- Stephen Fegan -- 2015-02-10

   Classes to sum up double precision values in a reliable way

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <limits>
#include <vector>
#include <array>

#include "calin_global_definitions.hpp"

namespace calin { namespace math { namespace accumulator {

CALIN_TYPEALIAS(accumulator_value_type, double);

class SimpleAccumulator
{
public:
  using value_type = accumulator_value_type;
  SimpleAccumulator(const value_type val=0): sum_{val} { /* nothing to see here */ }
  void reset(double val=0) { sum_ = val; }
  void accumulate(const value_type x) { sum_ += x; }
  value_type total() const { return sum_; };
  operator value_type() const { return total(); }
  bool operator== (const SimpleAccumulator& o) const { return sum_==o.sum_; }
private:
  value_type sum_;
};

template<typename value_type> class BasicKahanAccumulator
{
  // Kahan summation algorithm - summation with correction
  // Reference: http://en.wikipedia.org/wiki/Kahan_summation_algorithm
public:
  // using value_type = T;

  BasicKahanAccumulator(const value_type val=0): sum_{val}, cor_{0} { /* nothing to see here */ }
  void reset(double val=0) { sum_ = val; cor_ = 0; }
  void accumulate(const value_type x)
  {
    value_type Y { x-cor_ };
    value_type T { sum_+Y };
    cor_=(T-sum_)-Y;
    sum_=T;
  }
  value_type total() const { return sum_; }
  value_type correction() const { return -cor_; }
  std::vector<value_type> corrections() const { return { correction() }; }
  operator value_type() const { return total(); }
  bool operator== (const BasicKahanAccumulator<value_type>& o) const
  {
    return sum_==o.sum_ and cor_==o.cor_;
  }
private:
  value_type sum_;
  value_type cor_;
};

CALIN_TYPEALIAS(KahanAccumulator, BasicKahanAccumulator<double>);

#if 0

class CascadingIntegerAccumulator
{
 public:
  using value_type = accumulator_value_type;

  CascadingIntegerAccumulator(double val=0) { /* nothing to see here */ }
  void reset(double val=0) { acc_.fill(0); }
  void accumulate(const value_type x);
  value_type total() const;
  std::vector<value_type> corrections() const;
  operator value_type() const { return total(); }

 private:
  using acc_type = int64_t;

  static constexpr int mantissa_digits { std::numeric_limits<value_type>::digits-1 }; // Subtract one as leading "1" is suppressed
  static constexpr int exponent_digits { sizeof(value_type)*8-mantissa_digits-1 };
  static constexpr int num_accumulators { 1<<exponent_digits };

  static constexpr uint64_t mantissa_mask { (uint64_t{1}<<mantissa_digits)-1 };
  static constexpr uint64_t exponent_mask { (uint64_t{1}<<exponent_digits)-1 };

  mutable std::array<acc_type,num_accumulators> acc_ {};
};

#endif

class Accumulator
{
public:
  virtual ~Accumulator() = 0;
  virtual void reset(double val=0);
  virtual void accumulate(const double x) = 0;
  virtual double total() = 0;
};

template<typename Acc> class BasicAccumulator: public Accumulator
{
public:
  BasicAccumulator(double val=0): Accumulator{}, acc_{val}
  { /* nothing to see here */ }
  ~BasicAccumulator() { /* nothing to see here */ }
  virtual void reset(double val=0) override { acc_.reset(val); }
  virtual void accumulate(const double x) override { acc_.accumulate(x); }
  virtual double total() override { return acc_.total(); }
private:
  Acc acc_;
};

CALIN_TYPEALIAS(RecommendedAccumulator, KahanAccumulator);
CALIN_TYPEALIAS(LikelihoodAccumulator, KahanAccumulator);

} } } // namespace calin::math::accumulator
