/*

   calin/math/histogram.hpp -- Stephen Fegan -- 2015-02-16

   Simple histogramming classes

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

#include <string>
#include <algorithm>
#include <vector>
#include <deque>
#include <cmath>
#include <limits>
#include <iterator>

#include <Eigen/Core>

#include "calin_global_definitions.hpp"
#include "math/histogram.pb.h"
#include "math/accumulator.hpp"

// This class has some changes/improvements over the VERITAS ChiLA
// implementation.  They are motivated by our preference for
// simplicity over "speed at all costs", and desire for mathematical
// rigor.

// Major changes are:
// - Remove temlated type for binned variable - now only double is supported
// - Use deque instead of two vectors for simplicity (and to remove
//   inconsistency in rounding at the bin edges). If this turns out to not
//   be fast enough a custom double-ended vector would be the solution.
// - Support diffrent accumulators for bin weights. SimpleAccumulator is
//   sufficient for most cases, but KahanAccumulator may be used if necessary
// - Include limits as option in this class rather than in separate LimitedHist

// Design has regrettibly evolved from one simple class to include a
// multitude of base and support classes. It is in danger of becoming
// overly complexified

// 2016-04-18 : yes, I agree with my earlier self, it is much too complex

namespace calin { namespace math { namespace histogram {

// Base class for one-dimensional binned data -- separated from the
// actual historam class since it is in common with the integral
// distribution (CDF) class.

CALIN_TYPEALIAS(DefaultAccumulator,
                calin::math::accumulator::SimpleAccumulator);

template<typename T, typename Container = std::vector<T>> class BinnedData1D
{
public:
  CALIN_TYPEALIAS(data_type, T);
  CALIN_TYPEALIAS(data_container_type, Container);

  BinnedData1D(double dxval, double xval_align = 0.5,
               const std::string& xval_units = std::string()):
      dxval_{dxval}, dxval_inv_{1.0/dxval_}, xval_align_{xval_align}, xval_units_{xval_units}
    { /* nothing to see here */ }
  BinnedData1D(double dxval, double xval_limit_lo, double xval_limit_hi,
               double xval_align = 0.5,
               const std::string& xval_units = std::string()):
      dxval_{dxval}, dxval_inv_{1.0/dxval_}, xval_align_{xval_align}, limited_{true},
      xval_limit_lo_{xval_limit_lo}, xval_limit_hi_{xval_limit_hi},
      xval_units_{xval_units}
    { /* nothing to see here */ }
  BinnedData1D(const calin::ix::math::histogram::Histogram1DConfig& config):
      dxval_{config.dxval()}, dxval_inv_{1.0/dxval_}, xval_align_{config.xval_align()},
      limited_{config.limited()}, xval_limit_lo_{config.xval_limit_lo()},
      xval_limit_hi_{config.xval_limit_hi()}, xval_units_{config.xval_units()}
    { /* nothing to see here */ }
  BinnedData1D(const calin::ix::math::histogram::AccumulatedAndSerializedHistogram1DConfig& config):
      dxval_{config.dxval()}, dxval_inv_{1.0/dxval_}, xval_align_{config.xval_align()},
      limited_{config.limited()}, xval_limit_lo_{config.xval_limit_lo()},
      xval_limit_hi_{config.xval_limit_hi()}, xval_units_{config.xval_units()}
    { /* nothing to see here */ }

  BinnedData1D(const BinnedData1D& o) = default;
#ifndef SWIG
  BinnedData1D& operator=(const BinnedData1D& o) = default;
#endif

  // Getters and setters
  double dxval() const { return dxval_; }
  double xval_align() const { return xval_align_; }
  std::string xval_units() const { return xval_units_; }
  void set_xval_units(const std::string& units) { xval_units_=units; }
  double xval0() const { return xval0_; }

  // Functions to get size of histogram and clear it
  bool empty() const { return bins_.empty(); }
  int size() const { return static_cast<int>(bins_.size()); }
  int nbin() const { return size(); }
  void clear() { bins_.clear(); overflow_lo_=T{}; overflow_hi_=T{}; }

  // Bin to value and vice-versa
  double xval_left(int ibin) const { return ibin*dxval_+xval0_; }
  double xval_right(int ibin) const { return (ibin+1)*dxval_+xval0_; }
  double xval_center(int ibin) const { return (ibin+0.5)*dxval_+xval0_; }
  Eigen::VectorXd all_xval_left() const { Eigen::VectorXd x(bins_.size());
    for(int ibin=0;ibin<this->size();ibin++)x[ibin] = xval_left(ibin);
    return x; }
  Eigen::VectorXd all_xval_right() const { Eigen::VectorXd x(bins_.size());
    for(int ibin=0;ibin<this->size();ibin++)x[ibin] = xval_right(ibin);
    return x; }
  Eigen::VectorXd all_xval_center() const { Eigen::VectorXd x(bins_.size());
    for(int ibin=0;ibin<this->size();ibin++)x[ibin] = xval_center(ibin);
    return x; }
  int ibin(double x) const { return std::floor((x-xval0_)*dxval_inv_); }
  int ibin_and_rem(double x, double& dx) const { int ix = std::floor((x-xval0_)*dxval_inv_); dx = x-xval_left(ix); return ix; }
  bool has_ibin(int ibin) const {
    return ibin>=0 and (unsigned)ibin<bins_.size(); }
  bool has_xval(double x) const
  {
    if (x<xval0_)return false;
    return (unsigned)ibin(x) < bins_.size();
  }
  static double xalign_to_center_value(double xval, double dx)
  {
    return 0.5-xval+std::round(xval/dx)*dx;
  }

  // Limits
  bool is_limited() const { return limited_; }
  double xval_limit_lo() const { return xval_limit_lo_; }
  double xval_limit_hi() const { return xval_limit_hi_; }

  void trim_zeros() {
    if(bins_.empty())return;
    unsigned istart = 0;
    for(auto i=bins_.begin(); i!=bins_.end(); ++i, ++istart)
      if(*i != T())break;
    unsigned iend = bins_.size();
    for(auto i=bins_.end(); i!=bins_.begin(); --iend)
      if(*(--i) != T())break;
    xval0_ = xval_left(istart);
    Container c;
    while(istart != iend)c.emplace_back(bins_[istart++]);
    bins_ = std::move(c);
  }

protected:
  BinnedData1D(double dxval, double xval_align, double xval0,
               bool limited, double xval_limit_lo, double xval_limit_hi,
               const std::string& xval_units):
      dxval_(dxval), dxval_inv_(1.0/dxval_), xval_align_(xval_align), xval0_(xval0),
      xval_limit_lo_(xval_limit_lo), xval_limit_hi_(xval_limit_hi),
      xval_units_(xval_units) { /* nothing to see here */ }

#ifndef SWIG
  template<typename iterator>
  BinnedData1D(double dxval, double xval_align, double xval0,
               iterator bins_data_begin, iterator bins_data_end,
               bool limited, double xval_limit_lo, double xval_limit_hi,
               const T& overflow_lo, const T& overflow_hi,
               const std::string& xval_units):
      dxval_(dxval), dxval_inv_(1.0/dxval_), xval_align_(xval_align), xval0_(xval0),
      bins_(bins_data_begin, bins_data_end), limited_(limited),
      xval_limit_lo_(xval_limit_lo), xval_limit_hi_(xval_limit_hi),
      overflow_lo_(overflow_lo), overflow_hi_(overflow_hi),
      xval_units_(xval_units) { /* nothing to see here */ }
#endif

  // Retrieve value for bin
  T& bin(int ibin) { return bins_[ibin]; }
  T& checked_bin(int ibin) { return bins_.at(ibin); }
  T& overflow_lo() { return overflow_lo_; }
  T& overflow_hi() { return overflow_hi_; }

  const T& bin(int ibin) const { return bins_[ibin]; }
  const T& checked_bin(int ibin) const { return bins_.at(ibin); }
  const T& overflow_lo() const { return overflow_lo_; }
  const T& overflow_hi() const { return overflow_hi_; }

  T& bin_with_extend(const double x, bool& x_outside_limits)
  {
    if(!std::isfinite(x))throw std::out_of_range("bin_with_extend");
    if(limited_)
    {
      if(x < xval_limit_lo_) { x_outside_limits = true; return overflow_lo_; }
      if(x >= xval_limit_hi_) { x_outside_limits = true; return overflow_hi_; }
    }
    x_outside_limits = false;
    return unchecked_bin_with_extend(x);
  }

  T& unchecked_bin_with_extend(const double x)
  {
    int thebin = bins_.empty() ? -1 : ibin(x);
    if(thebin < 0)
    {
      do {
        bins_.emplace_front();
        ++thebin;
      }while(thebin<0);
      //xval0_ = std::floor((x+xval_align_)*dxval_inv_)*dxval_ - xval_align_;
      xval0_ = (std::floor(x*dxval_inv_+xval_align_) - xval_align_)*dxval_;
    }
    else if(static_cast<unsigned>(thebin) >= bins_.size())
    {
      do {
        bins_.emplace_back();
      }while(static_cast<unsigned>(thebin) >= bins_.size());
    }
    return bins_[thebin];
  }

  T& bin_with_extend_back(const double x, bool& x_outside_limits)
  {
    if(!std::isfinite(x))throw std::out_of_range("bin_with_extend_back");
    if(limited_)
    {
      if(x < xval_limit_lo_) { x_outside_limits = true; return overflow_lo_; }
      if(x >= xval_limit_hi_) { x_outside_limits = true; return overflow_hi_; }
    }
    x_outside_limits = false;
    return unchecked_bin_with_extend_back(x);
  }

  T& unchecked_bin_with_extend_back(const double x)
  {
    if(bins_.empty())
    {
      bins_.emplace_back();
      //xval0_ = std::floor((x+xval_align_)*dxval_inv_)*dxval_ - xval_align_;
      xval0_ = (std::floor(x*dxval_inv_+xval_align_) - xval_align_)*dxval_;
      return bins_.front();
    }
    int thebin { ibin(x) };
    while(thebin<=bins_.size())bins_.emplace_back();
    return bins_[thebin];
  }

  double dxval_           = 1.0;
  double dxval_inv_       = 1.0;
  double xval_align_      = 0.5;
  double xval0_           = 0;
  Container bins_;
  bool limited_           = false;
  double xval_limit_lo_   = -std::numeric_limits<double>::infinity();
  double xval_limit_hi_   = std::numeric_limits<double>::infinity();
  T overflow_lo_          = {};
  T overflow_hi_          = {};
  std::string xval_units_;
};

#ifdef SWIG
// What a mess!
%template (BinnedDataCDFBase) BinnedData1D<double>;
%template (BinnedDataHistogramBase) BinnedData1D<calin::math::accumulator::SimpleAccumulator,std::deque<calin::math::accumulator::SimpleAccumulator>>;
#endif

// ============================================================================
//
// Bin accessor and iterator -- the accessor provides default access to the
// x-value of the bin and its bin number. Sub-classes should provide access to
// the data itself in whatever way is meaningful
//
// ============================================================================

#ifndef SWIG

template<typename DataBinner> class basic_bin_accessor
{
 public:
  CALIN_TYPEALIAS(data_binner_type, DataBinner);
  CALIN_TYPEALIAS(data_type, typename DataBinner::data_type);

  basic_bin_accessor(DataBinner& binner, int ibin):
      binner_{&binner},ibin_{ibin} {}

  basic_bin_accessor(const basic_bin_accessor&) = default;
  basic_bin_accessor& operator=(const basic_bin_accessor&) = default;

  DataBinner* binner() { return binner_; }
  int ibin() { return ibin_; }
  double dxval() const { return binner_->dxval(); }
  double xval_left() const { return binner_->xval_left(ibin_); }
  double xval_right() const { return binner_->xval_right(ibin_); }
  double xval_center() const { return binner_->xval_center(ibin_); }
  int ibin() const { return ibin_; }
  bool has_bin() const { return binner_->has_ibin(ibin_); }
 protected:
  const data_type& data() const { return binner_->bin(ibin_); }
  data_type& data() { return binner_->bin(ibin_); }
  DataBinner* binner_;
  int ibin_;
};

template<typename DataBinner,
         typename bin_accessor = basic_bin_accessor<DataBinner> >
class basic_iterator:
      public std::iterator<std::bidirectional_iterator_tag,
                           bin_accessor, int, bin_accessor*, bin_accessor&>,
      protected bin_accessor
{
 public:
  CALIN_TYPEALIAS(bin_accessor_type, bin_accessor);
  CALIN_TYPEALIAS(data_binner_type, DataBinner);
  CALIN_TYPEALIAS(data_type, typename bin_accessor::data_type);

  basic_iterator(DataBinner& data, int ibin):
      bin_accessor {data,ibin} {}

  basic_iterator(const basic_iterator&) = default;
  basic_iterator& operator=(const basic_iterator&) = default;

  bin_accessor* operator->() { return this; }
  bin_accessor& operator*() { return *this; }

  basic_iterator& operator++() {
    ++this->ibin_; return *this; }
  basic_iterator operator++(int) {
    basic_iterator i=*this; ++this->ibin_; return i; }
  basic_iterator& operator--() {
    --this->ibin_; return *this; }
  basic_iterator operator--(int) {
    basic_iterator i=*this; --this->ibin_; return i; }
  basic_iterator& operator+=(int ibin) {
    this->ibin_+=ibin; return *this; }
  basic_iterator& operator-=(int ibin) {
    this->ibin_-=ibin; return *this; }
  basic_iterator operator+(int ibin) const {
    basic_iterator i(*this); i+=ibin; return i; }
  basic_iterator operator-(int ibin) const {
    basic_iterator i(*this); i-=ibin; return i; }
  int operator-(const basic_iterator& o) {
    return this->ibin_ - o.ibin_; }
  bool operator!=(const basic_iterator& o) const {
    return this->ibin_ != o.ibin_; }
  bool operator==(const basic_iterator& o) const {
    return this->ibin_ == o.ibin_; }
  bool operator<(const basic_iterator& o) const {
    return this->ibin_ < o.ibin_; }
  bool operator<=(const basic_iterator& o) const {
    return this->ibin_ <= o.ibin_; }
  bool operator>(const basic_iterator& o) const {
    return this->ibin_ > o.ibin_; }
  bool operator>=(const basic_iterator& o) const {
    return this->ibin_ >= o.ibin_; }
};

template<typename DataBinner,
         typename bin_accessor = basic_bin_accessor<DataBinner> >
class basic_reverse_iterator:
      public std::iterator<std::bidirectional_iterator_tag,
                           bin_accessor, int, bin_accessor*, bin_accessor&>,
      protected bin_accessor
{
 public:
  CALIN_TYPEALIAS(bin_accessor_type, bin_accessor);
  CALIN_TYPEALIAS(data_binner_type, DataBinner);
  CALIN_TYPEALIAS(data_type, typename bin_accessor::data_type);

  basic_reverse_iterator(DataBinner& data, int ibin):
      bin_accessor {data,data.size()-ibin-1} {}

  basic_reverse_iterator(const basic_reverse_iterator&) = default;
  basic_reverse_iterator& operator=(const basic_reverse_iterator&) = default;

  bin_accessor* operator->() { return this; }
  bin_accessor& operator*() { return *this; }

  basic_reverse_iterator& operator++() {
    --this->ibin_; return *this; }
  basic_reverse_iterator operator++(int) {
    basic_reverse_iterator i=*this; --this->ibin_; return i; }
  basic_reverse_iterator& operator--() {
    ++this->ibin_; return *this; }
  basic_reverse_iterator operator--(int) {
    basic_reverse_iterator i=*this; ++this->ibin_; return i; }
  basic_reverse_iterator& operator+=(int ibin) {
    this->ibin_-=ibin; return *this; }
  basic_reverse_iterator& operator-=(int ibin) {
    this->ibin_+=ibin; return *this; }
  basic_reverse_iterator operator+(int ibin) const {
    basic_reverse_iterator i(*this); i+=ibin; return i; }
  basic_reverse_iterator operator-(int ibin) const {
    basic_reverse_iterator i(*this); i-=ibin; return i; }
  int operator-(const basic_reverse_iterator& o) {
    return o.ibin_ - this->ibin_; }
  bool operator!=(const basic_reverse_iterator& o) const {
    return this->ibin_ != o.ibin_; }
  bool operator==(const basic_reverse_iterator& o) const {
    return this->ibin_ == o.ibin_; }
  bool operator<(const basic_reverse_iterator& o) const {
    return o.ibin_ < this->ibin_; }
  bool operator<=(const basic_reverse_iterator& o) const {
    return o.ibin_ <= this->ibin_; }
  bool operator>(const basic_reverse_iterator& o) const {
    return  o.ibin_ > this->ibin_; }
  bool operator>=(const basic_reverse_iterator& o) const {
    return o.ibin_ >= this->ibin_; }
};
#endif

// ============================================================================
//
// One-dimensional histogram class templated by accumulator
//
// ============================================================================

template<typename Acc> class BasicHistogram1D:
      public BinnedData1D<Acc, std::deque<Acc>>
{
#ifdef SWIG
  typedef BinnedData1D<Acc,std::deque<Acc>> Base;
#else
  using Base = BinnedData1D<Acc,std::deque<Acc>>;
#endif
 public:
  CALIN_TYPEALIAS(accumulator_type, Acc);

#ifndef SWIG
  class bin_accessor : public basic_bin_accessor<BasicHistogram1D>
  {
   public:
    using basic_bin_accessor<BasicHistogram1D>::basic_bin_accessor;
    bin_accessor& operator=(const bin_accessor&) = default;
    double weight() const { return this->binner_->weight(this->ibin_); }
  };

  class const_bin_accessor : public basic_bin_accessor<const BasicHistogram1D>
  {
   public:
    using basic_bin_accessor<const BasicHistogram1D>::basic_bin_accessor;
    double weight() const { return this->binner_->weight(this->ibin_); }
  };

  using iterator =
      basic_iterator<BasicHistogram1D,bin_accessor>;
  using reverse_iterator =
      basic_reverse_iterator<BasicHistogram1D,bin_accessor>;

  using const_iterator =
      basic_iterator<const BasicHistogram1D,const_bin_accessor>;
  using const_reverse_iterator =
      basic_reverse_iterator<const BasicHistogram1D,const_bin_accessor>;
#endif

  BasicHistogram1D(double dxval, double xval_align = 0.5,
                   const std::string& name = std::string(),
                   const std::string& xval_units = std::string(),
                   const std::string& weight_units = std::string()):
      Base(dxval, xval_align, xval_units),
      name_{name}, weight_units_{weight_units} { /* nothing to see here */ }

  BasicHistogram1D(double dxval, double xval_limit_lo, double xval_limit_hi,
                   double xval_align = 0.5,
                   const std::string& name = std::string(),
                   const std::string& xval_units = std::string(),
                   const std::string& weight_units = std::string()):
      Base(dxval, xval_limit_lo, xval_limit_hi, xval_align, xval_units),
      name_{name}, weight_units_{weight_units}
    { /* nothing to see here */ }

  BasicHistogram1D(const calin::ix::math::histogram::Histogram1DConfig& config):
      Base{config}, name_{config.name()}, weight_units_{config.weight_units()}
    { /* nothing to see here */ }

  BasicHistogram1D(const calin::ix::math::histogram::AccumulatedAndSerializedHistogram1DConfig& config):
      Base{config}, name_{config.name()}, weight_units_{config.weight_units()},
      compatify_on_serialize_{config.compactify_output()},
      serialize_max_dense_bins_in_output_ { config.max_dense_bins_in_output() },
      serialize_max_output_rebinning_ { config.max_output_rebinning() }
    { /* nothing to see here */ }

  BasicHistogram1D(const calin::ix::math::histogram::Histogram1DData& data);

  BasicHistogram1D(const BasicHistogram1D& o) = default;
#ifndef SWIG
  BasicHistogram1D& operator=(const BasicHistogram1D& o) = default;
#endif

#ifndef SWIG
  // Get all data as protobuf message
  calin::ix::math::histogram::Histogram1DData*
  dump_as_proto(calin::ix::math::histogram::Histogram1DData* data = nullptr)const;
#else
  calin::ix::math::histogram::Histogram1DData* dump_as_proto() const;
  void dump_as_proto(calin::ix::math::histogram::Histogram1DData* data) const;
#endif

#ifndef SWIG
  calin::ix::math::histogram::Histogram1DData*
  dump_as_compactified_proto(int max_dense_bins_in_output, int max_output_rebinning,
    calin::ix::math::histogram::Histogram1DData* data = nullptr) const;
  calin::ix::math::histogram::Histogram1DData*
  serialize(calin::ix::math::histogram::Histogram1DData* data = nullptr) const;
#else
  calin::ix::math::histogram::Histogram1DData*
  dump_as_compactified_proto(int max_dense_bins_in_output, int max_output_rebinning);
  void dump_as_compactified_proto(int max_dense_bins_in_output, int max_output_rebinning,
    calin::ix::math::histogram::Histogram1DData* data) const;

  calin::ix::math::histogram::Histogram1DData* serialize() const;
  void serialize(calin::ix::math::histogram::Histogram1DData* data) const;
#endif

  static BasicHistogram1D*
    create_from_proto(calin::ix::math::histogram::Histogram1DData& data) {
    return new BasicHistogram1D(data); }

  // Getters and setters
  std::string name() const { return name_; }
  std::string weight_units() const { return weight_units_; }
  void set_name(const std::string& name) { name_=name; }
  void set_weight_units(const std::string& units) { weight_units_=units; }

  // Functions to get size of histogram and clear it
  void clear() { Base::clear(); sum_w_={}; sum_wx_={}; sum_wxx_={};
    min_x_ = std::numeric_limits<double>::infinity();
    max_x_ = -std::numeric_limits<double>::infinity();
  }

  // Insert x value weight into histogram
  inline bool insert(const double x, const double w = 1.0);

  bool insert_hist(const BasicHistogram1D& h);

  unsigned insert_vec(const std::vector<double>& x, const double w = 1.0)
  {
    unsigned insert_count { 0 };
    for(auto ix : x)if(insert(ix,w))insert_count++;
    return insert_count;
  }

  unsigned insert_two_vec(const std::vector<double>& x,
                          const std::vector<double>& w)
  {
    unsigned insert_count { 0 };
    for(unsigned i=0; i<std::min(x.size(), w.size()); i++)
      if(insert(x[i],w[i]))insert_count++;
    return insert_count;
  }

  // Retrieve value for bin
  double weight(int ibin) const { return this->bin(ibin).total(); }
  double checked_weight(int ibin) const { return this->checked_bin(ibin).total(); }
  double weight_overflow_lo() const { return this->overflow_lo().total(); }
  double weight_overflow_hi() const { return this->overflow_hi().total(); }
  Eigen::VectorXd all_weight() const { Eigen::VectorXd w(this->size());
    for(int ibin=0;ibin<this->size();ibin++)w[ibin] = weight(ibin);
    return w; }

  // Raw access to the accumulators
  accumulator_type& accumulator(int ibin) { return this->bin(ibin); }
  const accumulator_type& accumulator(int ibin) const { return this->bin(ibin); }
  accumulator_type& checked_accumulator(int ibin) { return this->checked_bin(ibin); }
  const accumulator_type& checked_accumulator(int ibin) const { return this->checked_bin(ibin); }

#ifndef SWIG
  // Accessors
  bin_accessor accessor(int ibin) { return bin_accessor{*this,ibin}; }
  const const_bin_accessor accessor(int ibin) const { return const_bin_accessor{*this,ibin}; }

  // Iterator functions
  iterator begin() {
    return iterator{*this, 0}; }
  iterator end() {
    return iterator{*this, this->size()}; }
  reverse_iterator rbegin() {
    return reverse_iterator{*this, 0}; }
  reverse_iterator rend() {
    return reverse_iterator{*this, this->size()}; }

  const_iterator begin() const {
    return const_iterator{*this, 0}; }
  const_iterator end() const {
    return const_iterator{*this, this->size()}; }
  const_iterator cbegin() const {
    return const_iterator{*this, 0}; }
  const_iterator cend() const {
    return const_iterator{*this, this->size()}; }
  const_reverse_iterator crbegin() const {
    return const_reverse_iterator{*this, 0}; }
  const_reverse_iterator crend() const {
    return const_reverse_iterator{*this, this->size()}; }
#endif

  // Moments
  double sum_w() const { return sum_w_.total(); }
  double sum_wx() const { return sum_wx_.total(); }
  double sum_wxx() const { return sum_wxx_.total(); }
  double mean() const { return sum_wx()/sum_w(); }
  double var() const { return sum_wxx()/sum_w()-mean()*mean(); }
  double std() const { return std::sqrt(var()); }

  // Extreme values
  double min_xval() const { return min_x_; }
  double max_xval() const { return max_x_; }

  // Sum function over bins
  template<typename Fcn, typename IntAcc = Acc>
      double summation(const Fcn& fcn) const {
    IntAcc acc;
    for(auto& ibin : *this)
      acc.accumulate(fcn(ibin.xval_center(), ibin.weight()));
    return acc.total();
  }

  // Integrate function over bins
  template<typename Fcn, typename IntAcc = Acc>
      double integrate(const Fcn& fcn) const {
    return summation<Fcn,IntAcc>(fcn)*this->dxval_;
  }

  template<typename Fcn> void visit_bins(const Fcn& fcn) const {
    for(auto ibin : *this)fcn(ibin.xval_center(), ibin.weight()); }

  // Equality test
  bool operator==(const BasicHistogram1D& o) const;

 private:
  Acc sum_w_;
  Acc sum_wx_;
  Acc sum_wxx_;
  double min_x_ = std::numeric_limits<double>::infinity();
  double max_x_ = -std::numeric_limits<double>::infinity();
  std::string name_;
  std::string weight_units_;
  bool compatify_on_serialize_ = true;
  int serialize_max_dense_bins_in_output_ = 0;
  int serialize_max_output_rebinning_ = 1;
};

template<typename Acc> bool BasicHistogram1D<Acc>::
insert(const double x, const double w)
{
  if(!std::isfinite(w))return false;
  bool x_exceeds_limits { false };
  Acc& bin = this->bin_with_extend(x, x_exceeds_limits);
  if(!x_exceeds_limits)
  {
    double wxx { w };
    sum_w_.accumulate(wxx);
    sum_wx_.accumulate(wxx *= x);
    sum_wxx_.accumulate(wxx *= x);
    min_x_ = std::min(min_x_, x);
    max_x_ = std::max(max_x_, x);
  }

  bin.accumulate(w);
  return true;
}

template<typename Acc> bool BasicHistogram1D<Acc>::insert_hist(const BasicHistogram1D& h)
{
  if((h.dxval_ != this->dxval_) or (h.xval_align_ != this->xval_align_)
      or (h.limited_ != this->limited_)
      or (h.xval_limit_lo_ != this->xval_limit_lo_)
      or (h.xval_limit_hi_ != this->xval_limit_hi_)) {
    return false;
  }
  if(this->bins_.empty()) {
    this->xval0_ = h.xval0_;
  }
  while(h.xval0_ < this->xval0_) {
    this->bins_.emplace_front();
    this->xval0_ -= this->dxval_;
  }
  int offset = int(round((h.xval0_ - this->xval0_)*this->dxval_inv_));
  while(h.bins_.size()+offset > this->bins_.size()) {
    this->bins_.emplace_back();
  }
  for(unsigned ibin=0 ; ibin<h.bins_.size(); ibin++) {
    this->bins_.at(ibin + offset).accumulate(h.bins_[ibin]);
  }
  this->overflow_lo_.accumulate(h.overflow_lo_);
  this->overflow_hi_.accumulate(h.overflow_hi_);
  this->sum_w_.accumulate(h.sum_w_);
  this->sum_wx_.accumulate(h.sum_wx_);
  this->sum_wxx_.accumulate(h.sum_wxx_);
  this->min_x_ = std::min(min_x_, h.min_x_);
  this->max_x_ = std::max(max_x_, h.max_x_);
  return true;
}

template<typename Acc>
BasicHistogram1D<Acc>::
BasicHistogram1D(const calin::ix::math::histogram::Histogram1DData& data):
    Base(data.dxval(), data.xval_align(), data.xval0(),
         &(*data.bins().begin()),&(*data.bins().end()), data.limited(),
         data.xval_limit_lo(), data.xval_limit_hi(),
         data.overflow_lo(), data.overflow_hi(), data.xval_units()),
    sum_w_{data.sum_w()}, sum_wx_{data.sum_wx()}, sum_wxx_{data.sum_wxx()},
    min_x_{data.xval_min()}, max_x_{data.xval_max()},
    name_{data.name()}, weight_units_{data.weight_units()}
{
  int offset = 0;
  for(auto isparse: data.sparse_bins()) {
    int ibin = isparse.first - offset;
    while(ibin < 0) {
      this->bins_.emplace_front();
      ++ibin;
      --offset;
    }
    while(ibin >= int(this->bins_.size())) {
      this->bins_.emplace_back();
    }
    this->bins_[ibin].accumulate(isparse.second);
  }
  this->xval0_ += offset*this->dxval_;
}

template<typename Acc> calin::ix::math::histogram::Histogram1DData*
BasicHistogram1D<Acc>::
dump_as_proto(calin::ix::math::histogram::Histogram1DData* data) const
{
  if(data == nullptr)data = new calin::ix::math::histogram::Histogram1DData {};
  else data->Clear();
  data->set_dxval(this->dxval_);
  data->set_xval_align(this->xval_align_);
  data->set_xval0(this->xval0_);
  for(const auto& iacc : this->bins_)data->add_bins(iacc.total());
  data->set_limited(this->limited_);
  data->set_xval_limit_lo(this->xval_limit_lo_);
  data->set_xval_limit_hi(this->xval_limit_hi_);
  data->set_overflow_lo(this->overflow_lo_.total());
  data->set_overflow_hi(this->overflow_hi_.total());
  data->set_sum_w(sum_w_.total());
  data->set_sum_wx(sum_wx_.total());
  data->set_sum_wxx(sum_wxx_.total());
  data->set_xval_min(min_x_);
  data->set_xval_max(max_x_);
  data->set_name(name_);
  data->set_xval_units(this->xval_units_);
  data->set_weight_units(weight_units_);
  return data;
}


template<typename Acc>
bool BasicHistogram1D<Acc>::operator==(const BasicHistogram1D& o) const
{
  return
      this->dxval_ == o.dxval_ and
      this->xval_align_ == o.xval_align_ and
      this->xval0_ == o.xval0_ and
      this->bins_ == o.bins_ and
      this->limited_ == o.limited_ and
      this->xval_limit_lo_ == o.xval_limit_lo_ and
      this->xval_limit_hi_ == o.xval_limit_hi_ and
      this->overflow_lo_ == o.overflow_lo_ and
      this->overflow_hi_ == o.overflow_hi_ and
      sum_w_ == o.sum_w_ and
      sum_wx_ == o.sum_wx_ and
      sum_wxx_ == o.sum_wxx_ and
      min_x_ == o.min_x_ and
      max_x_ == o.max_x_ and
      name_ == o.name_ and
      this->xval_units_ == o.xval_units_ and
      weight_units_ == o.weight_units_;
}

// ============================================================================
//
// CDF
//
// ============================================================================

class BinnedCDF: public BinnedData1D<double>
{
#ifdef SWIG
  typedef BinnedData1D<double> Base;
#else
  using Base = BinnedData1D<double>;
#endif

 public:

#ifndef SWIG
  class const_bin_accessor : public basic_bin_accessor<const BinnedCDF>
  {
   public:
    using basic_bin_accessor<const BinnedCDF>::basic_bin_accessor;
    double density() const { return this->binner_->density(this->ibin_); }
    double cumulative_right() const { return this->binner_->cumulative_right(this->ibin_); }
    double cumulative_left() const { return this->binner_->cumulative_left(this->ibin_); }
  };
  using const_iterator =
      basic_iterator<const BinnedCDF,const_bin_accessor>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  using bin_accessor = const_bin_accessor;
  using iterator = const_iterator;
  using reverse_iterator = std::reverse_iterator<iterator>;
#endif

  BinnedCDF(double dxval, double xval_align = 0.5):
      Base{dxval, xval_align} { }

  template<typename Acc> BinnedCDF(const BasicHistogram1D<Acc>& hist):
      Base{hist.dxval(), hist.xval_align(), hist.xval0(),
        hist.is_limited(), hist.xval_limit_lo(), hist.xval_limit_hi(),
        hist.xval_units()}
  {
    Acc cumsum;
    for(const auto& ibin : hist)
    {
      cumsum.accumulate(ibin.weight());
      this->bins_.emplace_back(cumsum.total());
    }
    double total_weight = cumsum.total();
    for(auto& ibin : this->bins_)ibin /= total_weight;
    overflow_hi_ = 1.0 + hist.weight_overflow_hi() / total_weight;
    overflow_lo_ = -hist.weight_overflow_lo() / total_weight;
    set_name(hist.name());
  }

#ifdef SWIG
  BinnedCDF(const BasicHistogram1D<DefaultAccumulator>& hist):
      BinnedCDF<DefaultAccumulator>(hist) { }
#endif

  // Getters and setters
  std::string name() const { return name_; }
  void set_name(const std::string& name) { name_=name; }

  // Statistics
  template<typename IntAcc = DefaultAccumulator>
  void moments2(double& m1, double& m2) const {
    IntAcc ax {};
    IntAcc axx {};
    visit_delta([&ax,&axx](double xc, double w){ double wx {w*xc};
        ax.accumulate(wx); axx.accumulate(wx*xc); });
    m1 = ax.total();
    m2 = axx.total();
  }
  void mean_and_variance(double& mean, double var) const {
    moments2(mean,var);
    var -= mean*mean;
  }
  double mean() const {
    DefaultAccumulator ax {};
    visit_delta([&ax](double xc, double w){ ax.accumulate(xc*w); });
    return ax.total();
  }
  double variance() const {
    double mean = 0;
    double var = 0;
    mean_and_variance(mean, var);
    return var;
  }
  double rms() const { return std::sqrt(variance()); }

  double median() const { return quantile(0.5); }
  double quantile(double q) const {
    auto lhb = this->bins_.cbegin();
    return quantile_with_lhb(q, lhb);
  }

  // Retrieve cumulative at edges of bin and density in bin
  double cumulative_right(int ibin) const { return this->bin(ibin); }
  double checked_cumulative_right(int ibin) const { return this->checked_bin(ibin); }
  double cumulative_left(int ibin) const { return (ibin==0)?0.0:this->bin(ibin-1); }
  double checked_cumulative_left(int ibin) const { return (ibin==0)?0.0:this->checked_bin(ibin-1); }
  double cumulative_overflow_lo() const { return this->overflow_lo(); }
  double cumulative_overflow_hi() const { return this->overflow_hi(); }

  double density(int ibin) const {
    double cright = cumulative_right(ibin);
    double cleft = cumulative_left(ibin);
    return (cright-cleft)/this->dxval();
  }

  double checked_density(int ibin) const {
    double cright = checked_cumulative_right(ibin);
    double cleft = cumulative_left(ibin); // no need to recheck
    return (cright-cleft)/this->dxval();
  };

  double delta(int ibin) const {
    double cright = cumulative_right(ibin);
    double cleft = cumulative_left(ibin);
    return cright-cleft;
  }

  double checked_delta(int ibin) const {
    double cright = checked_cumulative_right(ibin);
    double cleft = cumulative_left(ibin); // no need to recheck
    return cright-cleft;
  };

  inline double interpolate_cumulative(double xval)
  {
    if(this->is_limited())
    {
      if(xval<this->xval_limit_lo())return 0.0;
      if(xval>=this->xval_limit_hi())return 1.0;
    }
    int thebin { this->ibin(xval) };
    if(thebin<0)return 0.0;
    if(thebin>=this->size())return 1.0;
    double y0 { thebin==0?0.0:cumulative_left(thebin) };
    double y1 { cumulative_right(thebin) };
    double dx = xval - xval_left(thebin);
    return y0*(1.0-dx) + y1*dx;
  }

  // Integration over bins

  template<typename Fcn> void visit_density(const Fcn& fcn) const
  {
    double cleft = 0;
    for(auto ibin : *this)
    {
      double cright = ibin.cumulative_right();
      fcn(ibin.xval_center(), (cright-cleft)/this->dxval_);
      cleft = cright;
    }
  }

  template<typename Fcn> void visit_delta(const Fcn& fcn) const
  {
    double cleft = 0;
    for(auto ibin : *this)
    {
      double cright = ibin.cumulative_right();
      fcn(ibin.xval_center(), cright-cleft);
      cleft = cright;
    }
  }

  template<typename Fcn, typename IntAcc = DefaultAccumulator>
  double integrate_density(const Fcn& fcn) const
  {
    IntAcc acc;
    double cleft = 0;
    for(auto ibin : *this)
    {
      double cright = ibin.cumulative_right();
      acc.accumulate(fcn(ibin.xval_center(), (cright-cleft)/this->dxval_));
      cleft = cright;
    }
    return acc.total()*this->dxval_;
  }

  template<typename Fcn, typename IntAcc = DefaultAccumulator>
  double integrate_delta(const Fcn& fcn) const
  {
    IntAcc acc;
    double cleft = 0;
    for(auto ibin : *this)
    {
      double cright = ibin.cumulative_right();
      acc.accumulate(fcn(ibin.xval_center(), cright-cleft));
      cleft = cright;
    }
    return acc.total();
  }

#ifndef SWIG
  // Iterator functions
  iterator begin() { return iterator{*this, 0}; }
  iterator end() { return iterator{*this, this->size()}; }
  reverse_iterator rbegin() { return reverse_iterator{end()}; }
  reverse_iterator rend() { return reverse_iterator{begin()}; }

  const_iterator begin() const { return const_iterator{*this, 0}; }
  const_iterator end() const { return const_iterator{*this, this->size()}; }
  const_iterator cbegin() const { return const_iterator{*this, 0}; }
  const_iterator cend() const { return const_iterator{*this, this->size()}; }
  const_reverse_iterator crbegin() const {
    return const_reverse_iterator{end()}; }
  const_reverse_iterator crend() const {
    return const_reverse_iterator{begin()}; }
#endif

 protected:
  double quantile_with_lhb(double q, data_container_type::const_iterator& lhb)
      const
  {
    lhb = std::lower_bound(lhb, this->bins_.cend(), q);
    int ibin = lhb-this->bins_.cbegin();
    double y0 = (ibin==0)?0:*(lhb-1);
    double y1 = *lhb;
    double x0 = xval_left(ibin);
    double x1 = xval_right(ibin);
    return (q-y0)*(x1-x0)/(y1-y0)+x0;
  }

 protected:
  std::string name_;
};

#ifdef SWIG
%template (SimpleHist) BasicHistogram1D<DefaultAccumulator>;
#endif

CALIN_TYPEALIAS(SimpleHist, BasicHistogram1D<DefaultAccumulator>);
CALIN_TYPEALIAS(Histogram1D, BasicHistogram1D<DefaultAccumulator>);
CALIN_TYPEALIAS(KahanHist, BasicHistogram1D<accumulator::KahanAccumulator>);

Histogram1D* new_histogram_if_enabled(
  const calin::ix::math::histogram::AccumulatedAndSerializedHistogram1DConfig& config);

#ifndef SWIG
calin::ix::math::histogram::Histogram1DData* compactify(
  const calin::ix::math::histogram::Histogram1DData& original_hist,
  int max_dense_bins_in_output, int max_output_rebinning,
  calin::ix::math::histogram::Histogram1DData* compactified_hist = nullptr);
#else
calin::ix::math::histogram::Histogram1DData* compactify(
  const calin::ix::math::histogram::Histogram1DData& original_hist,
  int max_dense_bins_in_output, int max_output_rebinning);
void compactify(
  const calin::ix::math::histogram::Histogram1DData& original_hist,
  int max_dense_bins_in_output, int max_output_rebinning,
  calin::ix::math::histogram::Histogram1DData* compactified_hist);
#endif

#ifndef SWIG
calin::ix::math::histogram::Histogram1DData*
rebin(const calin::ix::math::histogram::Histogram1DData& original_hist, int rebinning_factor,
  calin::ix::math::histogram::Histogram1DData* rebinned_hist = nullptr);
#else
calin::ix::math::histogram::Histogram1DData*
rebin(const calin::ix::math::histogram::Histogram1DData& original_hist, int rebinning_factor);
void  rebin(const calin::ix::math::histogram::Histogram1DData& original_hist, int rebinning_factor,
  calin::ix::math::histogram::Histogram1DData* rebinned_hist);
#endif

#ifndef SWIG
calin::ix::math::histogram::Histogram1DData*
sparsify(const calin::ix::math::histogram::Histogram1DData& original_hist,
  calin::ix::math::histogram::Histogram1DData* sparsified_hist = nullptr);
#else
calin::ix::math::histogram::Histogram1DData*
sparsify(const calin::ix::math::histogram::Histogram1DData& original_hist);
void sparsify(const calin::ix::math::histogram::Histogram1DData& original_hist,
  calin::ix::math::histogram::Histogram1DData* sparsified_hist);
#endif

#ifndef SWIG
calin::ix::math::histogram::Histogram1DData*
densify(const calin::ix::math::histogram::Histogram1DData& original_hist,
  calin::ix::math::histogram::Histogram1DData* densified_hist = nullptr);
#else
calin::ix::math::histogram::Histogram1DData*
densify(const calin::ix::math::histogram::Histogram1DData& original_hist);
void densify(const calin::ix::math::histogram::Histogram1DData& original_hist,
  calin::ix::math::histogram::Histogram1DData* densified_hist);
#endif

template<typename Acc>
calin::ix::math::histogram::Histogram1DData* BasicHistogram1D<Acc>::
dump_as_compactified_proto(int max_dense_bins_in_output, int max_output_rebinning,
  calin::ix::math::histogram::Histogram1DData* data) const
{
  std::unique_ptr<calin::ix::math::histogram::Histogram1DData> hist { dump_as_proto() };
  return compactify(*hist, max_dense_bins_in_output, max_output_rebinning, data);
}

template<typename Acc>
calin::ix::math::histogram::Histogram1DData* BasicHistogram1D<Acc>::
serialize(calin::ix::math::histogram::Histogram1DData* data) const
{
  if(compatify_on_serialize_) {
    return dump_as_compactified_proto(serialize_max_dense_bins_in_output_,
      serialize_max_output_rebinning_, data);
  } else {
    return dump_as_proto(data);
  }
}

} } } // namespace calin::math::histogram
