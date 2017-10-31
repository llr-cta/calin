/*

   calin/diagnostics/value_capture.hpp -- Stephen Fegan -- 2016-05-23

   Value capture diagnostics

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <diagnostics/value_capture.pb.h>
#include <iact_data/event_visitor.hpp>

namespace calin { namespace diagnostics { namespace value_capture {

template<typename T> class ValueSupplierVisitor:
  public iact_data::event_visitor::TelescopeEventVisitor
{
public:
  CALIN_TYPEALIAS(value_type, T);
  ValueSupplierVisitor():
    TelescopeEventVisitor() { /* nothing to see here */ }
  virtual ~ValueSupplierVisitor() { /* nothing to see here */ }
  virtual bool get_value(T& value) = 0;
};

template<typename Supplier, typename Results>
class SequentialValueCaptureVisitor:
  public iact_data::event_visitor::TelescopeEventVisitor
{
public:
  SequentialValueCaptureVisitor(Supplier* supplier,
      typename Supplier::value_type guard_value = 0,
      bool inherit_supplier = false):
    TelescopeEventVisitor(),
    supplier_(supplier), inherit_supplier_(inherit_supplier),
    guard_value_(guard_value), results_()
  {
    results_.set_guard_value(guard_value_);
  }

  virtual ~SequentialValueCaptureVisitor()
  {
    if(inherit_supplier_)delete supplier_;
  }

  bool demand_waveforms() override { return false; }

  bool leave_telescope_event() override
  {
    typename Supplier::value_type value = guard_value_;
    if(supplier_->get_value(value))results_.add_value(value);
    else results_.add_value(guard_value_);
    return true;
  }

  Results results() { return results_; }

private:
  Supplier* supplier_ = nullptr;
  bool inherit_supplier_ = false;
  typename Supplier::value_type guard_value_;
  Results results_;
};

} } } // namespace calin::diagnostics::value_capture

#ifndef SWIG
#ifndef CALIN_VALUE_CAPTURE_NO_EXTERN

extern template
class calin::diagnostics::value_capture::ValueSupplierVisitor<int32_t>;
extern template
class calin::diagnostics::value_capture::SequentialValueCaptureVisitor<
  calin::diagnostics::value_capture::ValueSupplierVisitor<int32_t>,
  calin::ix::diagnostics::value_capture::CapturedInt32Values>;

extern template
class calin::diagnostics::value_capture::ValueSupplierVisitor<int64_t>;
extern template
class calin::diagnostics::value_capture::SequentialValueCaptureVisitor<
  calin::diagnostics::value_capture::ValueSupplierVisitor<int64_t>,
  calin::ix::diagnostics::value_capture::CapturedInt64Values>;

extern template
class calin::diagnostics::value_capture::ValueSupplierVisitor<double>;
extern template
class calin::diagnostics::value_capture::SequentialValueCaptureVisitor<
  calin::diagnostics::value_capture::ValueSupplierVisitor<double>,
  calin::ix::diagnostics::value_capture::CapturedDoubleValues>;

#endif // #ifdef CALIN_VALUE_CAPTURE_NO_EXTERN
#endif // #ifdef SWIG
