/*

   calin/diagnostics/delta_t.hpp -- Stephen Fegan -- 2016-06-03

   Delta-T visitor

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

#include <iact_data/event_visitor.hpp>
#include <diagnostics/value_capture.hpp>

namespace calin { namespace diagnostics { namespace delta_t {

class OneModuleCounterDeltaTVisitor :
  public value_capture::ValueSupplierVisitor<double>
{
public:
  CALIN_TYPEALIAS(value_type, double);
  OneModuleCounterDeltaTVisitor(int module_id = 0, int counter_index = 0,
    double conversion_to_ns = 1):
      ValueSupplierVisitor<double>(), module_id_(module_id),
      counter_index_(counter_index), conversion_to_ns_(conversion_to_ns)
    { /* nothing to see here */ }
  virtual ~OneModuleCounterDeltaTVisitor() { /* nothing to see here */ }
  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;
  bool get_value(double& value) override;
protected:
  int module_id_ = 0;
  int counter_index_ = 0;
  double conversion_to_ns_ = 1;
  int64_t last_counter_val_ = 0;
  int64_t last_counter_event_number_ = -1;
  bool has_delta_t_val_ = false;
  double delta_t_val_ = 0;
};


} } } // namespace calin::diagnostics::delta_t
