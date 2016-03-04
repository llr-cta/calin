/*

   calin/diagnostics/counter.hpp -- Stephen Fegan -- 2016-03-04

   Counter diagnostics visitor

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#pragma once

#include <diagnostics/counter.pb.h>
#include <iact_data/event_visitor.hpp>

namespace calin { namespace diagnostics { namespace counter {

class CounterDiagnostics:
  public iact_data::event_visitor::TelescopeEventVisitor
{
public:
  CounterDiagnostics();
  virtual ~CounterDiagnostics();

  bool visit_telescope_event(
    calin:: ix::iact_data::telescope_event::TelescopeEvent* event) override;

  calin::ix::diagnostics::counter::EventCounterDiagnosticsData& event_data() {
    return event_data_; }
private:
  uint64_t next_event_number_ = 0;
  calin::ix::diagnostics::counter::EventCounterDiagnosticsData event_data_;
};

} } } /// namespace calin::diagnostics::counter_diagnostics
