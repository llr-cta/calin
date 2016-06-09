/*

   calin/diagnostics/functional_stats.cpp -- Stephen Fegan -- 2016-05-10

   Functional diagnostics visitor - integrated stats with histo and
   cross-channel covariance

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

#define CALIN_DIAGNOSTICS_FUNCTIONAL_NO_EXTERN

#include <algorithm>

#include <io/log.hpp>
#include <math/special.hpp>
#include <math/covariance_calc.hpp>
#include <diagnostics/functional.hpp>

using calin::math::special::SQR;
using namespace calin::io::log;
using namespace calin::diagnostics::functional;
using calin::iact_data::functional_event_visitor::
  DualGainInt32FunctionalTelescopeEventVisitor;
using calin::ix::diagnostics::functional::
  FunctionalStatsVisitorConfig;
using calin::math::covariance_calc::cov_i64_gen;

SingleFunctionalValueSupplierVisitor::SingleFunctionalValueSupplierVisitor(
    calin::iact_data::functional_event_visitor::
      DualGainInt32FunctionalTelescopeEventVisitor* value_supplier,
    unsigned chan, bool low_gain):
  ValueSupplierVisitor<int32_t>(), value_supplier_(value_supplier),
  chan_(chan), low_gain_(low_gain)
{
  // nothing to see here
}

SingleFunctionalValueSupplierVisitor::~SingleFunctionalValueSupplierVisitor()
{
  // nothing to see here
}

bool SingleFunctionalValueSupplierVisitor::
visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  has_value_ = false;
  return true;
}

bool SingleFunctionalValueSupplierVisitor::visit_waveform(unsigned ichan,
  calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
  calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain)
{
  if(ichan==chan_)
  {
    if(low_gain_) {
      if(low_gain)
        has_value_ = true, value_ = value_supplier_->low_gain_value();
    } else {
      if(high_gain)
        has_value_ = true, value_ = value_supplier_->high_gain_value();
    }
  }
  return true;
}

bool SingleFunctionalValueSupplierVisitor::get_value(int32_t& value)
{
  if(has_value_)value = value_;
  return has_value_;
}

template
class calin::diagnostics::functional::FunctionalStatsVisitor<
  calin::iact_data::functional_event_visitor::
    DualGainInt32FunctionalTelescopeEventVisitor,
  calin::ix::diagnostics::functional::CameraIntFunctionalRawStats>;

template
class calin::diagnostics::functional::FunctionalStatsVisitor<
  calin::iact_data::functional_event_visitor::
    DualGainDoubleFunctionalTelescopeEventVisitor,
  calin::ix::diagnostics::functional::CameraDoubleFunctionalRawStats>;
