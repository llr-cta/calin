/*

   calin/diagnostics/functional_stats.cpp -- Stephen Fegan -- 2016-05-10

   Functional diagnostics visitor - integrated stats with histo and
   cross-channel covariance

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

#define CALIN_DIAGNOSTICS_FUNCTIONAL_NO_EXTERN

#include <algorithm>

#include <util/log.hpp>
#include <math/special.hpp>
#include <math/covariance_calc.hpp>
#include <diagnostics/functional.hpp>

using calin::math::special::SQR;
using namespace calin::util::log;
using namespace calin::diagnostics::functional;
using calin::iact_data::functional_event_visitor::
  DualGainInt32FunctionalTelescopeEventVisitor;
using calin::ix::diagnostics::functional::
  FunctionalStatsVisitorConfig;
using calin::math::covariance_calc::cov_i64_gen;

template
class calin::diagnostics::functional::
  BasicSingleFunctionalValueSupplierVisitor<int32_t>;

template
class calin::diagnostics::functional::
  BasicSingleFunctionalValueSupplierVisitor<double>;

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
