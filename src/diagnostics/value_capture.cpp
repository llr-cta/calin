/*

   calin/diagnostics/value_capture.cpp -- Stephen Fegan -- 2016-05-23

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

#include <diagnostics/value_capture.hpp>

template
class calin::diagnostics::value_capture::ValueSupplierVisitor<int32_t>;
template
class calin::diagnostics::value_capture::SequentialValueCaptureVisitor<
  calin::diagnostics::value_capture::ValueSupplierVisitor<int32_t>,
  calin::ix::diagnostics::value_capture::CapturedInt32Values>;

template
class calin::diagnostics::value_capture::ValueSupplierVisitor<int64_t>;
template
class calin::diagnostics::value_capture::SequentialValueCaptureVisitor<
  calin::diagnostics::value_capture::ValueSupplierVisitor<int64_t>,
  calin::ix::diagnostics::value_capture::CapturedInt64Values>;

template
class calin::diagnostics::value_capture::ValueSupplierVisitor<double>;
template
class calin::diagnostics::value_capture::SequentialValueCaptureVisitor<
  calin::diagnostics::value_capture::ValueSupplierVisitor<double>,
  calin::ix::diagnostics::value_capture::CapturedDoubleValues>;
