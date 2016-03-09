/*

   calin/diagnostics/event_number.hpp -- Stephen Fegan -- 2016-03-09

   Event number diagnostics visitors

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

#include <diagnostics/event_number.pb.h>
#include <iact_data/event_visitor.hpp>

namespace calin { namespace diagnostics { namespace event_number {

class SequentialEventNumberGlitchDetector:
  public iact_data::event_visitor::TelescopeEventVisitor
{
public:
  SequentialEventNumberGlitchDetector(bool test_local_event_number = false);
  virtual ~SequentialEventNumberGlitchDetector();

  bool visit_telescope_event(
    calin:: ix::iact_data::telescope_event::TelescopeEvent* event) override;

  calin::ix::diagnostics::event_number::
  SequentialEventNumberGlitchDetectorData& glitch_data() {
    return glitch_data_; }

private:
  bool test_local_event_number_ = false;
  int64_t last_event_number_ = -1;
  calin::ix::diagnostics::event_number::
  SequentialEventNumberGlitchDetectorData glitch_data_;
};

} } } /// namespace calin::diagnostics::event_number
