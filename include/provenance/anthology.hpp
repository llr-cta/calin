/*

   calin/provenance/anthology.hpp -- Stephen Fegan -- 2017-10-30

   Anthology of all provenance system info

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <provenance/anthology.pb.h>

namespace calin { namespace provenance { namespace anthology {

#ifndef SWIG
calin::ix::provenance::anthology::Anthology* get_current_anthology(calin::ix::provenance::anthology::Anthology* x = nullptr);
#else
calin::ix::provenance::anthology::Anthology* get_current_anthology();
void get_current_anthology(calin::ix::provenance::anthology::Anthology* x);
#endif

} } } // namespace calin::provenance::anthology
