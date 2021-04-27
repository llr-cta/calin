/*

   calin/simulation/sct_ray_tracer.hpp -- Stephen Fegan -- 2021-04-27

   Class for SCT ray tracing

   Copyright 2021, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <cmath>
#include <Eigen/Dense>

#include <simulation/sct_optics.pb.h>
#include <simulation/sct_facet_scheme.hpp>

namespace calin { namespace simulation { namespace sct_optics {

calin::ix::simulation::sct_optics::SCTArray*
make_sct_array(calin::ix::simulation::sct_optics::SCTRandomArrayParameters& param);

} } } // namespace calin::simulations::sct_optics
