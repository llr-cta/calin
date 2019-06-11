/*

   calin/iact_data/lstcam_layout.hpp -- Stephen Fegan -- 2018-10-05

   Camera layout for LSTCam

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

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

#include <calin_global_definitions.hpp>
#include <calin_global_config.hpp>
#include <iact_data/instrument_layout.pb.h>

namespace calin { namespace iact_data { namespace lstcam_layout {

calin::ix::iact_data::instrument_layout::CameraLayout*
lstcam_layout(
  calin::ix::iact_data::instrument_layout::CameraLayout* layout = nullptr);

} } } // namespace calin::iact_data::lstcam_layout
