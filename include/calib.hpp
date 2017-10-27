/* 

   calin/math.hpp -- Stephen Fegan -- 2015-04-17

   Header file to simplify including of calin math components. Loads
   compoenent header files and imports various symbols into calin::math
   namespace for ease of use.

   THIS FILE SHOULD NOT BE INCLUDED BY ANY CALIN HPP OR CPP FILE - IT
   IS ONLY FOR USE BY END USER CODE.

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include "calin_global_definitions.hpp"

#include "calib/spe_fit.hpp"

namespace calin { namespace calib {

using spe_fit::SingleElectronSpectrum;
using spe_fit::MultiElectronSpectrum;
using spe_fit::PoissonGaussianMES;
using spe_fit::PoissonGaussianMES_HighAccuracy;
using spe_fit::GeneralPoissonMES;
using spe_fit::SPELikelihood;

} } // namespace calin::calib
