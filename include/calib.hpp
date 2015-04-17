/* 

   calin/math.hpp -- Stephen Fegan -- 2015-04-17

   Header file to simplify including of calin math components. Loads
   compoenent header files and imports various symbols into calin::math
   namespace for ease of use.

   THIS FILE SHOULD NOT BE INCLUDED BY ANY CALIN HPP OR CPP FILE - IT
   IS ONLY FOR USE BY END USER CODE.

*/

#include "package_wide_definitions.hpp"

#include "calib/spe_fit.hpp"

namespace calin { namespace calib {

using spe_fit::SingleElectronSpectrum;
using spe_fit::MultiElectronSpectrum;
using spe_fit::PoissonGaussianMES;
using spe_fit::PoissonGaussianMES_HighAccuracy;
using spe_fit::GeneralPoissonMES;
using spe_fit::SPELikelihood;

} } // namespace calin::calib
