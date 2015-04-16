/* 

   calin/math/function.hpp -- Stephen Fegan -- 2015-02-24

   Base classes for functions and general parameterizable objects that
   can be used with optimizers, root finders, the MCMC algorithm etc.

*/

#include "package_wide_definitions.hpp"

#include "math/accumulator.hpp"
#include "math/historgram.hpp"

namespace calin { namespace math {

using accumulator::SimpleAccumulator;
using accumulator::KahanAccumulator;
using accumulator::Accumulator;
using accumulator::BasicAccumulator;
using accumulator::KahanAccumulator;

using histogram::BasicHistogram1D;
using histogram::Histogram1D;
using histogram::SimpleHist;
using histogram::BinnedCDF;

} } // namespace calin::math
