/* 

   calin/math/histogram.cpp -- Stephen Fegan -- 2015-02-16

*/

#include "math/accumulator.hpp"
#include "math/histogram.hpp"

namespace calin { namespace math { namespace histogram {

template class BasicHistogram1D<accumulator::SimpleAccumulator>;

} } } // namespace calin::math::histogram
