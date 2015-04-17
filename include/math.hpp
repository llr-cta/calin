/* 

   calin/math.hpp -- Stephen Fegan -- 2015-04-17

   Header file to simplify including of calin math components. Loads
   compoenent header files and imports various symbols into calin::math
   namespace for ease of use.

   THIS FILE SHOULD NOT BE INCLUDED BY ANY CALIN HPP OR CPP FILE - IT
   IS ONLY FOR USE BY END USER CODE.

*/

#include "package_wide_definitions.hpp"

#include "math/accumulator.hpp"
#include "math/historgram.hpp"
#include "math/function.hpp"
#include "math/pdf_1d.hpp"
#include "math/optimizer.hpp"
#include "math/hessian.hpp"
#include "math/nlopt_optimizer.hpp"

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

using function::ParameterAxis;
using function::DomainAxis;
using function::Parameterizable;
using function::MultiAxisFunction;
using function::SingleAxisFunction;
using function::ParameterizableMultiAxisFunction;
using function::ParameterizableSingleAxisFunction;

using pdf_1d::Parameterizable1DPDF;
using pdf_1d::GaussianPDF;
using pdf_1d::LimitedGaussianPDF;
using pdf_1d::LimitedExponentialPDF;

using optimizer::Optimizer;

using hessian::calculate_hessian;
using hessian::hessian_to_error_matrix;

using optimizer::NLOptOptimizer;

} } // namespace calin::math
