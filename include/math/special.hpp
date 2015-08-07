/* 

   calin/math/special.hpp -- Stephen Fegan -- 2015-08-06

   Various special functions, which could be implemented as interfaces 
   to GSL or directly

*/

#pragma once

#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_erf.h>

namespace calin { namespace math { namespace special {

inline double SQR(double x) { return x*x; }
inline double CUBE(double x) { return x*x*x; }
inline double QUAD(double x) { return SQR(SQR(x)); }

inline double dawson(double x)
{
  return gsl_sf_dawson(x);
}

inline double lerfc(double x)
{
  return gsl_sf_log_erfc(x);
}
  
} } } // namespace calin::math::special
