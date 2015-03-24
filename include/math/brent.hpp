/* 

   calin/math/brent.hpp -- Stephen Fegan -- 2015-03-24

   One dimensional root finder using Brent's method. This code is a
   wrapper around the C++ functions by John Burkardt

   http://people.sc.fsu.edu/~jburkardt/cpp_src/brent/brent.html

*/

#pragma once

namespace calin { namespace math { namespace brent {

class func_base{
public:
  virtual double operator() (double) = 0;
};

double glomin ( double a, double b, double c, double m, double e, double t,
  func_base& f, double &x );
double local_min ( double a, double b, double t, func_base& f,
  double &x );
double local_min_rc ( double &a, double &b, int &status, double value );
double zero ( double a, double b, double t, func_base& f );
void zero_rc ( double a, double b, double t, double &arg, int &status,
  double value );

// === simple wrapper functions
// === for convenience and/or compatibility
double glomin ( double a, double b, double c, double m, double e, double t,
  double f ( double x ), double &x );
double local_min ( double a, double b, double t, double f ( double x ),
  double &x );
double zero ( double a, double b, double t, double f ( double x ) );

template<typename F_of_X> class func_adapter: public func_base
{
 public:
  func_adapter(F_of_X f): func_base(), f_(f) { }
  virtual ~func_adapter() { }
  double operator() (double x) override { return f_(x); }
 private:
  F_of_X f_;
};

template<typename F_of_X>
double brent_zero(double xlo, double xhi, F_of_X f, double tol=0.001)
{
  func_adapter<F_of_X> adapter(f);
  return zero(xlo,xhi,tol,adapter);
}

} // namespace brent

using brent::brent_zero;

} } // namespace calin:math
