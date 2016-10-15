/*

   calin/math/interpolation_1d.hpp -- Stephen Fegan -- 2016-10-14

   Simple class to do 1D linear or exponential interpolation and extrapolation.
   Originally from EGS5 ACT code, see header below.

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

// Interpolation1D.hpp - Simple class to do 1D linear interpolation and
// - extrapolation
// Stephen Fegan - sfegan@llr.in2p3.fr - Novermber 2012
// $Id: Interpolation1D.hpp 5422 2013-06-26 14:01:03Z sfegan $

#pragma once

#include <cmath>
#include <cassert>
#include <vector>
#include <string>
#include <algorithm>
#include <utility>

namespace calin { namespace math { namespace interpolation_1d {

template<typename T> class LinearInterpolator
{
public:
  inline
  T interpolate(double x, double x0, const T& y0, double x1, const T& y1) const
  {
    double z = (x-x0)/(x1-x0);
    return y0*(1.0-z) + y1*z;
  }

  inline
  T integrate(double x0, const T& y0, double x1, const T& y1) const
  {
    return (y1+y0)*(0.5*(x1-x0));
  }
};

template<typename T> class ExpInterpolator
{
public:
  inline
  T interpolate(double x, double x0, const T& y0, double x1, const T& y1) const
  {
    double z = (x-x0)/(x1-x0);
    return std::exp(std::log(y0)*(1.0-z) + std::log(y1)*z);
  }

  inline
  T integrate(double x0, const T& y0, double x1, const T& y1) const
  {
    if(y0 == y1)return y0*(x1-x0);
    return ((y1-y0)*(x1-x0))/(std::log(y1)-std::log(y0));
  }
};

template<typename T, class Interpolator> class Interpolation1D
{
public:
  typedef std::pair<double, T>            xy_type;
  typedef std::vector<xy_type>            xy_vec_type;
  typedef unsigned                        size_type;
  //typedef typename xy_vec_type::size_type size_type;
  Interpolation1D(T ydef=T(), Interpolator interp = Interpolator()):
    m_interp(interp), m_ydef(ydef), m_xy() { }

  Interpolation1D(const std::vector<double>& x, const std::vector<T>& y,
		  T ydef=T(), Interpolator interp = Interpolator())
    : m_interp(interp), m_ydef(ydef), m_xy()
  {
    if(x.size() != y.size())
      throw std::range_error("Interpolation1D: x and y vectors must have same size");
    for(unsigned ixy=0;ixy<x.size();ixy++)
      m_xy.emplace_back(x[ixy],y[ixy]);
    std::sort(m_xy.begin(), m_xy.end());
  }

  template<class BinOp>
  Interpolation1D(const Interpolation1D& f1, const Interpolation1D& f2,
		  BinOp op, double epsilon = 0,
		  Interpolator interp = Interpolator()):
    m_interp(interp), m_ydef(op(0,f1.m_ydef,f2.m_ydef)), m_xy()
  {
    for(typename xy_vec_type::const_iterator i1 =
      f1.m_xy.begin(), i2=f2.m_xy.begin();
      i1!=f1.m_xy.end() || i2!=f2.m_xy.end();)
    {
    	typename xy_vec_type::const_iterator i;
    	if(i2==f2.m_xy.end() || (i1!=f1.m_xy.end() && i1->first<i2->first))
    	  i = i1++;
    	else
    	  i = i2++;
    	if(m_xy.empty() || fabs(i->first - m_xy.back().first)>epsilon)
  	  {
  	    xy_type xy(i->first, op(i->first, f1(i->first), f2(i->first)));
  	    m_xy.push_back(xy);
  	  }
    }
  }

  class BinAddOp {
    public: T operator() (double x, const T& y1, const T& y2) { return y1+y2; }
  };
  class BinSubOp {
    public: T operator() (double x, const T& y1, const T& y2) { return y1-y2; }
  };
  class BinMulOp {
    public: T operator() (double x, const T& y1, const T& y2) { return y1*y2; }
  };
  class BinDivOp {
    public: T operator() (double x, const T& y1, const T& y2) { return y1/y2; }
  };

  template<class UnOp> void transform(UnOp op)
  {
    for(typename xy_vec_type::const_iterator i = m_xy.begin();
        i != m_xy.end(); i++)
      i->second = op(i->first, i->second);
  }

  class UnNegOp {
    public: T operator() (double x, const T& y) { return -y; }
  };
  class UnExpOp {
    public: T operator() (double x, const T& y) { return std::exp(y); }
  };

  void insert(double x, const T& y)
  {
    m_xy.emplace_back(x,y);
    std::sort(m_xy.begin(), m_xy.end());
  }

  void clear() { m_xy.clear(); }

  size_type nXY() const { return m_xy.size(); }
  const xy_type& xyi(size_type ixy) { return m_xy[ixy]; }
  const double& xi(size_type ixy) { return m_xy[ixy].first; }
  const T& yi(size_type ixy) { return m_xy[ixy].second; }

  double xmin() const { return m_xy.front().first; }
  double xmax() const { return m_xy.back().first; }

  T y(double x) const
  {
    if(m_xy.size() == 0)
      return m_ydef;
    else if(m_xy.size() == 1)
      return m_xy.front().second;

    typename xy_vec_type::const_iterator xy1 =
      std::upper_bound(m_xy.begin(), m_xy.end(), std::make_pair(x, T()));
    typename xy_vec_type::const_iterator xy0;

    if(xy1 == m_xy.begin()) {
      xy0=xy1; xy1++;
    } else if(xy1 == m_xy.end()) {
      xy1--; xy0=xy1; xy0--;
    } else {
      xy0=xy1; xy0--;
    }

    return m_interp.interpolate(x, xy0->first, xy0->second,
		  xy1->first, xy1->second);
  }

  T operator() (double x) const { return y(x); }

  T integrate() const
  {
    if(m_xy.size() < 2)return T();
    T sum = T();
    typename xy_vec_type::const_iterator ixy=m_xy.begin();
    double x0 = ixy->first;
    T y0 = ixy->second;
    ixy++;
    while(ixy!=m_xy.end())
    {
      double x1 = ixy->first;
      T y1 = ixy->second;
      sum += m_interp.integrate(x0,y0,x1,y1);
      x0 = x1;
      y0 = y1;
      ixy++;
    }
    return sum;
  }

  T integrate(double xlo, double xhi)
  {
    assert(xlo <= xhi);
    if(m_xy.size() == 0)
      return m_ydef * (xhi-xlo);
    else if(m_xy.size() == 1)
      return m_xy.front().second * (xhi-xlo);
    T sum = T();
    double x0 = xlo;
    T y0 = y(x0);
    typename xy_vec_type::const_iterator ixy =
      std::upper_bound(m_xy.begin(), m_xy.end(), xy_type(xlo, T()));
    while(ixy!=m_xy.end() && ixy->first<xhi)
    {
      double x1 = ixy->first;
      T y1 = ixy->second;
      sum += m_interp.integrate(x0,y0,x1,y1);
      x0 = x1;
      y0 = y1;
      ixy++;
    }
    double x1 = xhi;
    T y1 = y(x1);
    sum += m_interp.integrate(x0,y0,x1,y1);
    return sum;
  }

  Interpolation1D& operator +=(const Interpolation1D& o) {
    *this = Interpolation1D(*this, o, BinAddOp()); return *this;
  }
  Interpolation1D& operator -=(const Interpolation1D& o) {
    *this = Interpolation1D(*this, o, BinSubOp()); return *this;
  }
  Interpolation1D& operator *=(const Interpolation1D& o) {
    *this = Interpolation1D(*this, o, BinMulOp()); return *this;
  }
  Interpolation1D& operator /=(const Interpolation1D& o) {
    *this = Interpolation1D(*this, o, BinDivOp()); return *this;
  }

  Interpolation1D operator+ (const Interpolation1D& o) {
    return Interpolation1D(*this, o, BinAddOp());
  }
  Interpolation1D operator- (const Interpolation1D& o) {
    return Interpolation1D(*this, o, BinSubOp());
  }
  Interpolation1D operator* (const Interpolation1D& o) {
    return Interpolation1D(*this, o, BinMulOp());
  }
  Interpolation1D operator/ (const Interpolation1D& o) {
    return Interpolation1D(*this, o, BinDivOp());
  }

private:
  Interpolator m_interp;
  T            m_ydef;
  xy_vec_type  m_xy;
};

typedef Interpolation1D<double, LinearInterpolator<double> > InterpLinear1D;
typedef Interpolation1D<double, ExpInterpolator<double> > InterpExp1D;

} } } // namespace calin::math::interpolation_1d
