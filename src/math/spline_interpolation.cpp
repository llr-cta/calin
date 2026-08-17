/*

   calin/math/spline_interpolation.hpp -- Stephen Fegan -- 2018-12-17

   Spline interpolation functions.

   Some of this code from my Fermi tools code "RegularSpline.hpp",
   These portions are copyright Stephen Fegan, 2010.

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <stdexcept>
#include <iostream>
#include <Eigen/Dense>

#include <math/spline_interpolation.hpp>
#include <math/special.hpp>
#include <util/log.hpp>

using namespace calin::util::log;
using namespace calin::math::spline_interpolation;

using calin::math::special::SQR;
using calin::math::special::CUBE;
using calin::math::special::QUAD;

InterpolationIntervals
calin::math::spline_interpolation::make_intervals(const std::vector<double>& x)
{
  if(x.size() < 2)
    throw std::runtime_error("make_intervals: need at least two x values");
  InterpolationIntervals intervals;
  intervals.xmin = x[0];

  intervals.regular_xmax = x[1];
  intervals.regular_dx = x[1] - x[0];
  intervals.regular_dx_inv = 1/intervals.regular_dx;
  intervals.irregular_begin = 1;

  if(intervals.regular_dx <= 0)
    std::runtime_error("make_intervals: x values must be monotonously increasing");
  for(unsigned i=2; i<x.size(); i++) {
    double dx = x[i] - x[i-1];
    if(dx <= 0)
      std::runtime_error("make_intervals: x values must be monotonously increasing");
    if(intervals.irregular_begin == (i-1) and
        std::abs(dx - intervals.regular_dx)/intervals.regular_dx < 1e-6) {
      intervals.irregular_begin = i;
      intervals.regular_xmax = x[i];
    }
  }
  intervals.xmax = x.back();
  intervals.x = x;
  intervals.dx.resize(x.size()-1);
  std::transform(x.begin(), x.end()-1, x.begin()+1, intervals.dx.begin(),
    [](double xlo, double xhi){ return xhi-xlo; });
  intervals.dx_inv.resize(x.size()-1);
  std::transform(intervals.dx.begin(), intervals.dx.end(), intervals.dx_inv.begin(),
    [](double dx){ return 1.0/dx; });
  if(intervals.irregular_begin == 1) {
    intervals.regular_xmax = x[0];
    intervals.regular_dx = 0;
    intervals.regular_dx_inv = 0;
    intervals.irregular_begin = 0;
  }
  return intervals;
}

std::vector<double> calin::math::spline_interpolation::
generate_cubic_spline_interpolation(
  const std::vector<double>& x, const std::vector<double>& y,
  BoundaryConitions bc_lhs, double bc_lhs_val,
  BoundaryConitions bc_rhs, double bc_rhs_val)
{
  unsigned np = y.size();

  if(x.size() != np)
    throw std::runtime_error("generate_cubic_spline_interpolation: x and y must "
      "have equal number of points");

  if(np<3)
    throw std::runtime_error("generate_cubic_spline_interpolation: must "
      "have at least three points to interpolate");

  std::vector<double> dx_inv(np-1, 0); // 1/dx
  for(unsigned i=0; i<np-1; i++)dx_inv[i] = 1.0/(x[i+1]-x[i]);

  std::vector<double> dy(np-1, 0);
  for(unsigned i=0; i<np-1; i++)dy[i] = y[i+1]-y[i];

  // First solve the tridiagonal matrix equation for the slopes at
  // each of the points.. do this in a manner so that we can later
  // calculate the value of the spline and the derivative with respect
  // to the input points and slopes

  std::vector<double> a(np, 0); // Sub diagnoal elements   [1,np-1]
  std::vector<double> b(np, 0); // Diagonal elements       [0,np-1]
  std::vector<double> c(np, 0); // Super diagonal elements [0,np-2]
  std::vector<double> d(np, 0); // RHS of equation         [0,np-1]
  std::vector<double> m(np, 0); // Slopes we are seeking   [0,np-1]

  // Enter primary spline constraints into matrix
  for(unsigned i=1; i<np-1; i++) {
    b[i] = 2*(dx_inv[i]+dx_inv[i-1]);
    a[i] = dx_inv[i-1];
    c[i] = dx_inv[i];
    d[i] = 3*(dy[i]*SQR(dx_inv[i]) + dy[i-1]*SQR(dx_inv[i-1]));
  }

  // We cannot use the default NOT_A_KNOT algorithm with only 2 or 3
  // points (i.e. 1 or 2 spline segments), since there is no "knot" in
  // the first case, giving no constraints at all, and only one "knot"
  // in the second case giving a single constraint.

  // In the 2-point case, if either of the BCS is set to NOT_A_KNOT we
  // zero out the cubic term. If both are set to zero then we also set
  // the quadratic term to zero (by changing bcsN to BCS_NATURAL).
  // This has the pleasent effect of making the spline into a line
  // segment if either is NOT_A_KNOT and the other is either NATURAL
  // or NOT_A_KNOT. In the case of NOT_A_KNOT combined with CLAMPED
  // (preassigned slope) we get a quadratic.

  // In the 3-point case, if either of the BCS is set to NOT_A_KNOT we
  // set the quadratic term to be constant across the "knot" as
  // expected. The other constraint then depends on what the other
  // BCS is. For NATURAL and CLAMPED the other constraint is obvious.
  // If the other is also NOT_A_KNOT we set the cubic term to zero to
  // get a quadratic that fits the three points.

  switch(bc_lhs)
  {
  case BC_NATURAL:
    b[0] = 2;
    c[0] = 1;
    d[0] = 3*dy[0]*dx_inv[0];
    break;
  case BC_CLAMPED_SLOPE:
    b[0] = 1;
    c[0] = 0;
    d[0] = bc_lhs_val;
    break;
  case BC_NOT_A_KNOT:
    b[0] = dx_inv[0]*(dx_inv[1]+dx_inv[0]);
    c[0] = SQR(dx_inv[1]+dx_inv[0]);
    d[0] = dy[1]*CUBE(dx_inv[1]) + dy[0]*(3*dx_inv[1]+2*dx_inv[0])*SQR(dx_inv[0]);
    break;
  }

  switch(bc_rhs)
  {
  case BC_NATURAL:
    b[np-1] = 2;
    a[np-1] = 1;
    d[np-1] = 3*dy[np-2]*dx_inv[np-2];
    break;
  case BC_CLAMPED_SLOPE:
    b[np-1] = 1;
    a[np-1] = 0;
    d[np-1] = bc_rhs_val;
    break;
  case BC_NOT_A_KNOT:
    if((np==3)and(bc_lhs==BC_NOT_A_KNOT)) {
      b[2] = 1;
      a[2] = 1;
      d[2] = 2*dy[1]*dx_inv[1];
    } else {
      b[np-1] = dx_inv[np-2]*(dx_inv[np-2]+dx_inv[np-3]);
      a[np-1] = SQR(dx_inv[np-2]+dx_inv[np-3]);
      d[np-1] = dy[np-2]*(2*dx_inv[np-2]+3*dx_inv[np-3])*SQR(dx_inv[np-2]) + dy[np-3]*CUBE(dx_inv[np-3]);
    }
    break;
  }

  c[0] /= b[0];
  d[0] /= b[0];

  for(unsigned ip = 1; ip<np; ip++) {
    double id = 1.0/(b[ip] - c[ip-1]*a[ip]);
    c[ip] *= id;
    d[ip] = (d[ip] - d[ip-1] * a[ip])*id;
  }

  m[np-1] = d[np-1];
  for(int ip=np-2; ip>=0; ip--) {
    m[ip] = d[ip]-c[ip]*m[ip+1];
  }

  return m;

  // std::vector<double> ypp(np);
  // for(unsigned ip=0; ip<np-1; ip++) {
  //   ypp[ip] = 6*dy[ip]*SQR(dx_inv[ip]) - (4*m[ip]+2*m[ip+1])*dx_inv[ip];
  // }
  // ypp[np-1] = -6*dy[np-2]*SQR(dx_inv[np-2]) + (2*m[np-2]+4*m[np-1])*dx_inv[np-2];
  //
  // return ypp;
}

Eigen::VectorXd calin::math::spline_interpolation::generate_cubic_spline_interpolation_eigen(
  const Eigen::VectorXd& x, const Eigen::VectorXd& y,
  BoundaryConitions bc_lhs, double bc_lhs_val,
  BoundaryConitions bc_rhs, double bc_rhs_val)
{
  return std_to_eigenvec(
    generate_cubic_spline_interpolation(eigen_to_stdvec(x), eigen_to_stdvec(y),
      bc_lhs, bc_lhs_val, bc_rhs, bc_rhs_val));
}

std::vector<double> calin::math::spline_interpolation::
fit_spline(
  const std::vector<double>& xknots,
  const std::vector<double>& xdata,
  const std::vector<double>& ydata,
  const std::vector<double>& wdata,
  BoundaryConitions bc_lhs, double bc_lhs_val,
  BoundaryConitions bc_rhs, double bc_rhs_val,
  double lambda, RegularizationType reg_type)
{
  unsigned N = xknots.size();
  unsigned M = xdata.size();
  if (ydata.size() != M) {
    throw std::runtime_error("fit_spline: xdata and ydata must have the same size");
  }
  if (!wdata.empty() && wdata.size() != M) {
    throw std::runtime_error("fit_spline: wdata must have the same size as xdata");
  }
  if (N < 3) {
    throw std::runtime_error("fit_spline: need at least 3 knots");
  }

  unsigned P = 0;
  if (lambda > 0.0) {
    if (reg_type == REG_FIRST_DIFFERENCE && N >= 2) {
      P = N - 1;
    } else if (reg_type == REG_SECOND_DIFFERENCE && N >= 3) {
      P = N - 2;
    }
  }

  Eigen::MatrixXd A(M + P, N);
  Eigen::VectorXd b(M + P);

  A.setZero();
  b.setZero();

  CubicMultiSpline multi_basis(xknots);

  std::vector<double> zero_y(N, 0.0);
  unsigned idx_inhom = multi_basis.add_spline(zero_y, "inhom", bc_lhs, bc_lhs_val, bc_rhs, bc_rhs_val);

  std::vector<unsigned> basis_idx(N);
  for (unsigned i = 0; i < N; ++i) {
    std::vector<double> basis_y(N, 0.0);
    basis_y[i] = 1.0;
    basis_idx[i] = multi_basis.add_spline(basis_y, "", bc_lhs, 0.0, bc_rhs, 0.0);
  }

  for (unsigned j = 0; j < M; ++j) {
    double x = xdata[j];
    std::vector<double> vals = multi_basis.value(x);
    b(j) = ydata[j] - vals[idx_inhom];
    for (unsigned i = 0; i < N; ++i) {
      A(j, i) = vals[basis_idx[i]];
    }
  }

  if (!wdata.empty()) {
    for (unsigned j = 0; j < M; ++j) {
      double w = std::sqrt(wdata[j]);
      A.row(j) *= w;
      b(j) *= w;
    }
  }

  if (lambda > 0.0 && P > 0) {
    double sqrt_lambda = std::sqrt(lambda);
    if (reg_type == REG_FIRST_DIFFERENCE) {
      for (unsigned i = 0; i < N - 1; ++i) {
        A(M + i, i)     = -sqrt_lambda;
        A(M + i, i + 1) =  sqrt_lambda;
      }
    } else if (reg_type == REG_SECOND_DIFFERENCE) {
      for (unsigned i = 0; i < N - 2; ++i) {
        A(M + i, i)     =  sqrt_lambda;
        A(M + i, i + 1) = -2.0 * sqrt_lambda;
        A(M + i, i + 2) =  sqrt_lambda;
      }
    }
  }

  Eigen::VectorXd y = A.colPivHouseholderQr().solve(b);

  std::vector<double> result(N);
  for (unsigned i = 0; i < N; ++i) {
    result[i] = y(i);
  }
  return result;
}

void calin::math::spline_interpolation::
fit_regular_spline(
  std::vector<double>& xknots_out,
  std::vector<double>& yknots_out,
  unsigned nknot,
  const std::vector<double>& xdata,
  const std::vector<double>& ydata,
  const std::vector<double>& wdata,
  BoundaryConitions bc_lhs, double bc_lhs_val,
  BoundaryConitions bc_rhs, double bc_rhs_val,
  double lambda, RegularizationType reg_type)
{
  if (xdata.empty()) {
    throw std::runtime_error("fit_regular_spline: xdata is empty");
  }
  if (nknot < 3) {
    throw std::runtime_error("fit_regular_spline: need at least 3 knots");
  }

  auto minmax = std::minmax_element(xdata.begin(), xdata.end());
  double xmin = *minmax.first;
  double xmax = *minmax.second;

  xknots_out.resize(nknot);
  double dx = (xmax - xmin) / (nknot - 1);
  for (unsigned i = 0; i < nknot; ++i) {
    xknots_out[i] = xmin + i * dx;
  }
  xknots_out.back() = xmax;

  yknots_out = fit_spline(xknots_out, xdata, ydata, wdata, bc_lhs, bc_lhs_val, bc_rhs, bc_rhs_val, lambda, reg_type);
  return;
}

double calin::math::spline_interpolation::
cubic_solve(double y, double dx, double dx_inv, double y0, double y1, double D0, double D1)
{
  // Roots of cubic from Numerical Recipies

  using calin::math::special::CUBE;
  using calin::math::special::SQR;

  D0 *= dx;
  D1 *= dx;
  double dy = y1-y0;

  double f = 1/(-2*dy + (D0 + D1));

  double c = f * (y0-y);
  double b = f * D0;
  double a = f * (3*dy - (2*D0 + D1));

  double Q1 = (SQR(a) - 3*b)/9;
  double R1 = (2*CUBE(a) - 9*a*b + 27*c)/54;

  double Q3 = CUBE(Q1);
  double R2 = SQR(R1);

  double t;
  if(R2 < Q3) {
    double theta = std::acos(R1/std::sqrt(Q3));
    t = -2*std::sqrt(Q1)*std::cos((theta-2*M_PI)/3)-a/3;
    if(t<0 or t>1) {
      t = -2*std::sqrt(Q1)*std::cos(theta/3)-a/3;
      if(t<0 or t>1) {
        t = -2*std::sqrt(Q1)*std::cos((theta+2*M_PI)/3)-a/3;
      }
    }
  } else {
    double A = -((R1<0)?-1.0:1.0) * std::cbrt(std::abs(R1) + std::sqrt(R2-Q3));
    double B = (A==0)?0:Q1/A;
    t = (A+B)-a/3;
  }

  return t;
}

CubicSpline::
CubicSpline(const std::vector<double>& x, const std::vector<double>& y,
    BoundaryConitions bc_lhs, double bc_lhs_val,
    BoundaryConitions bc_rhs, double bc_rhs_val):
  s_(), I_(), bc_lhs_(bc_lhs), bc_lhs_val_(bc_lhs_val), bc_rhs_(bc_rhs), bc_rhs_val_(bc_rhs_val)
{
  *static_cast<InterpolationIntervals*>(&s_) = make_intervals(x);
  s_.y = y;
  s_.dy_dx = generate_cubic_spline_interpolation(x, y, bc_lhs, bc_lhs_val, bc_rhs, bc_rhs_val);
  init();
}

CubicSpline::
CubicSpline(const std::vector<double>& x, const std::vector<double>& y,
    const std::vector<double>& dy_dx):
  s_(), I_()
{
  *static_cast<InterpolationIntervals*>(&s_) = make_intervals(x);
  s_.y = y;
  s_.dy_dx = dy_dx;
  init();
  bc_lhs_ = BC_CLAMPED_SLOPE;
  bc_lhs_val_ = derivative(xmin());
  bc_rhs_ = BC_CLAMPED_SLOPE;
  bc_rhs_val_ = derivative(xmax());
}

void CubicSpline::init()
{
  double I=0;
  I_.push_back(I);
  for(unsigned i=0; i<s_.x.size()-1; i++) {
    double dx = s_.x[i+1]-s_.x[i];
    double dx_inv = 1.0/dx;
    I += cubic_integral(1.0, dx, dx_inv, s_.y[i], s_.y[i+1], s_.dy_dx[i], s_.dy_dx[i+1]);
    I_.push_back(I);
  }
  for(unsigned iy=1;iy<s_.y.size();iy++) {
    if(s_.y[iy]<s_.y[iy-1])y_is_monotonic_inc_ = false;
    if(s_.y[iy]>s_.y[iy-1])y_is_monotonic_dec_ = false;
  }
}

CubicSpline* CubicSpline::new_regularized_spline(double dx) const
{
  std::vector<double> x_knots;
  if(dx <= 0.0) {
    throw std::out_of_range("CubicSpline::new_regularized_spline: dx must be positive");
  }
  for(double x=s_.xmin; x<s_.xmax+dx*0.01; x+=dx) {
    x_knots.push_back(x);
  }

  std::vector<double> y_knots(x_knots.size());
  std::transform(x_knots.begin(), x_knots.end(), y_knots.begin(),
    [this](double x){ return this->value(x); });

  // CubicSpline* new_spline = new CubicSpline(x_knots, y_knots,
  //   bc_lhs_, bc_lhs_val_, bc_rhs_, bc_rhs_val_);

  CubicSpline* new_spline = new CubicSpline(x_knots, y_knots,
    BC_CLAMPED_SLOPE, this->derivative(x_knots.front()),
    BC_CLAMPED_SLOPE, this->derivative(x_knots.back()));

  return new_spline;
}

void CubicSpline::rescale(double scale)
{
  for(auto& y : s_.y) { y *= scale; }
  for(auto& dy_dx : s_.dy_dx) { dy_dx *= scale; }
}

void CubicSpline::extend_linear_rhs(double dx)
{
  if(dx<=0) {
    dx = s_.dx.back();
  }
  if(dx==s_.regular_dx and s_.regular_xmax==s_.xmax) {
    s_.regular_xmax += dx;
    s_.irregular_begin += 1;
  }
  
  double new_y = s_.y.back() + s_.dy_dx.back()*dx;
  double dI = dx * (s_.y.back() + new_y) / 2.0;

  s_.xmax += dx;
  s_.x.push_back(s_.x.back() + dx);
  s_.dx.push_back(dx);
  s_.dx_inv.push_back(1/dx);
  s_.y.push_back(new_y);
  s_.dy_dx.push_back(s_.dy_dx.back());

  if (!I_.empty()) {
    I_.push_back(I_.back() + dI);
  }
}

void CubicSpline::extend_linear_lhs(double dx)
{
  if(dx<=0) {
    dx = s_.dx.front();
  }
  
  double new_y = s_.y.front() - s_.dy_dx.front()*dx;
  double dI = dx * (s_.y.front() + new_y) / 2.0;

  if (s_.irregular_begin > 0 && std::abs(dx - s_.regular_dx)/s_.regular_dx < 1e-6) {
    s_.irregular_begin += 1;
  } else if (s_.irregular_begin == 0 && std::abs(dx - s_.dx.front())/s_.dx.front() < 1e-6) {
    s_.regular_dx = dx;
    s_.regular_dx_inv = 1.0/dx;
    s_.regular_xmax = s_.x[1];
    s_.irregular_begin = 2;
  } else {
    s_.regular_dx = 0;
    s_.regular_dx_inv = 0;
    s_.regular_xmax = s_.xmin - dx;
    s_.irregular_begin = 0;
  }

  s_.xmin -= dx;
  s_.x.insert(s_.x.begin(), s_.xmin);
  s_.dx.insert(s_.dx.begin(), dx);
  s_.dx_inv.insert(s_.dx_inv.begin(), 1/dx);
  s_.y.insert(s_.y.begin(), new_y);
  s_.dy_dx.insert(s_.dy_dx.begin(), s_.dy_dx.front());

  if (!I_.empty()) {
    for (auto& integral : I_) {
      integral += dI;
    }
    I_.insert(I_.begin(), 0.0);
  }
}

void CubicSpline::make_monotonic()
{
  for (unsigned i = 0; i < s_.y.size() - 1; ++i) {
    double m = (s_.y[i+1] - s_.y[i]) / s_.dx[i];
    if (m == 0.0) {
      s_.dy_dx[i] = 0.0;
      s_.dy_dx[i+1] = 0.0;
    } else {
      double alpha = s_.dy_dx[i] / m;
      double beta = s_.dy_dx[i+1] / m;
      if (alpha < 0.0) { s_.dy_dx[i] = 0.0; alpha = 0.0; }
      if (beta < 0.0) { s_.dy_dx[i+1] = 0.0; beta = 0.0; }
      double norm = alpha * alpha + beta * beta;
      if (norm > 9.0) {
        double tau = 3.0 / std::sqrt(norm);
        s_.dy_dx[i] = tau * alpha * m;
        s_.dy_dx[i+1] = tau * beta * m;
      }
    }
  }
  I_.clear();
  init();
}

double CubicSpline::ymax() const
{
  double ymax = s_.y.front();
  for(unsigned iknot=0; iknot<s_.y.size()-1; ++iknot) {
    double dx = s_.dx[iknot];
    double dx_inv = 1.0/dx;
    double y0 = s_.y[iknot];
    double y1 = s_.y[iknot + 1];
    double D0 = s_.dy_dx[iknot];
    double D1 = s_.dy_dx[iknot + 1];
    double t0;
    double t1;
    cubic_extrema(t0, t1, dx, dx_inv, y0, y1, D0, D1);
    if(not std::isnan(t0) and t0>0 and t0 <1) {
      ymax = std::max(ymax, cubic_value(t0, dx, dx_inv, y0, y1, D0, D1));
    }
    if(not std::isnan(t1) and t1>0 and t1<1) {
      ymax = std::max(ymax, cubic_value(t1, dx, dx_inv, y0, y1, D0, D1));
    }
    ymax = std::max(ymax, y1);
  }
  return ymax;
}

double CubicSpline::ymin() const
{
  double ymin = s_.y.front();
  for(unsigned iknot=0; iknot<s_.y.size()-1; ++iknot) {
    double dx = s_.dx[iknot];
    double dx_inv = 1.0/dx;
    double y0 = s_.y[iknot];
    double y1 = s_.y[iknot + 1];
    double D0 = s_.dy_dx[iknot];
    double D1 = s_.dy_dx[iknot + 1];
    double t0;
    double t1;
    cubic_extrema(t0, t1, dx, dx_inv, y0, y1, D0, D1);
    if(not std::isnan(t0) and t0>0 and t0 <1) {
      ymin = std::min(ymin, cubic_value(t0, dx, dx_inv, y0, y1, D0, D1));
    }
    if(not std::isnan(t1) and t1>0 and t1<1) {
      ymin = std::min(ymin, cubic_value(t1, dx, dx_inv, y0, y1, D0, D1));
    }
    ymin = std::min(ymin, y1);
  }
  return ymin;
}

double CubicSpline::x_at_ymin() const
{
  double xmin = s_.x.front();
  double ymin = s_.y.front();
  for(unsigned iknot=0; iknot<s_.y.size()-1; ++iknot) {
    double dx = s_.dx[iknot];
    double dx_inv = 1.0/dx;
    double y0 = s_.y[iknot];
    double y1 = s_.y[iknot + 1];
    double D0 = s_.dy_dx[iknot];
    double D1 = s_.dy_dx[iknot + 1];
    double t0;
    double t1;
    cubic_extrema(t0, t1, dx, dx_inv, y0, y1, D0, D1);
    if(not std::isnan(t0) and t0>0 and t0<1) {
      double yval = cubic_value(t0, dx, dx_inv, y0, y1, D0, D1);
      if(yval < ymin) {
        ymin = yval;
        xmin = t0 * dx + s_.x[iknot];
      }
    }
    if(not std::isnan(t1) and t1>0 and t1<1) {
      double yval = cubic_value(t1, dx, dx_inv, y0, y1, D0, D1);
      if(yval < ymin) {
        ymin = yval;
        xmin = t1 * dx + s_.x[iknot];
      }
    }
    if(y1 < ymin) {
      ymin = y1;
      xmin = s_.x[iknot+1];
    }
  }
  return xmin;
}

double CubicSpline::x_at_ymax() const
{
  double xmax = s_.x.front();
  double ymax = s_.y.front();
  for(unsigned iknot=0; iknot<s_.y.size()-1; ++iknot) {
    double dx = s_.dx[iknot];
    double dx_inv = 1.0/dx;
    double y0 = s_.y[iknot];
    double y1 = s_.y[iknot + 1];
    double D0 = s_.dy_dx[iknot];
    double D1 = s_.dy_dx[iknot + 1];
    double t0;
    double t1;
    cubic_extrema(t0, t1, dx, dx_inv, y0, y1, D0, D1);
    if(not std::isnan(t0) and t0>0 and t0<1) {
      double yval = cubic_value(t0, dx, dx_inv, y0, y1, D0, D1);
      if(yval > ymax) {
        ymax = yval;
        xmax = t0 * dx + s_.x[iknot];
      }
    }
    if(not std::isnan(t1) and t1>0 and t1<1) {
      double yval = cubic_value(t1, dx, dx_inv, y0, y1, D0, D1);
      if(yval > ymax) {
        ymax = yval;
        xmax = t1 * dx + s_.x[iknot];
      }
    }
    if(y1 > ymax) {
      ymax = y1;
      xmax = s_.x[iknot+1];
    }
  }
  return xmax;
}

double CubicSpline::value(double x) const
{
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  return cubic_value(t, dx, dx_inv, s_.y[i], s_.y[i+1], s_.dy_dx[i], s_.dy_dx[i+1]);
}

double CubicSpline::value_with_linear_extrapolation(double x) const
{
  if (x < s_.xmin) {
    return s_.y.front() + s_.dy_dx.front() * (x - s_.xmin);
  }
  if (x > s_.xmax) {
    return s_.y.back() + s_.dy_dx.back() * (x - s_.xmax);
  }
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  return cubic_value(t, dx, dx_inv, s_.y[i], s_.y[i+1], s_.dy_dx[i], s_.dy_dx[i+1]);
}

double CubicSpline::derivative(double x) const
{
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  return cubic_1st_derivative(t, dx, dx_inv, s_.y[i], s_.y[i+1], s_.dy_dx[i], s_.dy_dx[i+1]);
}

double CubicSpline::derivative_and_value(double x, double& value) const
{
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  return cubic_1st_derivative_and_value(value, t, dx, dx_inv, s_.y[i], s_.y[i+1], s_.dy_dx[i], s_.dy_dx[i+1]);
}

double CubicSpline::second_derivative(double x) const
{
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  return cubic_2nd_derivative(t, dx, dx_inv, s_.y[i], s_.y[i+1], s_.dy_dx[i], s_.dy_dx[i+1]);
}

double CubicSpline::third_derivative(double x) const
{
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  return cubic_3rd_derivative(t, dx, dx_inv, s_.y[i], s_.y[i+1], s_.dy_dx[i], s_.dy_dx[i+1]);
}

double CubicSpline::integral(double x) const
{
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  return cubic_integral(t, dx, dx_inv, s_.y[i], s_.y[i+1], s_.dy_dx[i], s_.dy_dx[i+1], I_[i]);
}

double CubicSpline::invert(double y) const
{
  std::vector<double>::const_iterator ifind = s_.y.end();
  if(y_is_monotonic_inc_) {
    ifind = std::upper_bound(s_.y.begin(), s_.y.end(), y);
  } else if(y_is_monotonic_dec_) {
    ifind = std::upper_bound(s_.y.begin(), s_.y.end(), y, std::greater<double>());
  } else {
    throw std::runtime_error("invert: y values must be monotonic");
  }
  if(ifind == s_.y.end())return NAN;
  if(ifind == s_.y.begin())return NAN;
  unsigned i = ifind - s_.y.begin() - 1;
  double dx = s_.x[i+1]-s_.x[i];
  double dx_inv = 1.0/dx;
  double t;
  if(y == s_.y[i+1]) {
    t = 1.0;
  } else if(y == s_.y[i]) {
    t = 0.0;
  } else {
    t = cubic_solve(y, dx, dx_inv, s_.y[i], s_.y[i+1], s_.dy_dx[i], s_.dy_dx[i+1]);
  }
  //std::cout << i << ' ' << t << '\n';
  return t*dx + s_.x[i];
}

double CubicSpline::find(double y, double xmin) const
{
  double xfound = std::numeric_limits<double>::infinity();
  for(unsigned i = find_interval(xmin, s_); i<s_.x.size()-1 and
    xfound == std::numeric_limits<double>::infinity(); ++i)
  {
    double x0 = s_.x[i];
    double x1 = s_.x[i+1];
    double dx = x1-x0;
    double y0 = s_.y[i];
    double y1 = s_.y[i+1];
    double dy = y1-y0;
    double D0 = s_.dy_dx[i]*dx;
    double D1 = s_.dy_dx[i+1]*dx;

    double c = 3*dy - (2*D0 + D1);
    double d = -2*dy + (D0 + D1);
    double d_inv = 1/d;

    double A = c * d_inv;
    double B = D0 * d_inv;
    double C = (y0-y) * d_inv;

    double t0,t1,t2;
    unsigned nroot = calin::math::special::solve_cubic_equation(t0,t1,t2,A,B,C);

    if(t0>=0 and t0<=1 and t0*dx+x0>xmin) {
      xfound = std::min(xfound, t0*dx+x0);
    }
    if(nroot>1 and t1>=0 and t1<=1 and t1*dx+x0>xmin) {
      xfound = std::min(xfound, t1*dx+x0);
    }
    if(nroot>2 and t2>=0 and t2<=1 and t2*dx+x0>xmin) {
      xfound = std::min(xfound, t2*dx+x0);
    }
  }
  return xfound;
}

Eigen::VectorXd CubicSpline::turning_points(bool include_extrapolation) const
{
  std::vector<double> tps;
  for (unsigned i = 0; i < s_.y.size() - 1; ++i) {
    double dx = s_.dx[i];
    double y0 = s_.y[i];
    double y1 = s_.y[i+1];
    double D0 = s_.dy_dx[i] * dx;
    double D1 = s_.dy_dx[i+1] * dx;
    double dy = y1 - y0;
    double c = 3*dy - (2*D0 + D1);
    double d = -2*dy + (D0 + D1);

    std::vector<double> roots;
    if (std::abs(d) < 1e-14 * std::max({std::abs(dy), std::abs(D0), std::abs(D1), 1.0})) {
      if (std::abs(c) > 1e-14 * std::max({std::abs(D0), std::abs(D1), 1.0})) {
        roots.push_back(-D0 / (2*c));
      }
    } else {
      double disc = c*c - 3*D0*d;
      if (disc >= 0.0) {
        double sq = std::sqrt(disc);
        double q = -c - std::copysign(sq, c);
        roots.push_back(q / (3*d));
        if (q != 0.0) {
          roots.push_back(D0 / q);
        }
      }
    }

    for (double t : roots) {
      bool keep = false;
      if (t >= 0.0 && t <= 1.0) {
        keep = true;
      } else if (include_extrapolation) {
        if (i == 0 && t < 0.0) keep = true;
        if (i == s_.y.size() - 2 && t > 1.0) keep = true;
      }
      if (keep) {
        tps.push_back(s_.x[i] + t * dx);
      }
    }
  }

  std::sort(tps.begin(), tps.end());
  
  if (!tps.empty()) {
    std::vector<double> unique_tps;
    unique_tps.push_back(tps[0]);
    for (size_t i = 1; i < tps.size(); ++i) {
      if (std::abs(tps[i] - unique_tps.back()) > 1e-9) {
        unique_tps.push_back(tps[i]);
      }
    }
    tps = std::move(unique_tps);
  }

  Eigen::VectorXd result(tps.size());
  for (size_t i = 0; i < tps.size(); ++i) {
    result(i) = tps[i];
  }
  return result;
}

CubicMultiSpline::
CubicMultiSpline(const std::vector<double>& x)
{
  s_ = make_intervals(x);
}

unsigned CubicMultiSpline::add_spline(const std::vector<double>& y,
  const std::string& name,
  BoundaryConitions bc_lhs, double bc_lhs_val,
  BoundaryConitions bc_rhs, double bc_rhs_val)
{
  unsigned ispline = y_.size();
  y_.emplace_back(y);
  dy_dx_.emplace_back(generate_cubic_spline_interpolation(s_.x, y, bc_lhs, bc_lhs_val, bc_rhs, bc_rhs_val));
  info_.emplace_back(name, bc_lhs, bc_lhs_val, bc_rhs, bc_rhs_val);
  return ispline;
}

void CubicMultiSpline::replace_spline(unsigned ispline,
  const std::vector<double>& y, const std::string& name,
  BoundaryConitions bc_lhs, double bc_lhs_val,
  BoundaryConitions bc_rhs, double bc_rhs_val)
{
  if(ispline >= y_.size()) {
    throw std::out_of_range("replace_spline: ispline out of range");
  }
  y_[ispline] = y;
  dy_dx_[ispline] = generate_cubic_spline_interpolation(s_.x, y, bc_lhs, bc_lhs_val, bc_rhs, bc_rhs_val);
  if(not name.empty()) { info_[ispline].name = name; }
  info_[ispline].bc_lhs = bc_lhs;
  info_[ispline].bc_lhs_val = bc_lhs_val;
  info_[ispline].bc_rhs = bc_rhs;
  info_[ispline].bc_rhs_val = bc_rhs_val;
}

CubicMultiSpline* CubicMultiSpline::new_regularized_multi_spline(double dx) const
{
  std::vector<double> x_knots;
  if(dx <= 0.0) {
    dx = s_.regular_dx;
  }
  for(double x=s_.xmin; x<s_.xmax; x+=dx) {
    x_knots.push_back(x);
  }
  if(x_knots.back() < s_.xmax) {
    x_knots.push_back(s_.xmax);
  }
  CubicMultiSpline* new_spline = new CubicMultiSpline(x_knots);
  for(unsigned ispline=0; ispline<info_.size(); ++ispline) {
    std::vector<double> y_knots(x_knots.size());
    std::transform(x_knots.begin(), x_knots.end(), y_knots.begin(),
      [this,ispline](double x){ return this->value(x, ispline); });
    new_spline->add_spline(y_knots, info_[ispline].name,
      info_[ispline].bc_lhs, info_[ispline].bc_lhs_val,
      info_[ispline].bc_rhs, info_[ispline].bc_rhs_val);
  }
  return new_spline;
}

double CubicMultiSpline::value(double x, unsigned ispline) const
{
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  return cubic_value(t, dx, dx_inv, y_[ispline][i], y_[ispline][i+1],
    dy_dx_[ispline][i], dy_dx_[ispline][i+1]);
}

void CubicMultiSpline::value(double x, unsigned ispline0, double& value0,
  unsigned ispline1, double& value1) const
{
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  value0 = cubic_value(t, dx, dx_inv, y_[ispline0][i], y_[ispline0][i+1],
    dy_dx_[ispline0][i], dy_dx_[ispline0][i+1]);
  value1 = cubic_value(t, dx, dx_inv, y_[ispline1][i], y_[ispline1][i+1],
    dy_dx_[ispline1][i], dy_dx_[ispline1][i+1]);
}

void CubicMultiSpline::value(double x, unsigned ispline0, double& value0,
  unsigned ispline1, double& value1, unsigned ispline2, double& value2) const
{
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  value0 = cubic_value(t, dx, dx_inv, y_[ispline0][i], y_[ispline0][i+1],
    dy_dx_[ispline0][i], dy_dx_[ispline0][i+1]);
  value1 = cubic_value(t, dx, dx_inv, y_[ispline1][i], y_[ispline1][i+1],
    dy_dx_[ispline1][i], dy_dx_[ispline1][i+1]);
  value2 = cubic_value(t, dx, dx_inv, y_[ispline2][i], y_[ispline2][i+1],
    dy_dx_[ispline2][i], dy_dx_[ispline2][i+1]);
}

void CubicMultiSpline::value(double x,
  unsigned ispline0, double& value0, unsigned ispline1, double& value1,
  unsigned ispline2, double& value2, unsigned ispline3, double& value3) const
{
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  value0 = cubic_value(t, dx, dx_inv, y_[ispline0][i], y_[ispline0][i+1],
    dy_dx_[ispline0][i], dy_dx_[ispline0][i+1]);
  value1 = cubic_value(t, dx, dx_inv, y_[ispline1][i], y_[ispline1][i+1],
    dy_dx_[ispline1][i], dy_dx_[ispline1][i+1]);
  value2 = cubic_value(t, dx, dx_inv, y_[ispline2][i], y_[ispline2][i+1],
    dy_dx_[ispline2][i], dy_dx_[ispline2][i+1]);
  value3 = cubic_value(t, dx, dx_inv, y_[ispline3][i], y_[ispline3][i+1],
    dy_dx_[ispline3][i], dy_dx_[ispline3][i+1]);
}

void CubicMultiSpline::value(double x,
  unsigned ispline0, double& value0, unsigned ispline1, double& value1,
  unsigned ispline2, double& value2, unsigned ispline3, double& value3,
  unsigned ispline4, double& value4) const
{
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  value0 = cubic_value(t, dx, dx_inv, y_[ispline0][i], y_[ispline0][i+1],
    dy_dx_[ispline0][i], dy_dx_[ispline0][i+1]);
  value1 = cubic_value(t, dx, dx_inv, y_[ispline1][i], y_[ispline1][i+1],
    dy_dx_[ispline1][i], dy_dx_[ispline1][i+1]);
  value2 = cubic_value(t, dx, dx_inv, y_[ispline2][i], y_[ispline2][i+1],
    dy_dx_[ispline2][i], dy_dx_[ispline2][i+1]);
  value3 = cubic_value(t, dx, dx_inv, y_[ispline3][i], y_[ispline3][i+1],
    dy_dx_[ispline3][i], dy_dx_[ispline3][i+1]);
  value4 = cubic_value(t, dx, dx_inv, y_[ispline4][i], y_[ispline4][i+1],
    dy_dx_[ispline4][i], dy_dx_[ispline4][i+1]);
}

void CubicMultiSpline::value(double x,
  unsigned ispline0, double& value0, unsigned ispline1, double& value1,
  unsigned ispline2, double& value2, unsigned ispline3, double& value3,
  unsigned ispline4, double& value4, unsigned ispline5, double& value5) const
{
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  value0 = cubic_value(t, dx, dx_inv, y_[ispline0][i], y_[ispline0][i+1],
    dy_dx_[ispline0][i], dy_dx_[ispline0][i+1]);
  value1 = cubic_value(t, dx, dx_inv, y_[ispline1][i], y_[ispline1][i+1],
    dy_dx_[ispline1][i], dy_dx_[ispline1][i+1]);
  value2 = cubic_value(t, dx, dx_inv, y_[ispline2][i], y_[ispline2][i+1],
    dy_dx_[ispline2][i], dy_dx_[ispline2][i+1]);
  value3 = cubic_value(t, dx, dx_inv, y_[ispline3][i], y_[ispline3][i+1],
    dy_dx_[ispline3][i], dy_dx_[ispline3][i+1]);
  value4 = cubic_value(t, dx, dx_inv, y_[ispline4][i], y_[ispline4][i+1],
    dy_dx_[ispline4][i], dy_dx_[ispline4][i+1]);
  value5 = cubic_value(t, dx, dx_inv, y_[ispline5][i], y_[ispline5][i+1],
    dy_dx_[ispline5][i], dy_dx_[ispline5][i+1]);
}

std::vector<double> CubicMultiSpline::value(double x) const
{
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  std::vector<double> values(y_.size());
  for(unsigned ispline=0;ispline<y_.size();ispline++) {
    values[ispline] =
      cubic_value(t, dx, dx_inv, y_[ispline][i], y_[ispline][i+1],
        dy_dx_[ispline][i], dy_dx_[ispline][i+1]);
  }
  return values;
}

Eigen::VectorXd CubicMultiSpline::value_as_eigen(double x) const
{
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  Eigen::VectorXd values(y_.size());
  for(unsigned ispline=0;ispline<y_.size();ispline++) {
    values(ispline) =
      cubic_value(t, dx, dx_inv, y_[ispline][i], y_[ispline][i+1],
        dy_dx_[ispline][i], dy_dx_[ispline][i+1]);
  }
  return values;
}

std::vector<double> CubicMultiSpline::derivative(double x) const
{
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  std::vector<double> derivative(y_.size());
  for(unsigned ispline=0;ispline<y_.size();ispline++) {
    derivative[ispline] =
      cubic_1st_derivative(t, dx, dx_inv, y_[ispline][i], y_[ispline][i+1],
        dy_dx_[ispline][i], dy_dx_[ispline][i+1]);
  }
  return derivative;
}

double CubicMultiSpline::derivative(double x, unsigned ispline) const
{
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  return cubic_1st_derivative(t, dx, dx_inv, y_[ispline][i], y_[ispline][i+1],
    dy_dx_[ispline][i], dy_dx_[ispline][i+1]);
}

double CubicMultiSpline::derivative_and_value(double x, unsigned ispline, double& value) const
{
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  return cubic_1st_derivative_and_value(value, t, dx, dx_inv, y_[ispline][i], y_[ispline][i+1],
    dy_dx_[ispline][i], dy_dx_[ispline][i+1]);
}

double CubicMultiSpline::second_derivative(double x, unsigned ispline) const
{
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  return cubic_2nd_derivative(t, dx, dx_inv, y_[ispline][i], y_[ispline][i+1],
    dy_dx_[ispline][i], dy_dx_[ispline][i+1]);
}

double CubicMultiSpline::second_derivative_and_value(double x, unsigned ispline,
  double& first_derivative, double& value) const
{
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  return cubic_2nd_derivative_and_value(first_derivative, value,
    t, dx, dx_inv, y_[ispline][i], y_[ispline][i+1],
    dy_dx_[ispline][i], dy_dx_[ispline][i+1]);
}

double CubicMultiSpline::integral(double x, unsigned ispline) const
{
  double dx;
  double dx_inv;
  unsigned i = find_interval(x, s_, dx, dx_inv);
  double t = (x-s_.x[i])*dx_inv;
  double I0 = 0;
  for(unsigned j=0;j<i;j++) {
    I0 = cubic_integral(1.0, dx, dx_inv, y_[ispline][j], y_[ispline][j+1],
      dy_dx_[ispline][j], dy_dx_[ispline][j+1], I0);
  }
  return cubic_integral(t, dx, dx_inv, y_[ispline][i], y_[ispline][i+1],
    dy_dx_[ispline][i], dy_dx_[ispline][i+1], I0);
}

double CubicMultiSpline::test_vcl_value(double x, unsigned ispline, unsigned ivec) const
{
  typedef typename calin::util::vcl::VCL256Architecture::double_vt double_vt;
  double_vt vy = vcl_value<calin::util::vcl::VCL256Architecture>(x, ispline);
  return vy[ivec];
}

std::vector<double> CubicMultiSpline::test_vcl_value(double x0, double x1, double x2, double x3,
  unsigned ispline) const
{
  typedef typename calin::util::vcl::VCL256Architecture::double_vt double_vt;
  std::vector<double> ay(4);
  double_vt vx(x0,x1,x2,x3);
  double_vt vy = vcl_value<calin::util::vcl::VCL256Architecture>(vx, ispline);
  vy.store(ay.data());
  return ay;
}
