/*

   calin/math/geometry.cpp -- Stephen Fegan -- 2021-04-30

   Functions for geometrical manuipulations

   Copyright 2021, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <numeric>
#include <cmath>
#include <random>

#include <math/geometry.hpp>
#include <util/log.hpp>

using namespace calin::util::log;
using namespace calin::math::geometry;

Eigen::Quaterniond calin::math::geometry::
euler_to_quaternion(const calin::ix::common_types::EulerAngles3D& euler)
{
  using Eigen::Quaterniond;
  using Eigen::AngleAxisd;
  using Eigen::Vector3d;
  Eigen::Quaterniond q;
  switch(euler.rotation_order()) {
  case calin::ix::common_types::EulerAngles3D::ZXZ:
    return AngleAxisd(euler.alpha()*(M_PI/180.0), Vector3d::UnitZ())
          * AngleAxisd(euler.beta()*(M_PI/180.0), Vector3d::UnitX())
          * AngleAxisd(euler.gamma()*(M_PI/180.0), Vector3d::UnitZ());
  case calin::ix::common_types::EulerAngles3D::XYX:
    return AngleAxisd(euler.alpha()*(M_PI/180.0), Vector3d::UnitX())
          * AngleAxisd(euler.beta()*(M_PI/180.0), Vector3d::UnitY())
          * AngleAxisd(euler.gamma()*(M_PI/180.0), Vector3d::UnitX());
  case calin::ix::common_types::EulerAngles3D::YZY:
    return AngleAxisd(euler.alpha()*(M_PI/180.0), Vector3d::UnitY())
          * AngleAxisd(euler.beta()*(M_PI/180.0), Vector3d::UnitZ())
          * AngleAxisd(euler.gamma()*(M_PI/180.0), Vector3d::UnitY());
  case calin::ix::common_types::EulerAngles3D::ZYZ:
    return AngleAxisd(euler.alpha()*(M_PI/180.0), Vector3d::UnitZ())
          * AngleAxisd(euler.beta()*(M_PI/180.0), Vector3d::UnitY())
          * AngleAxisd(euler.gamma()*(M_PI/180.0), Vector3d::UnitZ());
  case calin::ix::common_types::EulerAngles3D::XZX:
    return AngleAxisd(euler.alpha()*(M_PI/180.0), Vector3d::UnitX())
          * AngleAxisd(euler.beta()*(M_PI/180.0), Vector3d::UnitZ())
          * AngleAxisd(euler.gamma()*(M_PI/180.0), Vector3d::UnitX());
  case calin::ix::common_types::EulerAngles3D::YXY:
    return AngleAxisd(euler.alpha()*(M_PI/180.0), Vector3d::UnitY())
          * AngleAxisd(euler.beta()*(M_PI/180.0), Vector3d::UnitX())
          * AngleAxisd(euler.gamma()*(M_PI/180.0), Vector3d::UnitY());
  default:
    throw std::runtime_error("Unsupported rotation order");
  }
}

Eigen::Matrix3d calin::math::geometry::
euler_to_matrix(const calin::ix::common_types::EulerAngles3D& euler)
{
  return euler_to_quaternion(euler).toRotationMatrix();
}

void calin::math::geometry::quaternion_to_euler(
  calin::ix::common_types::EulerAngles3D* euler, const Eigen::Quaterniond& q)
{
  return matrix_to_euler(euler, q.toRotationMatrix());
}

void calin::math::geometry::matrix_to_euler(
  calin::ix::common_types::EulerAngles3D* euler, const Eigen::Matrix3d& m)
{
  Eigen::Vector3d v;
  switch(euler->rotation_order()) {
  case calin::ix::common_types::EulerAngles3D::ZXZ:
    v = m.eulerAngles(2,0,2); break;
  case calin::ix::common_types::EulerAngles3D::XYX:
    v = m.eulerAngles(0,1,0); break;
  case calin::ix::common_types::EulerAngles3D::YZY:
    v = m.eulerAngles(1,2,1); break;
  case calin::ix::common_types::EulerAngles3D::ZYZ:
    v = m.eulerAngles(2,1,2); break;
  case calin::ix::common_types::EulerAngles3D::XZX:
    v = m.eulerAngles(0,2,0); break;
  case calin::ix::common_types::EulerAngles3D::YXY:
    v = m.eulerAngles(1,0,1); break;
  default:
    throw std::runtime_error("Unsupported rotation order");
  }
  if(v(1) < 0) {
    v(0) = std::fmod(v(0) + 2*M_PI, 2*M_PI) - M_PI;
    v(1) = std::fabs(v(1));
    v(2) = std::fmod(v(2) + 2*M_PI, 2*M_PI) - M_PI;
  }
  euler->set_alpha(v(0) * 180.0/M_PI);
  euler->set_beta(v(1) * 180.0/M_PI);
  euler->set_gamma(v(2) * 180.0/M_PI);
}

void calin::math::geometry::scattering_euler(
  calin::ix::common_types::EulerAngles3D* euler, double dispersion, math::rng::RNG& rng,
  double twist_dispersion)
{
  double alpha = 0;
  double beta = 0;
  double gamma = 0;

  if(dispersion>0)
  {
    alpha = rng.uniform() * 360.0 - 180.0;
    beta = dispersion*std::sqrt(-2.0*std::log(rng.uniform()));
    gamma = -alpha;
  }

  if(twist_dispersion > 0) {
    gamma = std::fmod(std::fmod(gamma + rng.uniform()*twist_dispersion + 180, 360) + 360, 360) - 180;
  }

  euler->set_alpha(alpha);
  euler->set_beta(beta);
  euler->set_gamma(gamma);
}

calin::ix::common_types::EulerAngles3D
calin::math::geometry::scattering_euler(double dispersion, calin::math::rng::RNG& rng,
  double twist_dispersion, calin::ix::common_types::EulerAngles3D::RotationOrder rotation_order)
{
  calin::ix::common_types::EulerAngles3D euler;
  euler.set_rotation_order(rotation_order);
  scattering_euler(&euler, dispersion, rng, twist_dispersion);
  return euler;
}

bool calin::math::geometry::euler_is_zero(
  const calin::ix::common_types::EulerAngles3D& euler)
{
  return euler.alpha()==0 and euler.beta()==0 and euler.gamma()==0;
}

namespace {

  // Algorithm to find the smallest enclosing circle of circles in 2D from ChatGPT

  struct circle_t {
    Eigen::Vector2d c;
    double r;
  };

  inline bool contains(const circle_t& a,
                       const circle_t& b,
                       double eps = 1e-9) {
    return (a.c - b.c).norm() + b.r <= a.r + eps;
  }

  circle_t enclosing_from_one(const circle_t& a) {
    return a;
  }

  circle_t enclosing_from_two(const circle_t& a,
                              const circle_t& b) {
    Eigen::Vector2d d = b.c - a.c;
    double dist = d.norm();

    if (dist + a.r <= b.r + 1e-12) return b;
    if (dist + b.r <= a.r + 1e-12) return a;

    Eigen::Vector2d dir = d / dist;
    Eigen::Vector2d p1 = a.c - dir * a.r;
    Eigen::Vector2d p2 = b.c + dir * b.r;

    Eigen::Vector2d center = 0.5 * (p1 + p2);
    double radius = 0.5 * (p2 - p1).norm();

    return {center, radius};
  }

  circle_t enclosing_from_three(const circle_t& a,
                                const circle_t& b,
                                const circle_t& c) {
    // Try all 2-circle bases first (LP basis logic)
    circle_t ab = enclosing_from_two(a, b);
    if (contains(ab, c)) return ab;

    circle_t ac = enclosing_from_two(a, c);
    if (contains(ac, b)) return ac;

    circle_t bc = enclosing_from_two(b, c);
    if (contains(bc, a)) return bc;

    // Solve true nonlinear system:
    // |x - p1| + r1 = |x - p2| + r2 = |x - p3| + r3

    Eigen::Vector2d x = (a.c + b.c + c.c) / 3.0;

    for (int iter = 0; iter < 30; ++iter) {
      double d1 = (x - a.c).norm() + a.r;
      double d2 = (x - b.c).norm() + b.r;
      double d3 = (x - c.c).norm() + c.r;

      Eigen::Vector2d g1 = (x - a.c).normalized();
      Eigen::Vector2d g2 = (x - b.c).normalized();
      Eigen::Vector2d g3 = (x - c.c).normalized();

      Eigen::Vector2d f;
      f(0) = d1 - d2;
      f(1) = d1 - d3;

      Eigen::Matrix2d J;
      J.row(0) = (g1 - g2).transpose();
      J.row(1) = (g1 - g3).transpose();

      double det = J.determinant();
      if (std::abs(det) < 1e-14) {
        break;
      }

      Eigen::Vector2d delta = J.inverse() * f;
      x -= delta;

      if (delta.norm() < 1e-10) {
        break;
      }
    }

    double radius = (x - a.c).norm() + a.r;
    circle_t result = {x, radius};

    if (!contains(result, a) ||
        !contains(result, b) ||
        !contains(result, c)) {
      circle_t best = ab;
      if (ac.r < best.r) best = ac;
      if (bc.r < best.r) best = bc;
      return best;
    }

    return result;
  }

  circle_t solve_basis(const std::vector<circle_t>& basis) {
    if (basis.empty()) {
      return {{0.0, 0.0}, -1.0};
    }

    if (basis.size() == 1) {
      return enclosing_from_one(basis[0]);
    }

    if (basis.size() == 2) {
      return enclosing_from_two(basis[0],
                                basis[1]);
    }

    return enclosing_from_three(basis[0],
                                basis[1],
                                basis[2]);
  }

  circle_t lp_solve(std::vector<circle_t>& circles,
                    std::vector<circle_t>& basis,
                    int n) {
    if (n == 0 || basis.size() == 3) {
      return solve_basis(basis);
    }

    circle_t d = lp_solve(circles, basis, n - 1);

    if (d.r >= 0.0 &&
        contains(d, circles[n - 1])) {
      return d;
    }

    basis.push_back(circles[n - 1]);
    circle_t result = lp_solve(circles, basis, n - 1);
    basis.pop_back();

    return result;
  }

  circle_t smallest_enclosing_circle_internal(
    std::vector<circle_t>& circles) {

    std::shuffle(circles.begin(),
                 circles.end(),
                 std::mt19937{std::random_device{}()});

    std::vector<circle_t> basis;
    return lp_solve(circles,
                    basis,
                    static_cast<int>(circles.size()));
  }

  Eigen::VectorXd smallest_enclosing_circle_driver(
    const Eigen::MatrixXd& input) {

    const int n = static_cast<int>(input.rows());
    Eigen::Vector3d out;
    out.setZero();

    if (n == 0) {
      return out;
    }

    std::vector<circle_t> circles;
    circles.reserve(n);

    int max_idx = 0;
    double max_r = input(0, 2);

    for (int i = 0; i < n; ++i) {
      circle_t c;
      c.c = Eigen::Vector2d(input(i, 0),
                            input(i, 1));
      c.r = input(i, 2);
      circles.push_back(c);

      if (c.r > max_r) {
        max_r = c.r;
        max_idx = i;
      }
    }

    // Fast path: biggest circle already contains all
    const circle_t& biggest = circles[max_idx];
    bool all_contained = true;

    for (int i = 0; i < n; ++i) {
      if (!contains(biggest, circles[i])) {
        all_contained = false;
        break;
      }
    }

    circle_t result;
    if (all_contained) {
      result = biggest;
    } else {
      result = smallest_enclosing_circle_internal(circles);
    }

    out << result.c.x(),
           result.c.y(),
           result.r;

    return out;
  }

} // unnamed namespace

Eigen::Vector3d calin::math::geometry::
containing_circle_2d(const Eigen::MatrixXd& circles)
{
  if(circles.cols() != 3) {
    throw std::runtime_error("circles must have 3 columns");
  }
  return smallest_enclosing_circle_driver(circles);
}

Eigen::Vector4d calin::math::geometry::
containing_circle_sphere(const Eigen::MatrixXd& circles_on_sphere)
{
  // Thanks to Fergal Daly for telling me (in ~1999) that the stereographic 
  // projection of a circle on a sphere is a circle in 2D. Yay !
  if(circles_on_sphere.cols() != 4) {
    throw std::runtime_error("circles_on_sphere must have 4 columns");
  }
  Eigen::MatrixXd circles_2d(circles_on_sphere.rows(),3);
  for(int i=0; i<circles_on_sphere.rows(); ++i) {
    circles_2d.row(i) = project_circle_on_sphere_to_2d(circles_on_sphere.row(i));
  }
  Eigen::Vector3d containing_2d = containing_circle_2d(circles_2d);
  return deproject_circle_on_sphere_from_2d(containing_2d);
}
