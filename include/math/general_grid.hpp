/*

   calin/math/general_grid.hpp -- Stephen Fegan -- 2017-08-28

   Collection of functions which operate on general grids

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <Eigen/Core>

namespace calin { namespace math { namespace general_grid {

class Grid
{
public:
  virtual ~Grid();
  virtual unsigned num_neighbours() = 0;
  virtual void gridid_to_xy(unsigned gridid, double& x, double& y) = 0;
  virtual unsigned xy_to_gridid(double x, double y) = 0;
  virtual std::vector<unsigned> gridid_to_neighbour_gridids(unsigned gridid) = 0;
  virtual void gridid_to_vertexes_xy(unsigned gridid,
    Eigen::VectorXd& xv, Eigen::VectorXd& yv) = 0;
};

class HexGrid: public Grid
{
public:
  HexGrid(double scale = 1, double theta = 0, double xoffset = 0, double yoffset = 0):
    scale_(scale), ctheta_(std::cos(theta)), stheta_(std::sin(theta)),
    xoffset_(xoffset), yoffset_(yoffset) { /* nothing to see here */ }
  virtual ~HexGrid();
  unsigned num_neighbours() override;
  void gridid_to_xy(unsigned gridid, double& x, double& y) override;
  unsigned xy_to_gridid(double x, double y) override;
  std::vector<unsigned> gridid_to_neighbour_gridids(unsigned gridid) override;
  void gridid_to_vertexes_xy(unsigned gridid,
    Eigen::VectorXd& xv, Eigen::VectorXd& yv) override;
private:
  double scale_ = 1.0;
  double ctheta_ = 1.0;
  double stheta_ = 0.0;
  double xoffset_ = 0.0;
  double yoffset_ = 0.0;
};

std::vector<std::pair<unsigned,unsigned> > compute_region_boundary(
  const std::vector<unsigned>& region_gridids, Grid* grid);

unsigned num_bounday_curves(
  const std::vector<std::pair<unsigned,unsigned> >& boundary);

void extract_bounday_curve(
  Eigen::VectorXd& xv, Eigen::VectorXd& yv, Grid* grid,
  const std::vector<std::pair<unsigned,unsigned> >& boundary,
  unsigned icurve = 0, bool closed = false);

} } } // namespace calin::math::general_grid
