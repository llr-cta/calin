/*

   calin/math/general_grid.cpp -- Stephen Fegan -- 2017-08-28

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

#include <set>

#include <calin_global_definitions.hpp>
#include <math/general_grid.hpp>
#include <math/hex_array.hpp>

using namespace calin::math::general_grid;

Grid::~Grid()
{
  // nothing to see here
}

HexGrid::~HexGrid()
{
  // nothing to see here
}

unsigned HexGrid::num_neighbours()
{
  return 6;
}

void HexGrid::gridid_to_xy(unsigned gridid, double& x, double& y)
{
  return calin::math::hex_array::hexid_to_xy_trans(gridid, x, y,
    ctheta_, stheta_, scale_, xoffset_, yoffset_);
}

unsigned HexGrid::xy_to_gridid(double x, double y)
{
  return calin::math::hex_array::xy_trans_to_hexid(x, y,
    ctheta_, stheta_, scale_, xoffset_, yoffset_);
}

std::vector<unsigned> HexGrid::gridid_to_neighbour_gridids(unsigned gridid)
{
  return calin::math::hex_array::hexid_to_neighbor_hexids(gridid);
}

void HexGrid::gridid_to_vertexes_xy(unsigned gridid,
  Eigen::VectorXd& xv, Eigen::VectorXd& yv)
{
  std::vector<double> sxv;
  std::vector<double> syv;
  calin::math::hex_array::hexid_to_vertexes_xy_trans(gridid, sxv, syv,
    ctheta_, stheta_, scale_, xoffset_, yoffset_);
  xv = std_to_eigenvec(sxv);
  yv = std_to_eigenvec(syv);
}

namespace {

  std::pair<unsigned,unsigned> make_segment_key(unsigned id_in, unsigned id_out)
  {
    return std::make_pair(id_in, id_out);
  }

  static const std::pair<unsigned,unsigned> endtag(0U,0U);
}

std::vector<std::pair<unsigned,unsigned> > calin::math::general_grid::
compute_region_boundary(const std::vector<unsigned>& region_gridids, Grid* grid)
{
  assert(dynamic_cast<HexGrid>(grid)); // TEMPORARY

  unsigned nnbr = grid->num_neighbours();

  // First, assemble list of all boundary segments - we search for pairs of
  // grid cells with one inside the camera and the other one not
  std::set<std::pair<unsigned, unsigned> > boundary_segments;
  for(auto gridid : region_gridids)
  {
    auto nbr = grid->gridid_to_neighbour_gridids(gridid);
    for(auto nbrgridid : nbr)
    {
      auto segment_key = make_segment_key(gridid,nbrgridid);
      auto reciprocal_segment_key = make_segment_key(nbrgridid,gridid);
      auto ibs = boundary_segments.find(reciprocal_segment_key);
      if(ibs == boundary_segments.end())boundary_segments.insert(segment_key);
      else boundary_segments.erase(ibs);
    }
  }

  // Second, build boundary chain
  std::vector<std::pair<unsigned,unsigned> > boundary_chain;
  while(not boundary_segments.empty())
  {
    auto bs = *boundary_segments.begin();
    boundary_segments.erase(boundary_segments.begin());
    if(not boundary_chain.empty())boundary_chain.push_back(endtag);
    boundary_chain.push_back(bs);
    while(1)
    {
      unsigned gridid_i = bs.first;
      unsigned gridid_o = bs.second;

      auto nbr = grid->gridid_to_neighbour_gridids(gridid_i);
      unsigned idx = std::find(nbr.begin(), nbr.end(), gridid_o)-nbr.begin();
      unsigned testid = nbr[(idx+1)%nnbr];
      auto test_bs = make_segment_key(gridid_i,testid);
      auto test_ibs = boundary_segments.find(test_bs);
      if(test_ibs != boundary_segments.end())
      {
        boundary_segments.erase(test_ibs);
        boundary_chain.push_back(test_bs);
        bs = test_bs;
        continue;
      }

      nbr = grid->gridid_to_neighbour_gridids(gridid_o);
      idx = std::find(nbr.begin(), nbr.end(), gridid_i)-nbr.begin();
      testid = nbr[(idx+nnbr-1)%nnbr];
      test_bs = make_segment_key(testid,gridid_o);
      test_ibs = boundary_segments.find(test_bs);
      if(test_ibs != boundary_segments.end())
      {
        boundary_segments.erase(test_ibs);
        boundary_chain.push_back(test_bs);
        bs = test_bs;
        continue;
      }

      // This algorithm only works for HEX grids.. would have to add code here
      // for "straight" segments

      break;
    }
  }

  return boundary_chain;
}

unsigned calin::math::general_grid::num_bounday_curves(
  const std::vector<std::pair<unsigned,unsigned> >& boundary)
{
  if(boundary.empty())return 0;
  unsigned ncurve = 0;
  auto istart = boundary.begin();
  while(1)
  {
    auto iend = find(istart, boundary.end(), endtag);
    if(istart != iend)ncurve++;
    istart = iend;
    if(istart != boundary.end())++istart;
    else break;
  }
  return ncurve;
}

void calin::math::general_grid::extract_bounday_curve(
  Eigen::VectorXd& xv, Eigen::VectorXd& yv, Grid* grid,
  const std::vector<std::pair<unsigned,unsigned> >& boundary,
  unsigned icurve, bool closed)
{
  unsigned nnbr = grid->num_neighbours();

  auto istart = boundary.begin();
  auto iend = find(istart, boundary.end(), endtag);

  while(icurve)
  {
    if(iend == boundary.end())break;
    if(istart != iend) --icurve;
    istart = iend;
    ++istart;
    iend = find(istart, boundary.end(), endtag);
  }

  if(icurve != 0 or istart == iend)
  {
    xv.resize(0);
    yv.resize(0);
    return;
  }

  if(closed) {
    xv.resize(iend-istart+1);
    yv.resize(iend-istart+1);
  } else {
    xv.resize(iend-istart);
    yv.resize(iend-istart);
  }

  unsigned iv = 0;

  unsigned nidx;
  Eigen::VectorXd x;
  Eigen::VectorXd y;

  while(istart != iend)
  {
    unsigned gridid_i = istart->first;
    unsigned gridid_o = istart->second;
    auto nbr = grid->gridid_to_neighbour_gridids(gridid_i);
    nidx = std::find(nbr.begin(), nbr.end(), gridid_o)-nbr.begin();
    grid->gridid_to_vertexes_xy(gridid_i, x, y);
    xv[iv] = x[(nidx+nnbr-1)%nnbr];
    yv[iv] = y[(nidx+nnbr-1)%nnbr];
    ++iv;
    ++istart;
  }
  if(closed)
  {
    xv[iv] = xv[0];
    yv[iv] = yv[0];
  }
}
