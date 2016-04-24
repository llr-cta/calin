/*

   calin/math/covariance_calc.cpp -- Stephen Fegan -- 2016-04-24

   Utility functions for covariance calculation

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

#include <cstdlib>

#include <math/covariance_calc.hpp>

double calin::math::covariance_calc::cov_i64_gen(int64_t sij, int64_t nij,
  int64_t si, int64_t ni, int64_t sj, int64_t nj)
{
  auto div_ij = std::div(sij, nij);
  const int64_t qij = div_ij.quot;
  const int64_t rij = div_ij.rem;

  auto div_i = std::div(si, ni);
  const int64_t qi = div_i.quot;
  const int64_t ri = div_i.rem;

  auto div_j = std::div(sj, nj);
  const int64_t qj = div_j.quot;
  const int64_t rj = div_j.rem;

  auto div_pi = std::div(ri*qj, ni);
  const int64_t qpi = div_pi.quot;
  const int64_t rpi = div_pi.rem;

  auto div_pj = std::div(rj*qi, nj);
  const int64_t qpj = div_pj.quot;
  const int64_t rpj = div_pj.rem;

  double val = double(qij-qi*qj-qpi-qpj)
    + double(rij)/double(nij) - double(ri*rj)/double(ni*nj)
    - double(rpi)/double(ni) - double(rpj)/double(nj);

  return val;
}
