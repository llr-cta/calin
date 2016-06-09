/*

   calin/math/covariance_calc.hpp -- Stephen Fegan -- 2016-04-24

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

#pragma once

#include <cstdint>

namespace calin { namespace math { namespace covariance_calc {

double cov_i64_gen(int64_t sij, int64_t nij,
  int64_t si, int64_t ni, int64_t sj, int64_t nj);

double cov_double_gen(double sij, int64_t nij,
  double si, int64_t ni, double sj, int64_t nj);

inline int64_t cov_gen(int64_t sij, int64_t nij,
  int64_t si, int64_t ni, int64_t sj, int64_t nj) {
  return cov_i64_gen(sij, nij, si, ni, sj, nj); }

inline int64_t cov_gen(int64_t sij, int32_t nij,
  int64_t si, int32_t ni, int64_t sj, int32_t nj) {
  return cov_i64_gen(sij, nij, si, ni, sj, nj); }

inline int64_t cov_gen(int32_t sij, int32_t nij,
  int32_t si, int32_t ni, int32_t sj, int32_t nj) {
  return cov_i64_gen(sij, nij, si, ni, sj, nj); }

inline int64_t cov_gen(int64_t sij, uint64_t nij,
  int64_t si, uint64_t ni, int64_t sj, uint64_t nj) {
  return cov_i64_gen(sij, nij, si, ni, sj, nj); }

inline int64_t cov_gen(int64_t sij, uint32_t nij,
  int64_t si, uint32_t ni, int64_t sj, uint32_t nj) {
  return cov_i64_gen(sij, nij, si, ni, sj, nj); }

inline int64_t cov_gen(int32_t sij, uint32_t nij,
  int32_t si, uint32_t ni, int32_t sj, uint32_t nj) {
  return cov_i64_gen(sij, nij, si, ni, sj, nj); }

inline double cov_gen(double sij, int64_t nij,
  double si, int64_t ni, double sj, int64_t nj) {
  return cov_double_gen(sij, nij, si, ni, sj, nj); }

inline double cov_gen(double sij, int32_t nij,
  double si, int32_t ni, double sj, int32_t nj) {
  return cov_double_gen(sij, nij, si, ni, sj, nj); }

inline double cov_gen(double sij, uint64_t nij,
  double si, uint64_t ni, double sj, uint64_t nj) {
  return cov_double_gen(sij, nij, si, ni, sj, nj); }

inline double cov_gen(double sij, uint32_t nij,
  double si, uint32_t ni, double sj, uint32_t nj) {
  return cov_double_gen(sij, nij, si, ni, sj, nj); }

} } } // namespace calin::math::covariance_calc
