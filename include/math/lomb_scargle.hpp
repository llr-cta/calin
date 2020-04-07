/*

   calin/math/lomb_scargle.hpp -- Stephen Fegan -- 2017-12-05

   Calculate Lomb-Scargle periodogram

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#pragma once

#include <Eigen/Dense>

namespace calin { namespace math { namespace lomb_scargle {

std::pair<double, double> amplitudes(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti, double freq);

double power(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti, double freq);

Eigen::VectorXd periodogram(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti,
  double freq_lo, double freq_hi, double delta_freq);

Eigen::VectorXd frequencies(const Eigen::VectorXd& periodogram,
  double freq_lo, double delta_freq);

// -----------------------------------------------------------------------------
// =============================================================================
// -----------------------------------------------------------------------------
//
// Classic "O(N * M)" algorithms with various optimizations
//
// -----------------------------------------------------------------------------
// =============================================================================
// -----------------------------------------------------------------------------

Eigen::VectorXd periodogram_slow(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti,
  double freq_lo, double freq_hi, double delta_freq);

Eigen::VectorXd periodogram_fast(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti,
  double freq_lo, double freq_hi, double delta_freq, unsigned renormalize_nfreq = 0);

Eigen::VectorXd periodogram_vcl128(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti,
  double freq_lo, double freq_hi, double delta_freq, unsigned renormalize_nfreq = 0,
  unsigned unroll = 2);

Eigen::VectorXd periodogram_vcl256(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti,
  double freq_lo, double freq_hi, double delta_freq, unsigned renormalize_nfreq = 0,
  unsigned unroll = 2);

Eigen::VectorXd periodogram_vcl256_float(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti,
  double freq_lo, double freq_hi, double delta_freq, unsigned renormalize_nfreq = 0,
  unsigned unroll = 2);

Eigen::MatrixXd multi_periodogram_vcl128(const Eigen::MatrixXd& xi, const Eigen::VectorXd& ti,
  double freq_lo, double freq_hi, double delta_freq, unsigned renormalize_nfreq = 0,
  unsigned unroll = 4);

Eigen::MatrixXd multi_periodogram_vcl256(const Eigen::MatrixXd& xi, const Eigen::VectorXd& ti,
  double freq_lo, double freq_hi, double delta_freq, unsigned renormalize_nfreq = 0,
  unsigned unroll = 4);

Eigen::MatrixXd multi_periodogram_vcl256_float(const Eigen::MatrixXd& xi, const Eigen::VectorXd& ti,
  double freq_lo, double freq_hi, double delta_freq, unsigned renormalize_nfreq = 0,
  unsigned unroll = 4);

} } } // namespace calin::math::lomb_scargle
