/*

   calin/iact_data/waveform_treatment_event_visitor.cpp -- Stephen Fegan -- 2018-01-11

   Waveform treatment event data visitor - process waveforms

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
#include <algorithm>

#include <util/memory.hpp>
#include <iact_data/waveform_treatment_event_visitor.hpp>

using namespace calin::iact_data::waveform_treatment_event_visitor;
using calin::util::memory::safe_aligned_recalloc;

void calin::iact_data::waveform_treatment_event_visitor::
OptimalWindowSumWaveformTreatmentEventVisitor::
scalar_analyze_waveforms(const uint16_t*__restrict__ data)
{
  int imax = 0;
  int max = 0;
  int bkg = 0;
  int sig = 0;
  int opt = 0;
  int opt_qt = 0;
  int iopt = 0;
  int all = 0;
  int win = 0;
  int win_qt = 0;
  unsigned isamp = 0;
  int samp[nsamp_];

  for(unsigned ichan=0;ichan<nchan_;ichan++)
  {
    samp[0] = data[ichan*nsamp_];
    imax = 0;
    max = samp[0];
    win = max;
    win_qt = 0;
    for(isamp = 1;isamp<window_n_;isamp++) {
      const int _samp = data[ichan*nsamp_+isamp];
      samp[isamp] = _samp;
      win += _samp;
      if(_samp > max) {
        imax = isamp;
        max = _samp;
      }
      win_qt += _samp*isamp;
    }
    opt = win;
    opt_qt = win_qt;
    iopt = 0;
    all = win;
    while(isamp<nsamp_) {
      int iss = isamp-window_n_;
      if(bkg_window_0_ == iss)bkg = win;
      if(sig_window_0_[ichan] == iss)sig = win;
      const int _samp = data[ichan*nsamp_+isamp];
      samp[isamp] = _samp;
      all += _samp;
      win_qt += _samp*isamp - samp[iss]*iss;
      win += _samp - samp[iss];
      if(win>opt) {
        opt = win;
        opt_qt = win_qt;
        iopt = iss+1;
      }
      if(_samp > max) {
        imax = isamp;
        max = _samp;
      }
      ++isamp;
    }
    if(bkg_window_0_ == int(nsamp_-window_n_))bkg = win;
    if(sig_window_0_[ichan] == int(nsamp_-window_n_))sig = win;

    chan_max_index_[ichan] = imax;
    chan_max_[ichan] = max;
    chan_bkg_win_sum_[ichan] = bkg;
    chan_sig_win_sum_[ichan] = sig;
    chan_opt_win_sum_[ichan] = opt;
    chan_opt_win_sum_qt_[ichan] = opt_qt;
    chan_opt_win_index_[ichan] = iopt;
    chan_all_sum_[ichan] = all;
  }
}

// *****************************************************************************
// *****************************************************************************
//
// POSSIBLY OBSOLETE VERSIONS
//
// *****************************************************************************
// *****************************************************************************

void calin::iact_data::waveform_treatment_event_visitor::
SingleGainDualWindowWaveformTreatmentEventVisitor::
scalar_analyze_waveforms(const uint16_t*__restrict__ data)
{
  int imax = 0;
  int max = 0;
  int bkg = 0;
  int sig = 0;
  int sig_max = 0;
  int isig_max = 0;
  int sum_q = 0;
  int sum_qt = 0;
  int win = 0;
  unsigned isamp = 0;
  int samp[nsamp_];

  for(unsigned ichan=0;ichan<nchan_;ichan++)
  {
    samp[0] = data[ichan*nsamp_];
    imax = 0;
    max = samp[0];
    win = max;
    sum_qt = 0;
    for(isamp = 1;isamp<window_n_;isamp++) {
      const int _samp = data[ichan*nsamp_+isamp];
      samp[isamp] = _samp;
      win += _samp;
      if(_samp > max) {
        imax = isamp;
        max = _samp;
      }
      sum_qt += _samp*isamp;
    }
    sig_max = win;
    isig_max = 0;
    sum_q = win;
    while(isamp<nsamp_) {
      int iss = isamp-window_n_;
      if(bkg_window_0_ == iss)bkg = win;
      if(sig_window_0_[ichan] == iss)sig = win;
      const int _samp = data[ichan*nsamp_+isamp];
      samp[isamp] = _samp;
      sum_q += _samp;
      sum_qt += _samp*isamp;
      win += _samp - samp[iss];
      if(win>sig_max) {
        sig_max = win;
        isig_max = iss+1;
      }
      if(_samp > max) {
        imax = isamp;
        max = _samp;
      }
      ++isamp;
    }
    if(bkg_window_0_ == int(nsamp_-window_n_))bkg = win;
    if(sig_window_0_[ichan] == int(nsamp_-window_n_))sig = win;

    chan_max_index_[ichan] = imax;
    chan_max_[ichan] = max;
    chan_bkg_win_sum_[ichan] = bkg;
    chan_sig_win_sum_[ichan] = sig;
    chan_sig_max_sum_[ichan] = sig_max;
    chan_sig_max_sum_index_[ichan] = isig_max;
    chan_all_sum_q_[ichan] = sum_q;
    chan_all_sum_qt_[ichan] = sum_qt;
    if(chan_ped_est_[ichan] >= 0) {
      chan_ped_est_[ichan] = ped_iir_old_*chan_ped_est_[ichan] + ped_iir_new_*float(bkg);
    } else {
      chan_ped_est_[ichan] = float(bkg);
    }
    chan_sig_[ichan] = float(sig) - chan_ped_est_[ichan];
    chan_mean_t_[ichan] =
      (double(window_n_*sum_qt) - double(chan_ped_est_[ichan]*nsamp_*(nsamp_-1)/2))/
        (double(window_n_*sum_q) - double(chan_ped_est_[ichan]*nsamp_));
  }
}
