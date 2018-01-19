/*

   calin/iact_data/waveform_treatment_event_visitor.hpp -- Stephen Fegan -- 2018-01-11

   Waveform treatment event data visitor - process waveforms

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

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

inline void calin::iact_data::waveform_treatment_event_visitor::
SingleGainDualWindowWaveformTreatmentEventVisitor::
analyze_waveforms(
  const uint16_t*__restrict__ data, unsigned nchan, int nsamp,
  int window_n, int bkg_window_0, const int* sig_window_0,
  float*__restrict__ ped, float ped_iir_old, float ped_iir_new,
  int*__restrict__ chan_max_index, int*__restrict__ chan_max,
  int*__restrict__ chan_bkg_win_sum, int* chan_sig_win_sum,
  int*__restrict__ chan_sig_max_sum, int*__restrict__ chan_sig_max_sum_index,
  int*__restrict__ chan_all_sum_q, int*__restrict__ chan_all_sum_qt,
  float*__restrict__ chan_sig, float*__restrict__ chan_mean_t)
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
  int samp[nsamp];

  for(unsigned ichan=0;ichan<nchan;ichan++)
  {
    samp[0] = data[ichan*nsamp];
    imax = 0;
    max = samp[0];
    win = max;
    sum_qt = 0;
    for(isamp = 1;isamp<window_n;isamp++) {
      const int _samp = data[ichan*nsamp+isamp];
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
    while(isamp<nsamp) {
      int iss = isamp-16;
      if(bkg_window_0 == iss)bkg = win;
      if(sig_window_0[ichan] == iss)sig = win;
      const int _samp = data[ichan*nsamp+isamp];
      samp[isamp] = _samp;
      sum_q += _samp;
      sum_qt += _samp*isamp;
      win += _samp - samp[iss];
      if(win>sig_max) {
        sig_max = win;
        isig_max = iss;
      }
      if(_samp > max) {
        imax = isamp;
        max = _samp;
      }
      ++isamp;
    }
    if(bkg_window_0 == nsamp-window_n)bkg = win;
    if(sig_window_0[ichan] == nsamp-window_n)sig = win;

    chan_max_index[ichan] = imax;
    chan_max[ichan] = max;
    chan_bkg_win_sum[ichan] = bkg;
    chan_sig_win_sum[ichan] = sig;
    chan_sig_max_sum[ichan] = sig_max;
    chan_sig_max_sum_index[ichan] = isig_max;
    chan_all_sum_q[ichan] = sum_q;
    chan_all_sum_qt[ichan] = sum_qt;
    if(ped[ichan] >= 0) {
      ped[ichan] = ped_iir_old*ped[ichan] + ped_iir_new*float(bkg);
    } else {
      ped[ichan] = float(bkg);
    }
    chan_sig[ichan] = float(sig) - ped[ichan];
    chan_mean_t[ichan] =
      (double(window_n*sum_qt) - double(ped[ichan]*nsamp*(nsamp-1)/2))/
        (double(window_n*sum_q) - double(ped[ichan]*nsamp));
  }
}
