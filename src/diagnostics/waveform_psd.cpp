/*

   calin/diagnostics/waveform.cpp -- Stephen Fegan -- 2016-04-13

   Waveform diagnostics visitor

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

#include <math/special.hpp>
#include <math/fftw_util.hpp>
#include <io/log.hpp>
#include <diagnostics/waveform.hpp>

using calin::math::special::SQR;
using namespace calin::math::fftw_util;
using namespace calin::io::log;
using namespace calin::diagnostics::waveform;

WaveformPSDVisitor::WaveformPSDVisitor():
  TelescopeEventVisitor()
{
  // nothing to see here
}

WaveformPSDVisitor::~WaveformPSDVisitor()
{
  // nothing to see here
}

bool WaveformPSDVisitor::demand_waveforms()
{
  return true;
}

bool WaveformPSDVisitor::is_parallelizable()
{
  return true;
}

WaveformPSDVisitor* WaveformPSDVisitor::new_sub_visitor()
{
  auto* sub_visitor = new WaveformPSDVisitor;
  sub_visitor->parent_ = this;
  return sub_visitor;
}

bool WaveformPSDVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* run_config)
{
  run_config_ = run_config;
  results_.Clear();
  unsigned nsample = run_config->num_samples();
  fftw_data_ = fftw_alloc_real(nsample);
  assert(fftw_data_);
  fftw_plan_fwd_ =
    fftw_plan_r2r_1d(nsample, fftw_data_, fftw_data_, FFTW_R2HC, 0);
  assert(fftw_plan_fwd_);
  fftw_plan_bwd_ =
    fftw_plan_r2r_1d(nsample, fftw_data_, fftw_data_, FFTW_HC2R, 0);
  assert(fftw_plan_bwd_);
  unsigned nfreq = 1 + nsample/2; // 1->1 2->2 3->2 4->3 5->3 etc
  for(int ichan = 0; ichan<run_config->configured_channel_id_size(); ichan++)
  {
    auto* hg_wf = results_.add_high_gain();
    hg_wf->mutable_psd_sum()->Resize(nfreq,0);
    hg_wf->mutable_psd_sum_squared()->Resize(nfreq,0);
    hg_wf->mutable_corr_sum()->Resize(nsample,0);
    hg_wf->mutable_corr_sum_squared()->Resize(nsample,0);
    auto* lg_wf = results_.add_low_gain();
    lg_wf->mutable_psd_sum()->Resize(nfreq,0);
    lg_wf->mutable_psd_sum_squared()->Resize(nfreq,0);
    lg_wf->mutable_corr_sum()->Resize(nsample,0);
    lg_wf->mutable_corr_sum_squared()->Resize(nsample,0);
  }
  return true;
}

bool WaveformPSDVisitor::leave_telescope_run()
{
  run_config_ = nullptr;
  fftw_destroy_plan(fftw_plan_bwd_);
  fftw_plan_bwd_ = nullptr;
  fftw_destroy_plan(fftw_plan_fwd_);
  fftw_plan_fwd_ = nullptr;
  fftw_free(fftw_data_);
  fftw_data_ = nullptr;
  return true;
}

bool WaveformPSDVisitor::visit_telescope_event(
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  // nothing to see here
  return true;
}

bool WaveformPSDVisitor::visit_waveform(unsigned ichan,
  calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
  calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain)
{
  const int index = run_config_->configured_channel_index(ichan);
  if(high_gain)
    process_one_waveform(high_gain, results_.mutable_high_gain(index));
  if(low_gain)
    process_one_waveform(low_gain, results_.mutable_low_gain(index));
  return true;
}

void WaveformPSDVisitor::
process_one_waveform(
  const calin::ix::iact_data::telescope_event::ChannelWaveform* wf,
  ix::diagnostics::waveform::WaveformRawPSD* psd)
{
  const unsigned nsample = run_config_->num_samples();
  assert(wf->samples_size() == int(nsample));
  psd->set_num_entries(psd->num_entries()+1);
  const auto* sample = wf->samples().data();
  auto* psd_sum = psd->mutable_psd_sum()->mutable_data();
  auto* psd_sum_squared = psd->mutable_psd_sum_squared()->mutable_data();
  std::transform(sample, sample+nsample, fftw_data_, [](decltype(*sample) x) {
    return double(x); });
  fftw_execute(fftw_plan_fwd_);
  const double* ri = fftw_data_;
  const double* ci = fftw_data_ + nsample-1;
  double psdi = SQR(*ri++);
  (*psd_sum++) += psdi;
  (*psd_sum_squared++) += SQR(psdi);
  while(ri < ci)
  {
    psdi = SQR(*ri++) + SQR(*ci--);
    (*psd_sum++) += psdi;
    (*psd_sum_squared++) += SQR(psdi);
  }
  if(ri==ci)
  {
    double psdi = SQR(*ri);
    *psd_sum += psdi;
    *psd_sum_squared += SQR(psdi);
  }
  hcvec_scale_and_multiply_conj(fftw_data_, fftw_data_, fftw_data_, nsample);
  fftw_execute(fftw_plan_bwd_);
  auto* corr_sum = psd->mutable_corr_sum()->mutable_data();
  for(unsigned isample=0; isample<nsample; isample++)
    corr_sum[isample] += fftw_data_[isample];
  auto* corr_sum_squared = psd->mutable_corr_sum_squared()->mutable_data();
  for(unsigned isample=0; isample<nsample; isample++)
    corr_sum_squared[isample] += fftw_data_[isample] * fftw_data_[isample];
}

bool WaveformPSDVisitor::merge_results()
{
  if(parent_)
  {
    assert(results_.high_gain_size() == parent_->results_.high_gain_size());
    assert(results_.low_gain_size() == parent_->results_.low_gain_size());
    for(int ichan=0; ichan<results_.high_gain_size(); ichan++)
      merge_one_gain(&results_.high_gain(ichan),
        parent_->results_.mutable_high_gain(ichan));
    for(int ichan=0; ichan<results_.low_gain_size(); ichan++)
      merge_one_gain(&results_.low_gain(ichan),
        parent_->results_.mutable_low_gain(ichan));
  }
  return true;
}

void WaveformPSDVisitor::merge_one_gain(
  const ix::diagnostics::waveform::WaveformRawPSD* from,
  ix::diagnostics::waveform::WaveformRawPSD* to)
{
  assert(to->psd_sum_size() == from->psd_sum_size());
  assert(to->psd_sum_squared_size() == from->psd_sum_squared_size());
  to->set_num_entries(to->num_entries() + from->num_entries());
  for(int i=0; i<from->psd_sum_size(); i++)
    to->set_psd_sum(i,to->psd_sum(i) + from->psd_sum(i));
  for(int i=0; i<from->psd_sum_squared_size(); i++)
    to->set_psd_sum_squared(i,to->psd_sum_squared(i) + from->psd_sum_squared(i));
  for(int i=0; i<from->corr_sum_size(); i++)
    to->set_corr_sum(i,to->corr_sum(i) + from->corr_sum(i));
  for(int i=0; i<from->corr_sum_squared_size(); i++)
    to->set_corr_sum_squared(i,to->corr_sum_squared(i) + from->corr_sum_squared(i));
}

Eigen::VectorXd WaveformPSDVisitor::psd_mean(
  const ix::diagnostics::waveform::WaveformRawPSD* stat)
{
  const int N = stat->psd_sum_size();
  Eigen::VectorXd m(N);
  const double one_over_n = 1.0/double(stat->num_entries());
  for(int i=0; i<N; i++)
    m(i) = stat->psd_sum(i) * one_over_n;
  return m;
}

Eigen::VectorXd WaveformPSDVisitor::psd_var(
  const ix::diagnostics::waveform::WaveformRawPSD* stat)
{
  const int N = stat->psd_sum_size();
  Eigen::VectorXd v(N);
  const double one_over_n = 1.0/double(stat->num_entries());
  for(int i=0; i<N; i++)
    v(i) = double(stat->psd_sum_squared(i)) * one_over_n
      - SQR(stat->psd_sum(i) * one_over_n);
  return v;
}

Eigen::VectorXd WaveformPSDVisitor::corr_mean(
  const ix::diagnostics::waveform::WaveformRawPSD* psd_stat,
  const ix::diagnostics::waveform::WaveformRawStats* trace_stat)
{
  const int N = psd_stat->corr_sum_size();
  assert(trace_stat->sum_size() == N);

  Eigen::VectorXd m = WaveformStatsVisitor::waveform_mean(trace_stat);

  auto* fftw_data = fftw_alloc_real(N);
  assert(fftw_data);
  auto fftw_plan_fwd =
    fftw_plan_r2r_1d(N, fftw_data, fftw_data, FFTW_R2HC, 0);
  assert(fftw_plan_fwd);
  auto fftw_plan_bwd =
    fftw_plan_r2r_1d(N, fftw_data, fftw_data, FFTW_HC2R, 0);
  assert(fftw_plan_bwd);

  std::copy(m.data(), m.data() + N, fftw_data);
  fftw_execute(fftw_plan_fwd);
  hcvec_scale_and_multiply_conj(fftw_data, fftw_data, fftw_data, N);
  fftw_execute(fftw_plan_bwd);

  const double one_over_n = 1.0/double(trace_stat->num_entries());
  const double dbl_N = double(N);
  for(int i=0; i<N; i++)
    m(i) = (psd_stat->corr_sum(i)*one_over_n - fftw_data[i])/dbl_N;

  fftw_destroy_plan(fftw_plan_bwd);
  fftw_destroy_plan(fftw_plan_fwd);
  fftw_free(fftw_data);

  return m;
}

Eigen::VectorXd WaveformPSDVisitor::corr_mean_centered(
  const ix::diagnostics::waveform::WaveformRawPSD* psd_stat,
  const ix::diagnostics::waveform::WaveformRawStats* trace_stat,
  Eigen::VectorXd& h)
{
  Eigen::VectorXd mp = WaveformPSDVisitor::corr_mean(psd_stat, trace_stat);
  int N = mp.size();
  Eigen::VectorXd m(N);
  h.resize(N);
  int iz = (N-1)/2;
  for(int i=iz; i<N; i++)h(i) = i-iz, m(i) = mp(i-iz);
  for(int i=0;i<iz; i++)h(i) = i-iz, m(i) = mp(N+i-iz);
  return m;
}

Eigen::VectorXd WaveformPSDVisitor::corr_var(
  const ix::diagnostics::waveform::WaveformRawPSD* psd_stat,
  const ix::diagnostics::waveform::WaveformRawStats* trace_stat)
{
  const int N = psd_stat->corr_sum_size();
  Eigen::VectorXd v(N);
  const double one_over_n = 1.0/double(psd_stat->num_entries());
  for(int i=0; i<N; i++)
    v(i) = double(psd_stat->corr_sum_squared(i)) * one_over_n
      - SQR(psd_stat->corr_sum(i) * one_over_n);
  return v;
}
