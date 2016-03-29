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
#include <io/log.hpp>
#include <diagnostics/waveform.hpp>

using calin::math::special::SQR;
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
  fftw_plan_ = fftw_plan_r2r_1d(nsample, fftw_data_, fftw_data_, FFTW_R2HC, 0);
  assert(fftw_plan_);
  unsigned nfreq = 1 + nsample/2; // 1->1 2->2 3->2 4->3 5->3 etc
  for(int ichan = 0; ichan<run_config->configured_channel_id_size(); ichan++)
  {
    auto* hg_wf = results_.add_high_gain();
    hg_wf->mutable_sum()->Resize(nfreq,0);
    hg_wf->mutable_sum_squared()->Resize(nfreq,0);
    auto* lg_wf = results_.add_low_gain();
    lg_wf->mutable_sum()->Resize(nfreq,0);
    lg_wf->mutable_sum_squared()->Resize(nfreq,0);
  }
  return true;
}

bool WaveformPSDVisitor::leave_telescope_run()
{
  run_config_ = nullptr;
  fftw_destroy_plan(fftw_plan_);
  fftw_plan_ = nullptr;
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
  auto* sum = psd->mutable_sum()->mutable_data();
  auto* sum_squared = psd->mutable_sum_squared()->mutable_data();
  std::transform(sample, sample+nsample, fftw_data_, [](decltype(*sample) x) {
    return double(x); });
  fftw_execute(fftw_plan_);
  const double* ri = fftw_data_;
  const double* ci = fftw_data_ + nsample-1;
  double psdi = SQR(*ri++);
  (*sum++) += psdi;
  (*sum_squared++) += SQR(psdi);
  while(ri < ci)
  {
    psdi = SQR(*ri++) + SQR(*ci--);
    (*sum++) += psdi;
    (*sum_squared++) += SQR(psdi);
  }
  if(ri==ci)
  {
    double psdi = SQR(*ri);
    *sum += psdi;
    *sum_squared += SQR(psdi);
  }
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
  assert(to->sum_size() == from->sum_size());
  assert(to->sum_squared_size() == from->sum_squared_size());
  to->set_num_entries(to->num_entries() + from->num_entries());
  for(int i=0; i<from->sum_size(); i++)
    to->set_sum(i,to->sum(i) + from->sum(i));
  for(int i=0; i<from->sum_squared_size(); i++)
    to->set_sum_squared(i,to->sum_squared(i) + from->sum_squared(i));
}

Eigen::VectorXd WaveformPSDVisitor::psd_mean(
  const ix::diagnostics::waveform::WaveformRawPSD* stat)
{
  const int N = stat->sum_size();
  Eigen::VectorXd m(N);
  const double one_over_n = 1.0/double(stat->num_entries());
  for(int i=0; i<N; i++)
    m(i) = double(stat->sum(i)) * one_over_n;
  return m;
}

Eigen::VectorXd WaveformPSDVisitor::psd_var(
  const ix::diagnostics::waveform::WaveformRawPSD* stat)
{
  const int N = stat->sum_size();
  Eigen::VectorXd v(N);
  const double one_over_n = 1.0/double(stat->num_entries());
  for(int i=0; i<N; i++)
    v(i) = double(stat->sum_squared(i)) * one_over_n
      - SQR(double(stat->sum(i)) * one_over_n);
  return v;
}
