/*

   calin/diagnostics/waveform_integrated.cpp -- Stephen Fegan -- 2016-05-10

   Waveform diagnostics visitor - integrated stats with histo and
   cross-channel covariance

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

#include <io/log.hpp>
#include <math/special.hpp>
#include <diagnostics/waveform.hpp>

using calin::math::special::SQR;
using namespace calin::diagnostics::waveform;
using calin::iact_data::functional_event_visitor::
  DualValueInt32FunctionalTelescopeEventVisitor;
using calin::ix::diagnostics::waveform::
  FunctionalWaveformStatsVisitorConfig;

FunctionalWaveformStatsVisitor::
FunctionalWaveformStatsVisitor(
    DualValueInt32FunctionalTelescopeEventVisitor* value_supplier,
    const FunctionalWaveformStatsVisitorConfig& config):
  TelescopeEventVisitor(), value_supplier_(value_supplier), config_(config)
{
  // nothing to see here
}

FunctionalWaveformStatsVisitor::~FunctionalWaveformStatsVisitor()
{
  // nothing to see here
}

bool FunctionalWaveformStatsVisitor::demand_waveforms()
{
  return true;
}

bool FunctionalWaveformStatsVisitor::is_parallelizable()
{
  return true;
}

FunctionalWaveformStatsVisitor*
FunctionalWaveformStatsVisitor::new_sub_visitor(
  const std::map<TelescopeEventVisitor*,TelescopeEventVisitor*>&
    antecedent_visitors)
{
  auto i_sub_value_supplier = antecedent_visitors.find(value_supplier_);
  if(i_sub_value_supplier == antecedent_visitors.end())
    throw std::logic_error("No thread specific value suppler!");
  auto* sub_value_supplier =
    dynamic_cast<DualValueInt32FunctionalTelescopeEventVisitor*>(
      i_sub_value_supplier->second);
  auto* sub_visitor =
    new FunctionalWaveformStatsVisitor(sub_value_supplier, config_);
  sub_visitor->parent_ = this;
  return sub_visitor;
}

bool FunctionalWaveformStatsVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* run_config)
{
  run_config_ = run_config;

  results_.Clear();

  int nchan = run_config_->configured_channel_id_size();
  for(auto chan_id : run_config_->configured_channel_id())
    results_.add_channel_id(chan_id);

  high_gain_mask_.resize(nchan);
  high_gain_signal_.resize(nchan);
  auto* hg_results = results_.mutable_high_gain();
  hg_results->mutable_num_sum_entries()->Resize(nchan,0);
  hg_results->mutable_num_sum_product_entries()->Resize(nchan*(nchan-1)/2,0);
  hg_results->mutable_sum()->Resize(nchan,0);
  hg_results->mutable_sum_squared()->Resize(nchan,0);
  hg_results->mutable_sum_product()->Resize(nchan*(nchan-1)/2,0);

  low_gain_mask_.resize(nchan);
  low_gain_signal_.resize(nchan);
  auto* lg_results = results_.mutable_low_gain();
  lg_results->mutable_num_sum_entries()->Resize(nchan,0);
  lg_results->mutable_num_sum_product_entries()->Resize(nchan*(nchan-1)/2,0);
  lg_results->mutable_sum()->Resize(nchan,0);
  lg_results->mutable_sum_squared()->Resize(nchan,0);
  lg_results->mutable_sum_product()->Resize(nchan*(nchan-1)/2,0);

  return true;
}

bool FunctionalWaveformStatsVisitor::leave_telescope_run()
{
  run_config_ = nullptr;
  return true;
}

bool FunctionalWaveformStatsVisitor::visit_telescope_event(
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  for(auto& ichan : high_gain_mask_)ichan = 0;
  for(auto& ichan : low_gain_mask_)ichan = 0;
  return true;
}

bool FunctionalWaveformStatsVisitor::leave_telescope_event()
{
  if(results_.has_high_gain())
    process_one_gain(high_gain_mask_, high_gain_signal_,
      results_.mutable_high_gain());
  if(results_.has_low_gain())
    process_one_gain(low_gain_mask_, low_gain_signal_,
      results_.mutable_low_gain());
  return true;
}

bool FunctionalWaveformStatsVisitor::visit_waveform(unsigned ichan,
  calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
  calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain)
{
  const int index = run_config_->configured_channel_index(ichan);
  if(high_gain) {
    high_gain_mask_[index] = 1;
    high_gain_signal_[index] = value_supplier_->high_gain_value();
  }
  if(low_gain) {
    low_gain_mask_[index] = 1;
    low_gain_signal_[index] = value_supplier_->low_gain_value();
  }
  return true;
}

void FunctionalWaveformStatsVisitor::visit_one_waveform(
  const calin::ix::iact_data::telescope_event::ChannelWaveform* wf,
  unsigned index, std::vector<int>& mask, std::vector<int32_t>& signal)
{
  mask[index] = 1;
}

void FunctionalWaveformStatsVisitor::process_one_gain(
  const std::vector<int>& mask,
  const std::vector<int32_t>& signal,
  calin::ix::diagnostics::waveform::OneGainIntFunctionalWaveformRawStats* stats)
{
  auto* num_sum = stats->mutable_num_sum_entries()->mutable_data();
  auto* sum = stats->mutable_sum()->mutable_data();
  auto* sum_squared = stats->mutable_sum_squared()->mutable_data();
  const unsigned nchan = mask.size();
  for(unsigned ichan=0; ichan<nchan; ichan++)
    if(mask[ichan])
    {
      auto si = signal[ichan];
      num_sum[ichan]++;
      sum[ichan] += si;
      sum_squared[ichan] += SQR(si);
    }

  if(config_.calculate_covariance())
  {
    auto* num_sum_product =
      stats->mutable_num_sum_product_entries()->mutable_data();
    auto* sum_product = stats->mutable_sum_product()->mutable_data();
    unsigned idx = 0;
    for(unsigned ichan=0; ichan<nchan; ichan++)
      if(mask[ichan])
      {
        auto si = signal[ichan];
        for(unsigned jchan=ichan+1; jchan<nchan; jchan++, idx++)
          if(mask[jchan])
          {
            num_sum_product[idx]++;
            sum_product[idx] += si*signal[jchan];
          }

      }
    }
}

bool FunctionalWaveformStatsVisitor::merge_results()
{
  if(parent_)
  {
    if(results_.has_high_gain())
      merge_one_gain(&results_.high_gain(),
        parent_->results_.mutable_high_gain());
    if(results_.has_low_gain())
      merge_one_gain(&results_.low_gain(),
        parent_->results_.mutable_low_gain());
  }
  return true;
}

void FunctionalWaveformStatsVisitor::merge_one_gain(
  const ix::diagnostics::waveform::OneGainIntFunctionalWaveformRawStats* from,
  ix::diagnostics::waveform::OneGainIntFunctionalWaveformRawStats* to)
{
  assert(to->num_sum_entries_size() == from->num_sum_entries_size());
  assert(to->sum_size() == from->sum_size());
  assert(to->sum_squared_size() == from->sum_squared_size());
  assert(to->num_sum_product_entries_size() ==
    from->num_sum_product_entries_size());
  assert(to->sum_product_size() == from->sum_product_size());

  for(int i=0; i<from->num_sum_entries_size(); i++)
    to->set_num_sum_entries(i,
      to->num_sum_entries(i) + from->num_sum_entries(i));
  for(int i=0; i<from->sum_size(); i++)
    to->set_sum(i, to->sum(i) + from->sum(i));
  for(int i=0; i<from->sum_squared_size(); i++)
    to->set_sum_squared(i, to->sum_squared(i) + from->sum_squared(i));
  for(int i=0; i<from->num_sum_product_entries_size(); i++)
    to->set_num_sum_product_entries(i,
      to->num_sum_product_entries(i) + from->num_sum_product_entries(i));
  for(int i=0; i<from->sum_product_size(); i++)
    to->set_sum_product(i, to->sum_product(i) + from->sum_product(i));
}

Eigen::VectorXd FunctionalWaveformStatsVisitor::channel_mean(
  const ix::diagnostics::waveform::OneGainIntFunctionalWaveformRawStats* stat)
{
  const int N = stat->sum_size();
  assert(N == stat->num_sum_entries_size());
  Eigen::VectorXd m(N);
  for(int i=0; i<N; i++)
    m(i) = double(stat->sum(i)) / double(stat->num_sum_entries(i));
  return m;
}

Eigen::VectorXd FunctionalWaveformStatsVisitor::channel_var(
  const ix::diagnostics::waveform::OneGainIntFunctionalWaveformRawStats* stat)
{
  const int N = stat->sum_size();
  assert(N == stat->sum_squared_size());
  assert(N == stat->num_sum_entries_size());
  Eigen::VectorXd v(N);
  for(int i=0; i<N; i++)
    v(i) = double(stat->sum_squared(i)) / double(stat->num_sum_entries(i))
      - SQR(double(stat->sum(i)) / double(stat->num_sum_entries(i)));
  return v;
}

Eigen::MatrixXd FunctionalWaveformStatsVisitor::channel_cov(
  const ix::diagnostics::waveform::OneGainIntFunctionalWaveformRawStats* stat)
{
  const int N = stat->sum_size();
  assert(N == stat->sum_squared_size());
  assert(N == stat->num_sum_entries_size());
  assert(N*(N-1)/2 == stat->sum_product_size());
  assert(N*(N-1)/2 == stat->num_sum_product_entries_size());

  Eigen::MatrixXd c(N,N);
  for(int i=0; i<N; i++)
    c(i,i) = double(stat->sum_squared(i)) / double(stat->num_sum_entries(i))
      - SQR(double(stat->sum(i)) / double(stat->num_sum_entries(i)));
  for(int i=0; i<N; i++)
    for(int j=i+1; j<N; j++)
    {
      int idx = N*(N-1)/2-(N-i)*(N-i-1)/2 + j - i - 1;
      double cij =
        double(stat->sum_product(idx)) /
          double(stat->num_sum_product_entries(idx))
        - double(stat->sum(i))*double(stat->sum(j)) /
          (double(stat->num_sum_entries(i))*double(stat->num_sum_entries(j)));
      c(i,j) = cij;
      c(j,i) = cij;
    }
  return c;
}
