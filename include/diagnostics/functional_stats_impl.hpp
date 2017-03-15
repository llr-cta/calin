/*

   calin/diagnostics/functional_stats_impl.hpp -- Stephen Fegan -- 2016-05-10

   Functional diagnostics visitor - integrated stats with histo and
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

#include <math/special.hpp>
#include <math/covariance_calc.hpp>

namespace calin { namespace diagnostics { namespace functional {

template<typename DualGainFunctionalVisitor, typename Results>
FunctionalStatsVisitor<DualGainFunctionalVisitor, Results>::
FunctionalStatsVisitor(
    DualGainFunctionalVisitor* value_supplier,
    const ix::diagnostics::functional::FunctionalStatsVisitorConfig& config):
  TelescopeEventVisitor(), value_supplier_(value_supplier), config_(config)
{
  // nothing to see here
}

template<typename DualGainFunctionalVisitor, typename Results>
FunctionalStatsVisitor<DualGainFunctionalVisitor, Results>::
~FunctionalStatsVisitor()
{
  // nothing to see here
}

template<typename DualGainFunctionalVisitor, typename Results>
bool FunctionalStatsVisitor<DualGainFunctionalVisitor, Results>::
demand_waveforms()
{
  return true;
}

template<typename DualGainFunctionalVisitor, typename Results>
bool FunctionalStatsVisitor<DualGainFunctionalVisitor, Results>::
is_parallelizable()
{
  return true;
}

template<typename DualGainFunctionalVisitor, typename Results>
FunctionalStatsVisitor<DualGainFunctionalVisitor, Results>*
FunctionalStatsVisitor<DualGainFunctionalVisitor, Results>::new_sub_visitor(
  const std::map<TelescopeEventVisitor*,TelescopeEventVisitor*>&
    antecedent_visitors)
{
  auto i_sub_value_supplier = antecedent_visitors.find(value_supplier_);
  if(i_sub_value_supplier == antecedent_visitors.end())
    throw std::logic_error("No thread specific value suppler!");
  auto* sub_value_supplier =
    dynamic_cast<DualGainFunctionalVisitor*>(
      i_sub_value_supplier->second);
  auto* sub_visitor =
    new FunctionalStatsVisitor(sub_value_supplier, config_);
  sub_visitor->parent_ = this;
  return sub_visitor;
}

template<typename DualGainFunctionalVisitor, typename Results>
bool FunctionalStatsVisitor<DualGainFunctionalVisitor, Results>::
visit_telescope_run(
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
  high_gain_hist_.clear();
  for(int ichan=0; ichan<nchan; ichan++)
  {
    calin::ix::math::histogram::Histogram1DConfig hist_config(
      config_.hist_config());
    if(hist_config.name().empty())
      hist_config.set_name(
        std::string("Channel ")+std::to_string(ichan)+" (high gain)");
    else hist_config.set_name(
      hist_config.name() + " (channel "+std::to_string(ichan)+", high gain)");
    if(hist_config.weight_units().empty())
      hist_config.set_weight_units("events");
    high_gain_hist_.emplace_back(hist_config);
  }

  auto* hg_results = results_.mutable_high_gain();
  hg_results->mutable_num_sum_entries()->Resize(nchan,0);
  hg_results->mutable_num_sum_product_entries()->Resize(nchan*(nchan-1)/2,0);
  hg_results->mutable_sum()->Resize(nchan,0);
  hg_results->mutable_sum_squared()->Resize(nchan,0);
  hg_results->mutable_sum_product()->Resize(nchan*(nchan-1)/2,0);
  for(int ichan=0; ichan<nchan; ichan++)hg_results->add_value_hist();

  low_gain_mask_.resize(nchan);
  low_gain_signal_.resize(nchan);
  low_gain_hist_.clear();
  for(int ichan=0; ichan<nchan; ichan++)
  {
    calin::ix::math::histogram::Histogram1DConfig hist_config(
      config_.hist_config());
    if(hist_config.name().empty())
      hist_config.set_name(
        std::string("Channel ")+std::to_string(ichan)+" (high gain)");
    else hist_config.set_name(
      hist_config.name() + " (channel "+std::to_string(ichan)+", high gain)");
    if(hist_config.weight_units().empty())
      hist_config.set_weight_units("events");
    low_gain_hist_.emplace_back(hist_config);
  }

  auto* lg_results = results_.mutable_low_gain();
  lg_results->mutable_num_sum_entries()->Resize(nchan,0);
  lg_results->mutable_num_sum_product_entries()->Resize(nchan*(nchan-1)/2,0);
  lg_results->mutable_sum()->Resize(nchan,0);
  lg_results->mutable_sum_squared()->Resize(nchan,0);
  lg_results->mutable_sum_product()->Resize(nchan*(nchan-1)/2,0);
  for(int ichan=0; ichan<nchan; ichan++)lg_results->add_value_hist();

  return true;
}

template<typename DualGainFunctionalVisitor, typename Results>
bool FunctionalStatsVisitor<DualGainFunctionalVisitor, Results>::
leave_telescope_run()
{
  for(unsigned ichan=0; ichan<high_gain_hist_.size(); ichan++) {
    auto* hist = high_gain_hist_[ichan].dump_as_proto();
    calin::math::histogram::merge_histogram1d_data(
      results_.mutable_high_gain()->mutable_value_hist(ichan), *hist);
    delete hist;
  }
  for(unsigned ichan=0; ichan<low_gain_hist_.size(); ichan++) {
    auto* hist = low_gain_hist_[ichan].dump_as_proto();
    calin::math::histogram::merge_histogram1d_data(
      results_.mutable_low_gain()->mutable_value_hist(ichan), *hist);
    delete hist;
  }
  run_config_ = nullptr;
  return true;
}

template<typename DualGainFunctionalVisitor, typename Results>
bool FunctionalStatsVisitor<DualGainFunctionalVisitor, Results>::
visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  for(auto& imask : high_gain_mask_)imask = 0;
  for(auto& imask : low_gain_mask_)imask = 0;
  return true;
}

template<typename DualGainFunctionalVisitor, typename Results>
bool FunctionalStatsVisitor<DualGainFunctionalVisitor, Results>::
leave_telescope_event()
{
  if(results_.has_high_gain())
    process_one_gain(high_gain_mask_, high_gain_signal_, high_gain_hist_,
      results_.mutable_high_gain());
  if(results_.has_low_gain())
    process_one_gain(low_gain_mask_, low_gain_signal_, low_gain_hist_,
      results_.mutable_low_gain());
  return true;
}

template<typename DualGainFunctionalVisitor, typename Results>
bool FunctionalStatsVisitor<DualGainFunctionalVisitor, Results>::
visit_waveform(unsigned ichan,
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

template<typename DualGainFunctionalVisitor, typename Results>
void FunctionalStatsVisitor<DualGainFunctionalVisitor, Results>::
visit_one_waveform(
  const calin::ix::iact_data::telescope_event::ChannelWaveform* wf,
  unsigned index, std::vector<int>& mask,
  std::vector<functional_value_type>& signal)
{
  mask[index] = 1;
}

template<typename DualGainFunctionalVisitor, typename Results>
template<typename OneGainRawStats>
void FunctionalStatsVisitor<DualGainFunctionalVisitor, Results>::
process_one_gain(const std::vector<int>& mask,
  const std::vector<functional_value_type>& signal,
  std::vector<calin::math::histogram::SimpleHist>& hist,
  OneGainRawStats* stats)
{
  using calin::math::special::SQR;
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
      hist[ichan].insert(double(si));
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
      else
      {
        idx += nchan - ichan - 1;
      }
    }
}

template<typename DualGainFunctionalVisitor, typename Results>
bool FunctionalStatsVisitor<DualGainFunctionalVisitor, Results>::
merge_results()
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

template<typename DualGainFunctionalVisitor, typename Results>
template<typename OneGainRawStats>
void FunctionalStatsVisitor<DualGainFunctionalVisitor, Results>::
merge_one_gain(const OneGainRawStats* from, OneGainRawStats* to)
{
  assert(to->num_sum_entries_size() == from->num_sum_entries_size());
  assert(to->sum_size() == from->sum_size());
  assert(to->sum_squared_size() == from->sum_squared_size());
  assert(to->num_sum_product_entries_size() ==
    from->num_sum_product_entries_size());
  assert(to->sum_product_size() == from->sum_product_size());
  assert(to->value_hist_size() == from->value_hist_size());

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
  for(int i=0; i<from->value_hist_size(); i++)
    to->mutable_value_hist(i)->IntegrateFrom(from->value_hist(i));
}

template<typename OneGainRawStats>
Eigen::VectorXd channel_mean(const OneGainRawStats* stat)
{
  const int N = stat->sum_size();
  assert(N == stat->num_sum_entries_size());
  Eigen::VectorXd m(N);
  for(int i=0; i<N; i++)
    m(i) = double(stat->sum(i)) / double(stat->num_sum_entries(i));
  return m;
}

template<typename OneGainRawStats>
Eigen::VectorXd channel_var(const OneGainRawStats* stat)
{
  using calin::math::special::SQR;
  using calin::math::covariance_calc::cov_gen;
  const int N = stat->sum_size();
  assert(N == stat->sum_squared_size());
  assert(N == stat->num_sum_entries_size());
  Eigen::VectorXd v(N);
  for(int i=0; i<N; i++)
    v(i) = cov_gen(stat->sum_squared(i), stat->num_sum_entries(i),
      stat->sum(i), stat->num_sum_entries(i),
      stat->sum(i), stat->num_sum_entries(i));
  return v;
}

template<typename OneGainRawStats>
Eigen::MatrixXd channel_cov(const OneGainRawStats* stat)
{
  using calin::math::special::SQR;
  using calin::math::covariance_calc::cov_gen;
  const int N = stat->sum_size();
  assert(N == stat->sum_squared_size());
  assert(N == stat->num_sum_entries_size());
  assert(N*(N-1)/2 == stat->sum_product_size());
  assert(N*(N-1)/2 == stat->num_sum_product_entries_size());

  Eigen::MatrixXd c(N,N);
  for(int i=0; i<N; i++)
    c(i,i) = cov_gen(stat->sum_squared(i), stat->num_sum_entries(i),
      stat->sum(i), stat->num_sum_entries(i),
      stat->sum(i), stat->num_sum_entries(i));
  int idx = 0;
  for(int i=0; i<N; i++)
    for(int j=i+1; j<N; j++, idx++)
    {
      double cij =
        cov_gen(stat->sum_product(idx), stat->num_sum_product_entries(idx),
          stat->sum(i), stat->num_sum_entries(i),
          stat->sum(j), stat->num_sum_entries(j));
      c(i,j) = cij;
      c(j,i) = cij;
    }
  return c;
}

template<typename OneGainRawStats>
Eigen::MatrixXd channel_cov_frac(const OneGainRawStats* stat)
{
  Eigen::MatrixXd c = channel_cov(stat);
  const int N = stat->sum_size();
  for(int i=0; i<N; i++) {
    double scale = 1.0/std::sqrt(c(i,i));
    for(int j=0; j<N; j++){
      c(i,j) *= scale;
      c(j,i) *= scale;
    }
  }
  return c;
}

} } } // namespace calin::diagnostics::functional_diagnostics
