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

#include <type_traits>

#include <Eigen/Dense>
#include <Eigen/LU>

#include <math/special.hpp>
#include <math/covariance_calc.hpp>

namespace calin { namespace diagnostics { namespace functional {

template<typename DualGainFunctionalVisitor, typename Results>
FunctionalStatsVisitor<DualGainFunctionalVisitor, Results>::
FunctionalStatsVisitor(
    DualGainFunctionalVisitor* value_supplier,
    const ix::diagnostics::functional::FunctionalStatsVisitorConfig& config):
  TelescopeEventVisitor(), value_supplier_(value_supplier), config_(config),
  high_gain_mean_hist_(config.hist_config()),
  low_gain_mean_hist_(config.hist_config())
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

  calin::ix::math::histogram::Histogram1DConfig hg_mean_hist_config(
    config_.hist_config());
  if(hg_mean_hist_config.name().empty())
    hg_mean_hist_config.set_name("Mean (high gain)");
  else hg_mean_hist_config.set_name(
    hg_mean_hist_config.name() + " (mean high gain)");
  if(hg_mean_hist_config.weight_units().empty())
    hg_mean_hist_config.set_weight_units("events");
  if(config_.mean_hist_dxval_multiplier() == 0.0)
    hg_mean_hist_config.set_dxval(hg_mean_hist_config.dxval()
      / std::sqrt(nchan));
  else
    hg_mean_hist_config.set_dxval(hg_mean_hist_config.dxval()
      * config_.mean_hist_dxval_multiplier());
  high_gain_mean_hist_ = calin::math::histogram::SimpleHist(hg_mean_hist_config);

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

  calin::ix::math::histogram::Histogram1DConfig lg_mean_hist_config(
    config_.hist_config());
  if(lg_mean_hist_config.name().empty())
    lg_mean_hist_config.set_name("Mean (low gain)");
  else lg_mean_hist_config.set_name(
    lg_mean_hist_config.name() + " (mean low gain)");
  if(lg_mean_hist_config.weight_units().empty())
    lg_mean_hist_config.set_weight_units("events");
  if(config_.mean_hist_dxval_multiplier() == 0.0)
    lg_mean_hist_config.set_dxval(lg_mean_hist_config.dxval()
      / std::sqrt(nchan));
  else
    lg_mean_hist_config.set_dxval(lg_mean_hist_config.dxval()
      * config_.mean_hist_dxval_multiplier());
  low_gain_mean_hist_ = calin::math::histogram::SimpleHist(lg_mean_hist_config);

  low_gain_mask_.resize(nchan);
  low_gain_signal_.resize(nchan);
  low_gain_hist_.clear();
  for(int ichan=0; ichan<nchan; ichan++)
  {
    calin::ix::math::histogram::Histogram1DConfig hist_config(
      config_.hist_config());
    if(hist_config.name().empty())
      hist_config.set_name(
        std::string("Channel ")+std::to_string(ichan)+" (low gain)");
    else hist_config.set_name(
      hist_config.name() + " (channel "+std::to_string(ichan)+", low gain)");
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

  results_.mutable_num_sum_high_low_gain_product_entries()->Resize(nchan,0);
  results_.mutable_sum_high_low_gain_product()->Resize(nchan,0);

  return true;
}

template<typename DualGainFunctionalVisitor, typename Results>
bool FunctionalStatsVisitor<DualGainFunctionalVisitor, Results>::
leave_telescope_run()
{
  for(unsigned ichan=0; ichan<high_gain_hist_.size(); ichan++) {
    auto* hist = high_gain_hist_[ichan].dump_as_proto();
    results_.mutable_high_gain()->mutable_value_hist(ichan)->IntegrateFrom(*hist);
    delete hist;
  }
  {
    auto* hist = high_gain_mean_hist_.dump_as_proto();
    results_.mutable_high_gain()->mutable_mean_hist()->IntegrateFrom(*hist);
    delete hist;
  }
  for(unsigned ichan=0; ichan<low_gain_hist_.size(); ichan++) {
    auto* hist = low_gain_hist_[ichan].dump_as_proto();
    results_.mutable_low_gain()->mutable_value_hist(ichan)->IntegrateFrom(*hist);
    delete hist;
  }
  {
    auto* hist = low_gain_mean_hist_.dump_as_proto();
    results_.mutable_low_gain()->mutable_mean_hist()->IntegrateFrom(*hist);
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
      high_gain_mean_hist_, results_.mutable_high_gain());
  if(results_.has_low_gain())
    process_one_gain(low_gain_mask_, low_gain_signal_, low_gain_hist_,
      low_gain_mean_hist_, results_.mutable_low_gain());
  if(results_.has_high_gain() and results_.has_low_gain())
  {
    auto* num_sum = results_.mutable_num_sum_high_low_gain_product_entries()->mutable_data();
    auto* sum = results_.mutable_sum_high_low_gain_product()->mutable_data();
    const unsigned nchan = std::min(high_gain_mask_.size(), low_gain_mask_.size());
    for(unsigned ichan=0; ichan<nchan; ichan++)
      if(high_gain_mask_[ichan] and low_gain_mask_[ichan])
      {
        auto si = high_gain_signal_[ichan] * low_gain_signal_[ichan];
        num_sum[ichan]++;
        sum[ichan] += si;
      }
  }

  return true;
}

template<typename DualGainFunctionalVisitor, typename Results>
bool FunctionalStatsVisitor<DualGainFunctionalVisitor, Results>::
visit_waveform(unsigned ichan,
  calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
  calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain)
{
  const int index = ichan; // run_config_->configured_channel_index(ichan);
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
template<typename OneGainRawStats>
void FunctionalStatsVisitor<DualGainFunctionalVisitor, Results>::
process_one_gain(const std::vector<int>& mask,
  const std::vector<functional_value_type>& signal,
  std::vector<calin::math::histogram::SimpleHist>& hist,
  calin::math::histogram::SimpleHist& mean_hist,
  OneGainRawStats* stats)
{
  using calin::math::special::SQR;
  auto* num_sum = stats->mutable_num_sum_entries()->mutable_data();
  auto* sum = stats->mutable_sum()->mutable_data();
  auto* sum_squared = stats->mutable_sum_squared()->mutable_data();
  const unsigned nchan = mask.size();
  unsigned nmean = 0;
  typename std::remove_reference<decltype(*sum)>::type mean_sum = 0;
  for(unsigned ichan=0; ichan<nchan; ichan++)
    if(mask[ichan])
    {
      auto si = signal[ichan];
      num_sum[ichan]++;
      sum[ichan] += si;
      sum_squared[ichan] += SQR(si);
      hist[ichan].insert(double(si));
      nmean++;
      mean_sum += si;
    }

  if(nmean) {
    double mean = double(mean_sum)/double(nmean);
    stats->set_num_sum_mean_entries(stats->num_sum_mean_entries() + 1);
    stats->set_sum_mean(stats->sum_mean() + mean);
    stats->set_sum_mean_squared(stats->sum_mean_squared() + SQR(mean));
    mean_hist.insert(mean);
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
  if(parent_)parent_->results_.IntegrateFrom(results_);
  return true;
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

template<typename DualGainRawStats>
Eigen::VectorXd channel_high_to_low_gain_cov(const DualGainRawStats* stat)
{
  using calin::math::special::SQR;
  using calin::math::covariance_calc::cov_gen;
  const int N = stat->num_sum_high_low_gain_product_entries_size();
  assert(N == stat->sum_high_low_gain_product_size());
  assert(stat->has_high_gain());
  assert(N == stat->high_gain().num_sum_entries_size());
  assert(N == stat->high_gain().sum_size());
  assert(stat->has_low_gain());
  assert(N == stat->low_gain().num_sum_entries_size());
  assert(N == stat->low_gain().sum_size());
  Eigen::VectorXd v(N);
  for(int i=0; i<N; i++)
    v(i) = cov_gen(stat->sum_high_low_gain_product(i),
      stat->num_sum_high_low_gain_product_entries(i),
      stat->high_gain().sum(i), stat->high_gain().num_sum_entries(i),
      stat->low_gain().sum(i), stat->low_gain().num_sum_entries(i));
  return v;
}

template<typename DualGainRawStats>
Eigen::VectorXd channel_high_to_low_gain_cov_frac(const DualGainRawStats* stat)
{
  Eigen::VectorXd v = channel_high_to_low_gain_cov(stat);
  v = v.array() / channel_var(&stat->high_gain()).array().sqrt();
  v = v.array() / channel_var(&stat->low_gain()).array().sqrt();
  return v;
}

template<typename OneGainRawStats>
double mean_of_mean_over_channels(const OneGainRawStats* stat)
{
  return stat->sum_mean() / double(stat->num_sum_mean_entries());
}

template<typename OneGainRawStats>
double var_of_mean_over_channels(const OneGainRawStats* stat)
{
  using calin::math::covariance_calc::cov_gen;
  return cov_gen(stat->sum_mean_squared(), stat->num_sum_mean_entries(),
    stat->sum_mean(), stat->num_sum_mean_entries(),
    stat->sum_mean(), stat->num_sum_mean_entries());
}

template<typename OneGainRawStats>
Eigen::VectorXd decompose_channel_independent_and_common_var(
  const OneGainRawStats* stat, double& common_variance_out)
{
  using calin::math::special::SQR;
  const int N = stat->sum_size();
  Eigen::MatrixXd M = Eigen::MatrixXd::Identity(N+1,N+1);
  M.block(0,N,N,1) = Eigen::MatrixXd::Ones(N,1);
  M.block(N,0,1,N) = Eigen::MatrixXd::Ones(1,N);
  M(N,N) = SQR(double(N));
  Eigen::VectorXd V(N+1);
  V.head(N) = channel_var(stat);
  V(N) = var_of_mean_over_channels(stat) * SQR(double(N));
  V = M.inverse() * V;
  common_variance_out = V[N];
  return V.head(N);
}

} } } // namespace calin::diagnostics::functional_diagnostics
