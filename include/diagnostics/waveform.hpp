  /*

   calin/diagnostics/waveform.hpp -- Stephen Fegan -- 2016-02-10

   Waveform diagnostics visitor

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <fftw3.h>

#include <iact_data/event_visitor.hpp>
#include <iact_data/functional_event_visitor.hpp>
#include <diagnostics/waveform.pb.h>
#include <math/histogram.hpp>
// #include <math/simd.hpp>
// #include <math/fft_simd.hpp>

namespace calin { namespace diagnostics { namespace waveform {

class WaveformSumParallelEventVisitor:
public iact_data::event_visitor::ParallelEventVisitor
{
public:
  WaveformSumParallelEventVisitor();

  virtual ~WaveformSumParallelEventVisitor();

  WaveformSumParallelEventVisitor* new_sub_visitor(
    std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
        calin::iact_data::event_visitor::ParallelEventVisitor*>
      antecedent_visitors) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
    calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager) override;
  bool leave_telescope_run() override;

  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;

  bool merge_results() override;

#ifndef SWIG
  calin::ix::diagnostics::waveform::WaveformMean* mean_waveforms(
    calin::ix::diagnostics::waveform::WaveformMean* wf = nullptr) const;
#else
  calin::ix::diagnostics::waveform::WaveformMean* mean_waveforms() const;
  void mean_waveforms(calin::ix::diagnostics::waveform::WaveformMean* wf) const;
#endif

  void high_gain_mean_wf(Eigen::MatrixXd& mean_waveform_out);
  void low_gain_mean_wf(Eigen::MatrixXd& mean_waveform_out);
  void high_gain_event_count(Eigen::VectorXi& count_out);
  void low_gain_event_count(Eigen::VectorXi& count_out);

protected:
#ifndef SWIG
  WaveformSumParallelEventVisitor* parent_ = nullptr;

  unsigned partial_max_num_entries_ = 32768;

  bool has_dual_gain_ = false;
  unsigned nchan_ = 0;
  unsigned nsamp_ = 0;

  uint32_t num_entries_i32_ = 0;

  uint64_t*__restrict__ high_gain_count_i64_ = nullptr;
  uint64_t*__restrict__ low_gain_count_i64_ = nullptr;
  int64_t*__restrict__ high_gain_wf_sum_i64_ = nullptr;
  int64_t*__restrict__ low_gain_wf_sum_i64_ = nullptr;

  uint32_t*__restrict__ high_gain_count_i32_ = nullptr;
  uint32_t*__restrict__ low_gain_count_i32_ = nullptr;
  int32_t*__restrict__ high_gain_wf_sum_i32_ = nullptr;
  int32_t*__restrict__ low_gain_wf_sum_i32_ = nullptr;
#endif // ndef SWIG
};

class WaveformCodeHistParallelEventVisitor:
public iact_data::event_visitor::ParallelEventVisitor
{
public:
  WaveformCodeHistParallelEventVisitor(bool max_sample_only = false, uint16_t max_code = 4095);

  virtual ~WaveformCodeHistParallelEventVisitor();

  WaveformCodeHistParallelEventVisitor* new_sub_visitor(
    std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
        calin::iact_data::event_visitor::ParallelEventVisitor*>
      antecedent_visitors) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
    calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager) override;
  bool leave_telescope_run() override;

  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;

  bool merge_results() override;

  calin::ix::diagnostics::waveform::WaveformCodeHist* waveform_code_hist() const;
  calin::ix::diagnostics::waveform::CompactWaveformCodeHist* compact_waveform_code_hist() const;

protected:
#ifndef SWIG
  void analyze_wf_image(const calin::ix::iact_data::telescope_event::Waveforms& wf);
  void transfer_u32_to_u64();

  WaveformCodeHistParallelEventVisitor* parent_ = nullptr;

  bool max_sample_only_ = false;
  uint16_t max_code_ = 4095;

  calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config_ = nullptr;
  bool has_dual_gain_ = false;
  unsigned nchan_ = 0;
  unsigned nsamp_ = 0;

  unsigned max_nsamp_in_u32_ = 0;
  unsigned nsamp_in_u32_ = 0;

  uint32_t*__restrict__ high_gain_code_hist_u32_ = nullptr;
  uint32_t*__restrict__ low_gain_code_hist_u32_ = nullptr;
  uint64_t*__restrict__ high_gain_code_hist_u64_ = nullptr;
  uint64_t*__restrict__ low_gain_code_hist_u64_ = nullptr;
#endif // ndef SWIG
};



class WaveformStatsParallelVisitor:
  public iact_data::event_visitor::ParallelEventVisitor
{
public:
  WaveformStatsParallelVisitor(bool calculate_psd = true,
    bool calculate_covariance = true);

  virtual ~WaveformStatsParallelVisitor();

  WaveformStatsParallelVisitor* new_sub_visitor(
    std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
        calin::iact_data::event_visitor::ParallelEventVisitor*>
      antecedent_visitors) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
    calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager) override;
  bool leave_telescope_run() override;

  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;

  bool merge_results() override;

  calin::ix::diagnostics::waveform::CameraWaveformRawStats results()
  {
    return results_;
  }

  calin::ix::diagnostics::waveform::CameraWaveformRawPSD psd_results()
  {
    return psd_results_;
  }

  static Eigen::VectorXd waveform_mean(
    const ix::diagnostics::waveform::WaveformRawStats* stat);
  static Eigen::VectorXd waveform_var(
    const ix::diagnostics::waveform::WaveformRawStats* stat);
  static Eigen::MatrixXd waveform_cov(
    const ix::diagnostics::waveform::WaveformRawStats* stat);
  static Eigen::MatrixXd waveform_cov_frac(
    const ix::diagnostics::waveform::WaveformRawStats* stat);

#ifndef SWIG
protected:
  void process_one_waveform(const uint16_t*__restrict__ wf,
    ix::diagnostics::waveform::PartialWaveformRawStats* p_stat,
    ix::diagnostics::waveform::WaveformRawStats* r_stat,
    ix::diagnostics::waveform::WaveformRawPSD* psd = nullptr);

  void merge_partial(
    ix::diagnostics::waveform::PartialWaveformRawStats* p_stat,
    ix::diagnostics::waveform::WaveformRawStats* r_stat);

  WaveformStatsParallelVisitor* parent_ = nullptr;
  calin::ix::diagnostics::waveform::CameraWaveformRawStats results_;
  calin::ix::diagnostics::waveform::PartialCameraWaveformRawStats partial_;
  unsigned partial_max_num_entries_ = 256;
  const ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration*
    run_config_ = nullptr;

  calin::ix::diagnostics::waveform::CameraWaveformRawPSD psd_results_;
  float*__restrict__ waveform_t_ = nullptr;
  float*__restrict__ waveform_f_ = nullptr;
  fftwf_plan fftw_plan_fwd_;
  fftwf_plan fftw_plan_bwd_;

  bool calculate_psd_ = false;
  bool calculate_covariance_ = false;
#endif
};

class WaveformPSDVisitor:
  public iact_data::event_visitor::TelescopeEventVisitor
{
public:
  WaveformPSDVisitor();

  virtual ~WaveformPSDVisitor();

  bool demand_waveforms() override;
  bool is_parallelizable() override;
  WaveformPSDVisitor* new_sub_visitor(
    const std::map<calin::iact_data::event_visitor::TelescopeEventVisitor*,
        calin::iact_data::event_visitor::TelescopeEventVisitor*>&
      antecedent_visitors) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config) override;
  bool leave_telescope_run() override;

  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;

  bool visit_waveform(unsigned ichan,
    calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
    calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain) override;

  bool merge_results() override;

  calin::ix::diagnostics::waveform::CameraWaveformRawPSD results()
  {
    return results_;
  }

  static Eigen::VectorXd psd_mean(
    const ix::diagnostics::waveform::WaveformRawPSD* stat);
  static Eigen::VectorXd psd_var(
    const ix::diagnostics::waveform::WaveformRawPSD* stat);

  static Eigen::VectorXd corr_mean(
    const ix::diagnostics::waveform::WaveformRawPSD* psd_stat,
    const ix::diagnostics::waveform::WaveformRawStats* trace_stat);
  static Eigen::VectorXd corr_mean_centered(
    const ix::diagnostics::waveform::WaveformRawPSD* stat,
    const ix::diagnostics::waveform::WaveformRawStats* trace_stat,
    Eigen::VectorXd& h);
  static Eigen::VectorXd corr_var(
    const ix::diagnostics::waveform::WaveformRawPSD* stat,
    const ix::diagnostics::waveform::WaveformRawStats* trace_stat);

protected:
  void process_one_waveform(
    const calin::ix::iact_data::telescope_event::ChannelWaveform* wf,
    ix::diagnostics::waveform::WaveformRawPSD* psd);

  WaveformPSDVisitor* parent_ = nullptr;
  calin::ix::diagnostics::waveform::CameraWaveformRawPSD results_;
  const ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration*
    run_config_ = nullptr;
  double* fftw_data_ = nullptr;
  fftw_plan fftw_plan_fwd_ = nullptr;
  fftw_plan fftw_plan_bwd_ = nullptr;
};

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

// OBSOLETE version that uses TelescopeEventVisitor interface

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

class WaveformStatsVisitor:
  public iact_data::event_visitor::TelescopeEventVisitor
{
public:
  WaveformStatsVisitor(bool calculate_covariance = true);

  virtual ~WaveformStatsVisitor();

  bool demand_waveforms() override;
  bool is_parallelizable() override;
  WaveformStatsVisitor* new_sub_visitor(
    const std::map<calin::iact_data::event_visitor::TelescopeEventVisitor*,
        calin::iact_data::event_visitor::TelescopeEventVisitor*>&
      antecedent_visitors) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config) override;
  bool leave_telescope_run() override;

  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;

  bool visit_waveform(unsigned ichan,
    calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
    calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain) override;

  bool merge_results() override;

  calin::ix::diagnostics::waveform::CameraWaveformRawStats results()
  {
    return results_;
  }

  static Eigen::VectorXd waveform_mean(
    const ix::diagnostics::waveform::WaveformRawStats* stat);
  static Eigen::VectorXd waveform_var(
    const ix::diagnostics::waveform::WaveformRawStats* stat);
  static Eigen::MatrixXd waveform_cov(
    const ix::diagnostics::waveform::WaveformRawStats* stat);
  static Eigen::MatrixXd waveform_cov_frac(
    const ix::diagnostics::waveform::WaveformRawStats* stat);

protected:
  void process_one_waveform(
    const calin::ix::iact_data::telescope_event::ChannelWaveform* wf,
    ix::diagnostics::waveform::PartialWaveformRawStats* p_stat,
    ix::diagnostics::waveform::WaveformRawStats* r_stat);

  void merge_partial(
    ix::diagnostics::waveform::PartialWaveformRawStats* p_stat,
    ix::diagnostics::waveform::WaveformRawStats* r_stat);

  WaveformStatsVisitor* parent_ = nullptr;
  calin::ix::diagnostics::waveform::CameraWaveformRawStats results_;
  calin::ix::diagnostics::waveform::PartialCameraWaveformRawStats partial_;
  unsigned partial_max_num_entries_ = 256;
  const ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration*
    run_config_ = nullptr;
  bool calculate_covariance_ = false;
};


} } } // namespace calin::diagnostics::waveform
