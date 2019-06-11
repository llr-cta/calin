  /*

   calin/diagnostics/waveform.hpp -- Stephen Fegan -- 2016-02-10

   Waveform diagnostics visitor

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

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
#include <math/simd.hpp>
#include <math/fft_simd.hpp>

namespace calin { namespace diagnostics { namespace waveform {

class WaveformStatsParallelVisitor:
  public iact_data::event_visitor::ParallelEventVisitor
{
public:
  WaveformStatsParallelVisitor(bool calculate_psd = true,
    bool calculate_covariance = true);

  virtual ~WaveformStatsParallelVisitor();

  WaveformStatsParallelVisitor* new_sub_visitor(
    const std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
        calin::iact_data::event_visitor::ParallelEventVisitor*>&
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

class AVX2_Unroll8_WaveformStatsParallelVisitor:
  public iact_data::event_visitor::ParallelEventVisitor
{
public:
  AVX2_Unroll8_WaveformStatsParallelVisitor(bool high_gain = true,
    bool calculate_psd = true, bool calculate_covariance = true);

  virtual ~AVX2_Unroll8_WaveformStatsParallelVisitor();

  AVX2_Unroll8_WaveformStatsParallelVisitor* new_sub_visitor(
    const std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
        calin::iact_data::event_visitor::ParallelEventVisitor*>&
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
#if defined(__AVX2__)
    return results_;
#else // defined(__AVX2__)
    throw std::runtime_error("AVX2_Unroll8_WaveformStatsParallelVisitor: AVX2 not supported at compile time");
#endif // defined(__AVX2__)
  }

protected:
#if defined(__AVX2__) and not defined(SWIG)
  void do_8_dft(__m256d*__restrict__& psd_sum, __m256d*__restrict__& psd_sumsq,
    __m256d*__restrict__& ac_sum, __m256d*__restrict__& ac_sumsq);
  void process_8_events();
  void merge_partials();

  AVX2_Unroll8_WaveformStatsParallelVisitor* parent_ = nullptr;
  calin::ix::diagnostics::waveform::CameraWaveformRawStats results_;

  calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager_ = nullptr;

  unsigned nkept_events_ = 0;
  calin::ix::iact_data::telescope_event::TelescopeEvent* kept_events_[8] =
    { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr };

  __m256i* samples_[8] =
    { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr };

  __m256* dft_xt_[8] =
    { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr };

  __m256* dft_xf_[8] =
    { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr };

  calin::math::fft_simd::FixedSizeRealToHalfComplexDFT<__m256>* dft_ = nullptr;

  unsigned nchan_ = 0;
  unsigned nsamp_ = 0;

  unsigned partial_num_entries_ = 0;
  unsigned* partial_chan_nevent_ = nullptr;
  __m256i* partial_chan_sum_ = nullptr;
  __m256i* partial_chan_sum_squared_ = nullptr;
  __m256i* partial_chan_sum_cov_ = nullptr;
#endif // defined(__AVX2__) and not defined(SWIG)

  unsigned partial_max_num_entries_ = 256;
  const ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration*
    run_config_ = nullptr;

  bool high_gain_ = false;
  bool calculate_psd_ = false;
  bool calculate_covariance_ = false;
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
