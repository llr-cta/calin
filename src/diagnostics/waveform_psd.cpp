/*

   calin/diagnostics/waveform.cpp -- Stephen Fegan -- 2021-09-24

   Waveform PSD diagnostics visitor

   Copyright 2021, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <math/special.hpp>
#include <util/log.hpp>
#include <io/json.hpp>
#include <util/memory.hpp>
#include <diagnostics/waveform.hpp>
#include <math/covariance_calc.hpp>
#include <math/fftw_util.hpp>

using namespace calin::util::log;
using namespace calin::diagnostics::waveform;
using namespace calin::math::fftw_util;
using calin::math::special::SQR;
using calin::math::covariance_calc::cov_i64_gen;
using calin::util::memory::safe_aligned_recalloc_and_fill;

WaveformPSDParallelVisitor::WaveformPSDParallelVisitor():
  ParallelEventVisitor()
{
  // nothing to see here
}

WaveformPSDParallelVisitor::~WaveformPSDParallelVisitor()
{
  // fftwf_destroy_plan(fftw_plan_bwd_);
  fftwf_destroy_plan(fftw_plan_fwd_);
  fftwf_free(waveform_t_);
  fftwf_free(waveform_f_);
}

WaveformPSDParallelVisitor* WaveformPSDParallelVisitor::new_sub_visitor(
  std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
      calin::iact_data::event_visitor::ParallelEventVisitor*>
    antecedent_visitors)
{
  auto* sub_visitor = new WaveformPSDParallelVisitor();
  sub_visitor->parent_ = this;
  return sub_visitor;
}

bool WaveformPSDParallelVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
  calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager,
  calin::ix::provenance::chronicle::ProcessingRecord* processing_record)
{
  if(processing_record) {
    processing_record->set_type("WaveformPSDParallelVisitor");
    processing_record->set_description("Waveform PSD");
  }

  run_config_ = run_config;
  has_dual_gain_ = (run_config->camera_layout().adc_gains() !=
    calin::ix::iact_data::instrument_layout::CameraLayout::SINGLE_GAIN);

  results_.Clear();

  unsigned N = run_config->num_samples();
  unsigned nfreq = 1 + N/2; // 1->1 2->2 3->2 4->3 5->3 etc
  free(waveform_t_);
  free(waveform_f_);
  waveform_t_ = fftwf_alloc_real(N);
  waveform_f_ = fftwf_alloc_real(N);
  assert(waveform_t_);
  assert(waveform_f_);
  fftw_plan_fwd_ =
    fftwf_plan_r2r_1d(N, waveform_t_, waveform_f_, FFTW_R2HC, FFTW_MEASURE|FFTW_DESTROY_INPUT);
  assert(fftw_plan_fwd_);
  // fftw_plan_bwd_ =
  //   fftwf_plan_r2r_1d(N, waveform_f_, waveform_t_, FFTW_HC2R, FFTW_MEASURE|FFTW_DESTROY_INPUT);
  // assert(fftw_plan_bwd_);
  for(int ichan = 0; ichan<run_config->configured_channel_id_size(); ichan++)
  {
    auto* hg_wf = results_.add_high_gain();
    hg_wf->mutable_psd_sum()->Resize(nfreq,0);
    if(has_dual_gain_) {
      auto* lg_wf = results_.add_low_gain();
      lg_wf->mutable_psd_sum()->Resize(nfreq,0);
    }
  }

  return true;
}

bool WaveformPSDParallelVisitor::leave_telescope_run(
  calin::ix::provenance::chronicle::ProcessingRecord* processing_record)
{
  run_config_ = nullptr;
  // fftwf_destroy_plan(fftw_plan_bwd_);
  // fftw_plan_bwd_ = nullptr;
  fftwf_destroy_plan(fftw_plan_fwd_);
  fftw_plan_fwd_ = nullptr;
  fftwf_free(waveform_t_);
  waveform_t_ = nullptr;
  fftwf_free(waveform_f_);
  waveform_f_ = nullptr;
  return true;
}

bool WaveformPSDParallelVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  const int nchan = run_config_->configured_channel_id_size();
  const int nsamp = run_config_->num_samples();

  if(event->has_high_gain_image() and event->high_gain_image().has_camera_waveforms())
  {
    const ix::iact_data::telescope_event::Waveforms* wf =
      &event->high_gain_image().camera_waveforms();
    const uint16_t*__restrict__ wf_data = reinterpret_cast<const uint16_t*__restrict__>(
      wf->raw_samples_array().data());
    for(int ichan = 0; ichan<nchan; ichan++) {
      if(wf->channel_signal_type(ichan) != calin::ix::iact_data::telescope_event::SIGNAL_NONE) {
        process_one_waveform(wf_data, results_.mutable_high_gain(ichan));
      }
      wf_data += nsamp;
    }
  }

  if(event->has_low_gain_image() and event->low_gain_image().has_camera_waveforms())
  {
    const ix::iact_data::telescope_event::Waveforms* wf =
      &event->low_gain_image().camera_waveforms();
    const uint16_t*__restrict__ wf_data = reinterpret_cast<const uint16_t*__restrict__>(
      wf->raw_samples_array().data());
    for(int ichan = 0; ichan<nchan; ichan++) {
      if(wf->channel_signal_type(ichan) != calin::ix::iact_data::telescope_event::SIGNAL_NONE) {
        process_one_waveform(wf_data, results_.mutable_low_gain(ichan));
      }
      wf_data += nsamp;
    }
  }

  if(event->has_image() and event->image().has_camera_waveforms())
  {
    const ix::iact_data::telescope_event::Waveforms* wf =
      &event->image().camera_waveforms();
    const uint16_t*__restrict__ wf_data = reinterpret_cast<const uint16_t*__restrict__>(
      wf->raw_samples_array().data());
    for(int ichan = 0; ichan<nchan; ichan++) {
      switch(wf->channel_signal_type(ichan)) {
      case calin::ix::iact_data::telescope_event::SIGNAL_UNIQUE_GAIN:
      case calin::ix::iact_data::telescope_event::SIGNAL_HIGH_GAIN:
        process_one_waveform(wf_data, results_.mutable_high_gain(ichan));
        break;
      case calin::ix::iact_data::telescope_event::SIGNAL_LOW_GAIN:
        process_one_waveform(wf_data, results_.mutable_low_gain(ichan));
        break;
      default:
        break;
      }
      wf_data += nsamp;
    }
  }
  return true;
}

void WaveformPSDParallelVisitor::
process_one_waveform(const uint16_t*__restrict__ wf,
  calin::ix::diagnostics::waveform::WaveformSumPSD* psd)
{
  const unsigned nsample = run_config_->num_samples();
  int32_t wf_sum = 0;
  for(unsigned isample=0; isample<nsample; isample++) {
    wf_sum += wf[isample];
  }
  float wf_mean = float(wf_sum)/float(nsample);
  for(unsigned isample=0; isample<nsample; isample++) {
    waveform_t_[isample] = float(wf[isample]) - wf_mean;
  }
  fftwf_execute(fftw_plan_fwd_);
  psd->increment_num_entries();
  auto* psd_sum = psd->mutable_psd_sum()->mutable_data();
  const float*__restrict__ ri = waveform_f_;
  const float*__restrict__ ci = waveform_f_ + nsample-1;
  double psdi = SQR(double(*ri++) + wf_sum);
  (*psd_sum++) += psdi;
  while(ri < ci)
  {
    psdi = SQR(double(*ri++)) + SQR(double(*ci--));
    (*psd_sum++) += psdi;
  }
  if(ri==ci)
  {
    psdi = SQR(double(*ri));
    *psd_sum += psdi;
  }
}

bool WaveformPSDParallelVisitor::merge_results()
{
  if(parent_) {
    parent_->results_.IntegrateFrom(results_);
  }
  return true;
}
