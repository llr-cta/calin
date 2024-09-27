  /*

   calin/diagnostics/waveform_psd_vcl.hpp -- Stephen Fegan -- 2021-09-27

   Waveform PSD vectorized diagnostics visitor

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

#pragma once

#include <math/fftw_util.hpp>

#include <iact_data/event_visitor.hpp>
#include <diagnostics/waveform.pb.h>
#include <diagnostics/waveform.hpp>
#include <math/special.hpp>
#include <util/log.hpp>
#include <util/memory.hpp>

using namespace calin::util::log;

namespace calin { namespace diagnostics { namespace waveform {

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCL_WaveformPSDParallelVisitor:
  public WaveformPSDParallelVisitor
{
public:
  VCL_WaveformPSDParallelVisitor(): WaveformPSDParallelVisitor()
  {
    // nothing to see here
  }

  virtual ~VCL_WaveformPSDParallelVisitor()
  {
    ::free(samples_);
    ::free(waveform_x_);
    ::free(waveform_f_);
    ::free(waveform_p_);
    ::free(psd_count_hg_);
    ::free(psd_count_lg_);
    ::free(dc_sum_hg_);
    ::free(dc_sum_lg_);
    ::free(psd_sum_hg_);
    ::free(psd_sum_lg_);
  }

  VCL_WaveformPSDParallelVisitor* new_sub_visitor(
    std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
        calin::iact_data::event_visitor::ParallelEventVisitor*>
      antecedent_visitors) override
  {
    auto* sub_visitor = new VCL_WaveformPSDParallelVisitor();
    sub_visitor->parent_ = this;
    return sub_visitor;
  }

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
    calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager,
    calin::ix::provenance::chronicle::ProcessingRecord* processing_record = nullptr) override
  {
    has_codelet_ = codelet_.has_codelet(run_config->num_samples());
    if(has_codelet_) {
      if(processing_record) {
        if(processing_record->type().empty()) {
          processing_record->set_type(calin::util::vcl::templated_class_name<VCLArchitecture>("VCL_WaveformPSDParallelVisitor"));
        }
        // Rest of the entries are set by WaveformPSDParallelVisitor::visit_telescope_run
      }

      static constexpr unsigned num_int16 = VCLArchitecture::num_int16;

      nchan_ = run_config->configured_channel_id_size();
      nsamp_ = run_config->num_samples();
      nfreq_ = calin::math::fftw_util::hcvec_num_real(nsamp_);
      nchan_block_ = (nchan_+num_int16-1)/num_int16;
      nsamp_block_ = (nsamp_+num_int16-1)/num_int16;
      nsamp_inv_ = 1.0/nsamp_;

      calin::util::memory::safe_aligned_recalloc_and_fill(signal_type_, nchan_);
      calin::util::memory::safe_aligned_recalloc_and_fill(samples_, nsamp_block_*num_int16);
      calin::util::memory::safe_aligned_recalloc_and_fill(waveform_x_, nsamp_);
      calin::util::memory::safe_aligned_recalloc_and_fill(waveform_f_, nsamp_);
      calin::util::memory::safe_aligned_recalloc_and_fill(waveform_p_, nfreq_);
      calin::util::memory::safe_aligned_recalloc_and_fill(psd_count_hg_, nchan_block_*num_int16);
      calin::util::memory::safe_aligned_recalloc_and_fill(dc_sum_hg_, nchan_block_*num_int16);
      calin::util::memory::safe_aligned_recalloc_and_fill(psd_sum_hg_, nchan_block_*num_int16 * nfreq_);
      if(run_config->camera_layout().adc_gains() !=
          calin::ix::iact_data::instrument_layout::CameraLayout::SINGLE_GAIN) {
        calin::util::memory::safe_aligned_recalloc_and_fill(psd_count_lg_, nchan_block_*num_int16);
        calin::util::memory::safe_aligned_recalloc_and_fill(dc_sum_lg_, nchan_block_*num_int16);
        calin::util::memory::safe_aligned_recalloc_and_fill(psd_sum_lg_, nchan_block_*num_int16 * nfreq_);
      }
    } else {
      LOG(WARNING) << "VCL_WaveformPSDParallelVisitor : no codelet for array size : "
        << run_config->num_samples() << '\n'
        << "Falling back to non-vectorized implementation.";
    }

    return WaveformPSDParallelVisitor::visit_telescope_run(run_config, event_lifetime_manager, processing_record);
  }

  bool leave_telescope_run(
    calin::ix::provenance::chronicle::ProcessingRecord* processing_record = nullptr) override
  {
    static constexpr unsigned num_double = VCLArchitecture::num_double;
    if(has_codelet_) {
      for(unsigned ichan=0; ichan<nchan_; ++ichan) {
        unsigned ichan_block = ichan / num_double;
        unsigned ichan_offset = ichan - ichan_block*num_double;
        auto* chan_results = results_.mutable_high_gain(ichan);
        chan_results->set_num_entries(chan_results->num_entries() + psd_count_hg_[ichan]);
        chan_results->set_dc_sum(chan_results->dc_sum() + dc_sum_hg_[ichan]);
        for(unsigned ifreq=0; ifreq<nfreq_; ++ifreq) {
          chan_results->set_psd_sum(ifreq, chan_results->psd_sum(ifreq) +
            psd_sum_hg_[(ichan_block*nfreq_ + ifreq) * num_double + ichan_offset]);
        }
      }

      if(psd_sum_lg_ != nullptr) {
        for(unsigned ichan=0; ichan<nchan_; ++ichan) {
          unsigned ichan_block = ichan / num_double;
          unsigned ichan_offset = ichan - ichan_block*num_double;
          auto* chan_results = results_.mutable_low_gain(ichan);
          chan_results->set_num_entries(chan_results->num_entries() + psd_count_lg_[ichan]);
          chan_results->set_dc_sum(chan_results->dc_sum() + dc_sum_lg_[ichan]);
          for(unsigned ifreq=0; ifreq<nfreq_; ++ifreq) {
            chan_results->set_psd_sum(ifreq, chan_results->psd_sum(ifreq) +
              psd_sum_lg_[(ichan_block*nfreq_ + ifreq) * num_double + ichan_offset]);
          }
        }
      }

      ::free(signal_type_);
      signal_type_ = nullptr;
      ::free(samples_);
      samples_ = nullptr;

      ::free(waveform_x_);
      waveform_x_ = nullptr;
      ::free(waveform_f_);
      waveform_f_ = nullptr;
      ::free(waveform_p_);
      waveform_p_ = nullptr;

      ::free(psd_count_hg_);
      psd_count_hg_ = nullptr;
      ::free(psd_count_lg_);
      psd_count_lg_ = nullptr;

      ::free(dc_sum_hg_);
      dc_sum_hg_ = nullptr;
      ::free(dc_sum_lg_);
      dc_sum_lg_ = nullptr;

      ::free(psd_sum_hg_);
      psd_sum_hg_ = nullptr;
      ::free(psd_sum_lg_);
      psd_sum_lg_ = nullptr;
    }

    return WaveformPSDParallelVisitor::leave_telescope_run(processing_record);
  }

  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override
  {
    if(not has_codelet_) {
      return WaveformPSDParallelVisitor::visit_telescope_event(seq_index, event);
    }

    if(event->has_high_gain_image() and event->high_gain_image().has_camera_waveforms()) {
      const ix::iact_data::telescope_event::Waveforms* wf =
        &event->high_gain_image().camera_waveforms();
      const uint16_t*__restrict__ wf_data = reinterpret_cast<const uint16_t*__restrict__>(
        wf->raw_samples_array().data());
      for(unsigned ichan = 0; ichan<nchan_; ichan++) {
        signal_type_[ichan] = wf->channel_signal_type(ichan);
      }
      vcl_analyze_waveforms(wf_data, signal_type_);
    }

    if(event->has_low_gain_image() and event->low_gain_image().has_camera_waveforms()) {
      const ix::iact_data::telescope_event::Waveforms* wf =
        &event->low_gain_image().camera_waveforms();
      const uint16_t*__restrict__ wf_data = reinterpret_cast<const uint16_t*__restrict__>(
        wf->raw_samples_array().data());
      for(unsigned ichan = 0; ichan<nchan_; ichan++) {
        signal_type_[ichan] = wf->channel_signal_type(ichan);
      }
      vcl_analyze_waveforms(wf_data, signal_type_);
    }

    if(event->has_image() and event->image().has_camera_waveforms()) {
      const ix::iact_data::telescope_event::Waveforms* wf =
        &event->image().camera_waveforms();
      const uint16_t*__restrict__ wf_data = reinterpret_cast<const uint16_t*__restrict__>(
        wf->raw_samples_array().data());
      for(unsigned ichan = 0; ichan<nchan_; ichan++) {
        signal_type_[ichan] = wf->channel_signal_type(ichan);
      }
      vcl_analyze_waveforms(wf_data, signal_type_);
    }

    return true;
  }

#ifndef SWIG
protected:
  void vcl_analyze_waveforms(const uint16_t* __restrict__ data, const int32_t*__restrict__ signal_type)
  {
    typedef typename VCLArchitecture::int32_vt int32_vt;
    typedef typename VCLArchitecture::int64_vt int64_vt;
    typedef typename VCLArchitecture::int64_bvt int64_bvt;
    typedef typename VCLArchitecture::float_vt float_vt;
    typedef typename VCLArchitecture::double_vt double_vt;

    static constexpr unsigned num_int16 = VCLArchitecture::num_int16;
    static constexpr unsigned num_int32 = VCLArchitecture::num_int32;
    static constexpr unsigned num_int64 = VCLArchitecture::num_int64;
    static constexpr unsigned num_double = VCLArchitecture::num_double;

    using calin::math::special::SQR;

    unsigned nchan_left = nchan_;

    for(unsigned ichan_block=0;ichan_block<nchan_block_;++ichan_block)
    {
      unsigned nchan_proc = std::min(num_int16, nchan_left);
      nchan_left -= nchan_proc;

      for(unsigned isamp_block=0;isamp_block<nsamp_block_;++isamp_block)
      {
        for(unsigned ichan_proc=0;ichan_proc<nchan_proc;++ichan_proc) {
          samples_[isamp_block*num_int16 + ichan_proc].load(data
            + (ichan_block*num_int16 + ichan_proc)*nsamp_
            + isamp_block*num_int16);
        }
        calin::util::vcl::transpose(samples_ + isamp_block*num_int16);
      }

      int32_vt q_sum;
      float_vt q_mean;
      const float_vt*__restrict__ ri;
      const float_vt*__restrict__ ci;
      double_vt*__restrict__ pi;
      int32_vt type_32;
      int64_vt type_64;
      int64_bvt mask_hg;
      int64_bvt mask_lg;
      double_vt dc_val;

      // -----------------------------------------------------------------------
      // PASS 1 : 16-bit -> 32-bit for values [..XX]
      // -----------------------------------------------------------------------

      q_sum = 0;

      for(unsigned isamp=0; isamp<nsamp_; ++isamp) {
        int32_vt q = calin::util::vcl::extend_16_to_32_low(samples_[isamp]);
        q_sum += q;
      }

      q_mean = calin::util::vcl::to_float(q_sum) / float(nsamp_);

      for(unsigned isamp=0; isamp<nsamp_; ++isamp) {
        int32_vt q = calin::util::vcl::extend_16_to_32_low(samples_[isamp]);
        waveform_x_[isamp] = calin::util::vcl::to_float(q) - q_mean;
      }

      codelet_.r2hc(nsamp_, waveform_x_, waveform_f_);

      type_32.load(signal_type + ichan_block*num_int16);

      // -----------------------------------------------------------------------
      // PASS 1A : 32-bit -> 64-bit for values [...X]
      // -----------------------------------------------------------------------

      type_64 = vcl::extend_low(type_32);

      mask_hg = (type_64==1) | (type_64==2);
      mask_lg = (type_64==3);

      ri = waveform_f_;
      ci = waveform_f_ + nsamp_ - 1;
      pi = waveform_p_;

      dc_val = vcl::extend_low(*ri++) + vcl::to_double(vcl::extend_low(q_sum));
      *(pi++) = SQR(dc_val);
      while(ri < ci) {
        *(pi++) = SQR(vcl::extend_low(*ri++)) + SQR(vcl::extend_low(*ci--));
      }
      if(ri==ci) {
        *(pi++) = SQR(vcl::extend_low(*ri++));
      }

#define COUNT_INDEX(array, block, subblock) \
  (array + num_int64*(block*4 + subblock))
#define SUM_INDEX(array, block, subblock, freq) \
  (array + num_double*((block*4 + subblock)*nfreq_ + ifreq))

      if(vcl::horizontal_or(mask_hg)) {
        int64_vt count;
        count.load(COUNT_INDEX(psd_count_hg_, ichan_block, 0));
        count = vcl::if_add(mask_hg, count, 1);
        count.store(COUNT_INDEX(psd_count_hg_, ichan_block, 0));

        double_vt dc_sum;
        dc_sum.load(COUNT_INDEX(dc_sum_hg_, ichan_block, 0));
        dc_sum = vcl::if_add(mask_hg, dc_sum, dc_val);
        dc_sum.store(COUNT_INDEX(dc_sum_hg_, ichan_block, 0));

        for(unsigned ifreq=0; ifreq<nfreq_; ifreq++) {
          double_vt psd_sum;
          psd_sum.load(SUM_INDEX(psd_sum_hg_, ichan_block, 0, ifreq));
          psd_sum = vcl::if_add(mask_hg, psd_sum, waveform_p_[ifreq]);
          psd_sum.store(SUM_INDEX(psd_sum_hg_, ichan_block, 0, ifreq));
        }
      }

      if(vcl::horizontal_or(mask_lg)) {
        int64_vt count;
        count.load(COUNT_INDEX(psd_count_lg_, ichan_block, 0));
        count = vcl::if_add(mask_lg, count, 1);
        count.store(COUNT_INDEX(psd_count_lg_, ichan_block, 0));

        double_vt dc_sum;
        dc_sum.load(COUNT_INDEX(dc_sum_lg_, ichan_block, 0));
        dc_sum = vcl::if_add(mask_lg, dc_sum, dc_val);
        dc_sum.store(COUNT_INDEX(dc_sum_lg_, ichan_block, 0));

        for(unsigned ifreq=0; ifreq<nfreq_; ifreq++) {
          double_vt psd_sum;
          psd_sum.load(SUM_INDEX(psd_sum_lg_, ichan_block, 0, ifreq));
          psd_sum = vcl::if_add(mask_lg, psd_sum, waveform_p_[ifreq]);
          psd_sum.store(SUM_INDEX(psd_sum_lg_, ichan_block, 0, ifreq));
        }
      }

      // -----------------------------------------------------------------------
      // PASS 1B : 32-bit -> 64-bit for values [..X.]
      // -----------------------------------------------------------------------

      type_64 = vcl::extend_high(type_32);

      mask_hg = (type_64==1) | (type_64==2);
      mask_lg = (type_64==3);

      ri = waveform_f_;
      ci = waveform_f_ + nsamp_ - 1;
      pi = waveform_p_;

      dc_val = vcl::extend_high(*ri++) + vcl::to_double(vcl::extend_high(q_sum));
      *(pi++) = SQR(dc_val);
      while(ri < ci) {
        *(pi++) = SQR(vcl::extend_high(*ri++)) + SQR(vcl::extend_high(*ci--));
      }
      if(ri==ci) {
        *(pi++) = SQR(vcl::extend_high(*ri++));
      }

      if(vcl::horizontal_or(mask_hg)) {
        int64_vt count;
        count.load(COUNT_INDEX(psd_count_hg_, ichan_block, 1));
        count = vcl::if_add(mask_hg, count, 1);
        count.store(COUNT_INDEX(psd_count_hg_, ichan_block, 1));

        double_vt dc_sum;
        dc_sum.load(COUNT_INDEX(dc_sum_hg_, ichan_block, 1));
        dc_sum = vcl::if_add(mask_hg, dc_sum, dc_val);
        dc_sum.store(COUNT_INDEX(dc_sum_hg_, ichan_block, 1));

        for(unsigned ifreq=0; ifreq<nfreq_; ifreq++) {
          double_vt psd_sum;
          psd_sum.load(SUM_INDEX(psd_sum_hg_, ichan_block, 1, ifreq));
          psd_sum = vcl::if_add(mask_hg, psd_sum, waveform_p_[ifreq]);
          psd_sum.store(SUM_INDEX(psd_sum_hg_, ichan_block, 1, ifreq));
        }
      }

      if(vcl::horizontal_or(mask_lg)) {
        int64_vt count;
        count.load(COUNT_INDEX(psd_count_lg_, ichan_block, 1));
        count = vcl::if_add(mask_lg, count, 1);
        count.store(COUNT_INDEX(psd_count_lg_, ichan_block, 1));

        double_vt dc_sum;
        dc_sum.load(COUNT_INDEX(dc_sum_lg_, ichan_block, 1));
        dc_sum = vcl::if_add(mask_lg, dc_sum, dc_val);
        dc_sum.store(COUNT_INDEX(dc_sum_lg_, ichan_block, 1));

        for(unsigned ifreq=0; ifreq<nfreq_; ifreq++) {
          double_vt psd_sum;
          psd_sum.load(SUM_INDEX(psd_sum_lg_, ichan_block, 1, ifreq));
          psd_sum = vcl::if_add(mask_lg, psd_sum, waveform_p_[ifreq]);
          psd_sum.store(SUM_INDEX(psd_sum_lg_, ichan_block, 1, ifreq));
        }
      }

      // -----------------------------------------------------------------------
      // PASS 2 : 16-bit -> 32-bit for values [XX..]
      // -----------------------------------------------------------------------

      q_sum = 0;

      for(unsigned isamp=0; isamp<nsamp_; ++isamp) {
        int32_vt q = calin::util::vcl::extend_16_to_32_high(samples_[isamp]);
        q_sum += q;
      }

      q_mean = calin::util::vcl::to_float(q_sum) / float(nsamp_);

      for(unsigned isamp=0; isamp<nsamp_; ++isamp) {
        int32_vt q = calin::util::vcl::extend_16_to_32_high(samples_[isamp]);
        waveform_x_[isamp] = calin::util::vcl::to_float(q) - q_mean;
      }

      codelet_.r2hc(nsamp_, waveform_x_, waveform_f_);

      type_32.load(signal_type + ichan_block*num_int16 + num_int32);

      // -----------------------------------------------------------------------
      // PASS 2A : 32-bit -> 64-bit for values [.X..]
      // -----------------------------------------------------------------------

      type_64 = vcl::extend_low(type_32);

      mask_hg = (type_64==1) | (type_64==2);
      mask_lg = (type_64==3);

      ri = waveform_f_;
      ci = waveform_f_ + nsamp_ - 1;
      pi = waveform_p_;

      dc_val = vcl::extend_low(*ri++) + vcl::to_double(vcl::extend_low(q_sum));
      *(pi++) = SQR(dc_val);
      while(ri < ci) {
        *(pi++) = SQR(vcl::extend_low(*ri++)) + SQR(vcl::extend_low(*ci--));
      }
      if(ri==ci) {
        *(pi++) = SQR(vcl::extend_low(*ri++));
      }

      if(vcl::horizontal_or(mask_hg)) {
        int64_vt count;
        count.load(COUNT_INDEX(psd_count_hg_, ichan_block, 2));
        count = vcl::if_add(mask_hg, count, 1);
        count.store(COUNT_INDEX(psd_count_hg_, ichan_block, 2));

        double_vt dc_sum;
        dc_sum.load(COUNT_INDEX(dc_sum_hg_, ichan_block, 2));
        dc_sum = vcl::if_add(mask_hg, dc_sum, dc_val);
        dc_sum.store(COUNT_INDEX(dc_sum_hg_, ichan_block, 2));

        for(unsigned ifreq=0; ifreq<nfreq_; ifreq++) {
          double_vt psd_sum;
          psd_sum.load(SUM_INDEX(psd_sum_hg_, ichan_block, 2, ifreq));
          psd_sum = vcl::if_add(mask_hg, psd_sum, waveform_p_[ifreq]);
          psd_sum.store(SUM_INDEX(psd_sum_hg_, ichan_block, 2, ifreq));
        }
      }

      if(vcl::horizontal_or(mask_lg)) {
        int64_vt count;
        count.load(COUNT_INDEX(psd_count_lg_, ichan_block, 2));
        count = vcl::if_add(mask_lg, count, 1);
        count.store(COUNT_INDEX(psd_count_lg_, ichan_block, 2));

        double_vt dc_sum;
        dc_sum.load(COUNT_INDEX(dc_sum_lg_, ichan_block, 2));
        dc_sum = vcl::if_add(mask_lg, dc_sum, dc_val);
        dc_sum.store(COUNT_INDEX(dc_sum_lg_, ichan_block, 2));

        for(unsigned ifreq=0; ifreq<nfreq_; ifreq++) {
          double_vt psd_sum;
          psd_sum.load(SUM_INDEX(psd_sum_lg_, ichan_block, 2, ifreq));
          psd_sum = vcl::if_add(mask_lg, psd_sum, waveform_p_[ifreq]);
          psd_sum.store(SUM_INDEX(psd_sum_lg_, ichan_block, 2, ifreq));
        }
      }

      // -----------------------------------------------------------------------
      // PASS 1B : 32-bit -> 64-bit for values [X...]
      // -----------------------------------------------------------------------

      type_64 = vcl::extend_high(type_32);

      mask_hg = (type_64==1) | (type_64==2);
      mask_lg = (type_64==3);

      ri = waveform_f_;
      ci = waveform_f_ + nsamp_ - 1;
      pi = waveform_p_;

      dc_val = vcl::extend_high(*ri++) + vcl::to_double(vcl::extend_high(q_sum));
      *(pi++) = SQR(dc_val);
      while(ri < ci) {
        *(pi++) = SQR(vcl::extend_high(*ri++)) + SQR(vcl::extend_high(*ci--));
      }
      if(ri==ci) {
        *(pi++) = SQR(vcl::extend_high(*ri++));
      }

      if(vcl::horizontal_or(mask_hg)) {
        int64_vt count;
        count.load(COUNT_INDEX(psd_count_hg_, ichan_block, 3));
        count = vcl::if_add(mask_hg, count, 1);
        count.store(COUNT_INDEX(psd_count_hg_, ichan_block, 3));

        double_vt dc_sum;
        dc_sum.load(COUNT_INDEX(dc_sum_hg_, ichan_block, 3));
        dc_sum = vcl::if_add(mask_hg, dc_sum, dc_val);
        dc_sum.store(COUNT_INDEX(dc_sum_hg_, ichan_block, 3));

        for(unsigned ifreq=0; ifreq<nfreq_; ifreq++) {
          double_vt psd_sum;
          psd_sum.load(SUM_INDEX(psd_sum_hg_, ichan_block, 3, ifreq));
          psd_sum = vcl::if_add(mask_hg, psd_sum, waveform_p_[ifreq]);
          psd_sum.store(SUM_INDEX(psd_sum_hg_, ichan_block, 3, ifreq));
        }
      }

      if(vcl::horizontal_or(mask_lg)) {
        int64_vt count;
        count.load(COUNT_INDEX(psd_count_lg_, ichan_block, 3));
        count = vcl::if_add(mask_lg, count, 1);
        count.store(COUNT_INDEX(psd_count_lg_, ichan_block, 3));

        double_vt dc_sum;
        dc_sum.load(COUNT_INDEX(dc_sum_lg_, ichan_block, 3));
        dc_sum = vcl::if_add(mask_lg, dc_sum, dc_val);
        dc_sum.store(COUNT_INDEX(dc_sum_lg_, ichan_block, 3));

        for(unsigned ifreq=0; ifreq<nfreq_; ifreq++) {
          double_vt psd_sum;
          psd_sum.load(SUM_INDEX(psd_sum_lg_, ichan_block, 3, ifreq));
          psd_sum = vcl::if_add(mask_lg, psd_sum, waveform_p_[ifreq]);
          psd_sum.store(SUM_INDEX(psd_sum_lg_, ichan_block, 3, ifreq));
        }
      }
    }

#undef COUNT_INDEX
#undef SUM_INDEX
  }

  typename VCLArchitecture::int16_vt*__restrict__ samples_ = nullptr;
  typename VCLArchitecture::float_vt*__restrict__ waveform_x_ = nullptr;
  typename VCLArchitecture::float_vt*__restrict__ waveform_f_ = nullptr;
  typename VCLArchitecture::double_vt*__restrict__ waveform_p_ = nullptr;

  int32_t*__restrict__ signal_type_ = nullptr;

  int64_t*__restrict__ psd_count_hg_ = nullptr;
  int64_t*__restrict__ psd_count_lg_ = nullptr;

  double_t*__restrict__ dc_sum_hg_ = nullptr;
  double_t*__restrict__ dc_sum_lg_ = nullptr;

  double_t*__restrict__ psd_sum_hg_ = nullptr;
  double_t*__restrict__ psd_sum_lg_ = nullptr;

  unsigned nchan_ = 0;
  unsigned nsamp_ = 0;
  unsigned nfreq_ = 0;
  unsigned nchan_block_ = 0;
  unsigned nsamp_block_ = 0;
  double nsamp_inv_ = 1.0;

  calin::math::fftw_util::FFTWCodelet<calin::util::vcl::VCLFloatReal<VCLArchitecture> > codelet_;
  bool has_codelet_ = false;
#endif // defined SWIG
};

} } } // namespace calin::diagnostics::waveform
