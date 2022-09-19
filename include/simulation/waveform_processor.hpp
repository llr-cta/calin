/*

   calin/simulation/waveform_processor.hpp -- Stephen Fegan -- 2022-09-16

   Class to process waveforms

   Copyright 2022, Stephen Fegan <sfegan@llr.in2p3.fr>
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

// Units:
// Time, t:      ns
// Height, z:    cm
// Energy, e:    MeV

#include <vector>
#include <ostream>
#include <map>

#include <fftw3.h>

#include <Eigen/Core>
#include <simulation/pe_processor.hpp>
#include <math/fftw_util.hpp>

namespace calin { namespace simulation { namespace waveform_processor {

class WaveformProcessor
{
public:
  WaveformProcessor(unsigned npixels, double trace_sampling_ns, unsigned trace_nsamples,
    double trace_advance_time, calin::math::rng::RNG* rng = nullptr,
    unsigned fftw_flags = 0, bool adopt_rng = false);
  ~WaveformProcessor();
  void load_pes_from_processor(
    calin::simulation::pe_processor::UnbinnedWaveformPEProcessor* pe_processor,
    unsigned iscope, double trace_advance_time);
  void add_nsb(double nsb_rate_ghz,
    calin::simulation::detector_efficiency::PEAmplitudeGenerator* nsb_pegen = nullptr,
    bool ac_couple=true);
  void convolve_unity_impulse_response(double pedestal = 0);
  void convolve_impulse_response(const double* impulse_response_dft, double pedestal = 0);
  void add_electronics_noise(const double* noise_spectrum_amplitude);
  void clear_pes();
  Eigen::MatrixXd pe_waveform() const;
  Eigen::MatrixXd el_waveform() const;
  double wavewform_t0() { return wavewform_t0_; }
private:
  void compute_pe_waveform_dft();
  void compute_el_waveform();
  unsigned npixels_;
  double trace_sampling_ns_;
  double trace_sampling_inv_;
  unsigned trace_nsamples_;
  double trace_advance_time_;
  calin::math::rng::RNG* rng_ = nullptr;
  bool adopt_rng_;

  double wavewform_t0_ = std::numeric_limits<double>::quiet_NaN();
  double* pe_waveform_ = nullptr;
  double* pe_waveform_dft_ = nullptr;
  double* el_waveform_dft_ = nullptr;
  double* el_waveform_ = nullptr;
  fftw_plan fwd_plan_;
  fftw_plan rev_plan_;
  bool pe_waveform_dft_valid_ = false;
  bool el_waveform_dft_valid_ = false;
  bool el_waveform_valid_ = false;
  double ac_coupling_constant_ = 0;
};

} } } // namespace calin::simulation::waveform_processor
