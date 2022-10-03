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
#include <iact_data/instrument_layout.pb.h>
#include <math/fftw_util.hpp>
#include <math/rng_vcl.hpp>
#include <util/log.hpp>

namespace calin { namespace simulation { namespace waveform_processor {

class WaveformProcessor
{
public:
  WaveformProcessor(unsigned npixels, double trace_sampling_ns, unsigned trace_nsamples,
    double trace_advance_time, calin::math::rng::RNG* rng = nullptr,
    unsigned fftw_flags = 0, bool adopt_rng = false);
  WaveformProcessor(const calin::ix::iact_data::instrument_layout::CameraLayout* camera,
    double trace_sampling_ns, unsigned trace_nsamples,
    double trace_advance_time, calin::math::rng::RNG* rng = nullptr,
    unsigned fftw_flags = 0, bool adopt_rng = false);
  ~WaveformProcessor();
  void load_pes_from_processor(
    calin::simulation::pe_processor::UnbinnedWaveformPEProcessor* pe_processor,
    unsigned iscope);
  void add_nsb(double nsb_rate_ghz,
    calin::simulation::detector_efficiency::PEAmplitudeGenerator* nsb_pegen = nullptr,
    bool ac_couple=true);
  void convolve_unity_impulse_response(double pedestal = 0);
#ifndef SWIG
  void convolve_impulse_response(const double* impulse_response_dft, double pedestal = 0);
  void add_electronics_noise(const double* noise_spectrum_amplitude, double scale = 1.0);
#endif
  void convolve_impulse_response(const Eigen::VectorXd& impulse_response_dft, double pedestal = 0);
  void add_electronics_noise(const Eigen::VectorXd& noise_spectrum_amplitude, double scale = 1.0);
  void clear_pes();

  Eigen::MatrixXd pe_waveform() const;
  Eigen::MatrixXd el_waveform();

  int digital_multiplicity_trigger(double threshold,
    unsigned time_over_threshold_samples, unsigned coherence_time_samples,
    unsigned multiplicity_threshold, bool loud = false);
  int digital_multiplicity_trigger_alt(double threshold,
    unsigned time_over_threshold_samples, unsigned coherence_time_samples,
    unsigned multiplicity_threshold, bool loud = false);

  int digital_nn_trigger(double threshold,
    unsigned time_over_threshold_samples, unsigned coherence_time_samples,
    unsigned nn_threshold);
  int digital_nn_trigger_alt(double threshold,
    unsigned time_over_threshold_samples, unsigned coherence_time_samples,
    unsigned nn_threshold);

  double wavewform_t0() { return wavewform_t0_; }
  double ac_coupling_constant() { return ac_coupling_constant_; }

#ifndef SWIG
  template<typename VCLArchitecture> void vcl_add_nsb(
    calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng, double nsb_rate_ghz,
    calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
    bool ac_couple=true);

  template<typename VCLArchitecture> void vcl_add_nsb(
    calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_a,
    calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_b,
    double nsb_rate_ghz,
    calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
    bool ac_couple=true);

  template<typename VCLArchitecture> void vcl_add_nsb(
    calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_a,
    calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_b,
    calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_c,
    double nsb_rate_ghz,
    calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
    bool ac_couple=true);

  template<typename VCLArchitecture> void vcl_add_nsb(
    calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_a,
    calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_b,
    calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_c,
    calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_d,
    double nsb_rate_ghz,
    calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
    bool ac_couple=true);

  template<typename VCLArchitecture> void vcl_add_electronics_noise(
    calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng,
    const double* noise_spectrum_amplitude, double scale = 1.0);

  template<typename VCLArchitecture> void vcl_add_electronics_noise(
    calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_a,
    calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_b,
    const double* noise_spectrum_amplitude, double scale = 1.0);

  template<typename VCLArchitecture> int vcl_digital_multiplicity_trigger_alt(
      double threshold,
      unsigned time_over_threshold_samples, unsigned coherence_time_samples,
      unsigned multiplicity_threshold, bool loud = false);

#endif // SWIG

  void vcl128_add_nsb(calin::math::rng::VCLRNG<calin::util::vcl::VCL128Architecture>& vcl_rng, double nsb_rate_ghz,
      calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
      bool ac_couple=true);
  void vcl256_add_nsb(calin::math::rng::VCLRNG<calin::util::vcl::VCL256Architecture>& vcl_rng, double nsb_rate_ghz,
      calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
      bool ac_couple=true);
  void vcl512_add_nsb(calin::math::rng::VCLRNG<calin::util::vcl::VCL512Architecture>& vcl_rng, double nsb_rate_ghz,
      calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
      bool ac_couple=true);

  void vcl128_add_nsb(calin::math::rng::VCLRNG<calin::util::vcl::VCL128Architecture>& vcl_rng_a,
      calin::math::rng::VCLRNG<calin::util::vcl::VCL128Architecture>& vcl_rng_b, double nsb_rate_ghz,
      calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
      bool ac_couple=true);
  void vcl256_add_nsb(calin::math::rng::VCLRNG<calin::util::vcl::VCL256Architecture>& vcl_rng_a,
      calin::math::rng::VCLRNG<calin::util::vcl::VCL256Architecture>& vcl_rng_b, double nsb_rate_ghz,
      calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
      bool ac_couple=true);
  void vcl512_add_nsb(calin::math::rng::VCLRNG<calin::util::vcl::VCL512Architecture>& vcl_rng_a,
      calin::math::rng::VCLRNG<calin::util::vcl::VCL512Architecture>& vcl_rng_b, double nsb_rate_ghz,
      calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
      bool ac_couple=true);

  void vcl128_add_nsb(calin::math::rng::VCLRNG<calin::util::vcl::VCL128Architecture>& vcl_rng_a,
      calin::math::rng::VCLRNG<calin::util::vcl::VCL128Architecture>& vcl_rng_b,
      calin::math::rng::VCLRNG<calin::util::vcl::VCL128Architecture>& vcl_rng_c, double nsb_rate_ghz,
      calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
      bool ac_couple=true);
  void vcl256_add_nsb(calin::math::rng::VCLRNG<calin::util::vcl::VCL256Architecture>& vcl_rng_a,
      calin::math::rng::VCLRNG<calin::util::vcl::VCL256Architecture>& vcl_rng_b,
      calin::math::rng::VCLRNG<calin::util::vcl::VCL256Architecture>& vcl_rng_c, double nsb_rate_ghz,
      calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
      bool ac_couple=true);
  void vcl512_add_nsb(calin::math::rng::VCLRNG<calin::util::vcl::VCL512Architecture>& vcl_rng_a,
      calin::math::rng::VCLRNG<calin::util::vcl::VCL512Architecture>& vcl_rng_b,
      calin::math::rng::VCLRNG<calin::util::vcl::VCL512Architecture>& vcl_rng_c, double nsb_rate_ghz,
      calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
      bool ac_couple=true);

  void vcl128_add_nsb(calin::math::rng::VCLRNG<calin::util::vcl::VCL128Architecture>& vcl_rng_a,
      calin::math::rng::VCLRNG<calin::util::vcl::VCL128Architecture>& vcl_rng_b,
      calin::math::rng::VCLRNG<calin::util::vcl::VCL128Architecture>& vcl_rng_c,
      calin::math::rng::VCLRNG<calin::util::vcl::VCL128Architecture>& vcl_rng_d, double nsb_rate_ghz,
      calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
      bool ac_couple=true);
  void vcl256_add_nsb(calin::math::rng::VCLRNG<calin::util::vcl::VCL256Architecture>& vcl_rng_a,
      calin::math::rng::VCLRNG<calin::util::vcl::VCL256Architecture>& vcl_rng_b,
      calin::math::rng::VCLRNG<calin::util::vcl::VCL256Architecture>& vcl_rng_c,
      calin::math::rng::VCLRNG<calin::util::vcl::VCL256Architecture>& vcl_rng_d, double nsb_rate_ghz,
      calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
      bool ac_couple=true);
  void vcl512_add_nsb(calin::math::rng::VCLRNG<calin::util::vcl::VCL512Architecture>& vcl_rng_a,
      calin::math::rng::VCLRNG<calin::util::vcl::VCL512Architecture>& vcl_rng_b,
      calin::math::rng::VCLRNG<calin::util::vcl::VCL512Architecture>& vcl_rng_c,
      calin::math::rng::VCLRNG<calin::util::vcl::VCL512Architecture>& vcl_rng_d, double nsb_rate_ghz,
      calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
      bool ac_couple=true);

  void vcl128_add_electronics_noise(
    calin::math::rng::VCLRNG<calin::util::vcl::VCL128Architecture>& vcl_rng,
    const Eigen::VectorXd& noise_spectrum_amplitude, double scale = 1.0);
  void vcl256_add_electronics_noise(
    calin::math::rng::VCLRNG<calin::util::vcl::VCL256Architecture>& vcl_rng,
    const Eigen::VectorXd& noise_spectrum_amplitude, double scale = 1.0);
  void vcl512_add_electronics_noise(
    calin::math::rng::VCLRNG<calin::util::vcl::VCL512Architecture>& vcl_rng,
    const Eigen::VectorXd& noise_spectrum_amplitude, double scale = 1.0);

  void vcl128_add_electronics_noise(
    calin::math::rng::VCLRNG<calin::util::vcl::VCL128Architecture>& vcl_rng_a,
    calin::math::rng::VCLRNG<calin::util::vcl::VCL128Architecture>& vcl_rng_b,
    const Eigen::VectorXd& noise_spectrum_amplitude, double scale = 1.0);
  void vcl256_add_electronics_noise(
    calin::math::rng::VCLRNG<calin::util::vcl::VCL256Architecture>& vcl_rng_a,
    calin::math::rng::VCLRNG<calin::util::vcl::VCL256Architecture>& vcl_rng_b,
    const Eigen::VectorXd& noise_spectrum_amplitude, double scale = 1.0);
  void vcl512_add_electronics_noise(
    calin::math::rng::VCLRNG<calin::util::vcl::VCL512Architecture>& vcl_rng_a,
    calin::math::rng::VCLRNG<calin::util::vcl::VCL512Architecture>& vcl_rng_b,
    const Eigen::VectorXd& noise_spectrum_amplitude, double scale = 1.0);

  int vcl256_digital_multiplicity_trigger_alt(double threshold,
      unsigned time_over_threshold_samples, unsigned coherence_time_samples,
      unsigned multiplicity_threshold, bool loud = false);

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

  int* neighbour_map_ = nullptr;
  int max_num_neighbours_ = 0;
};

#include "waveform_processor_vcl_function_implemntations.hpp"

} } } // namespace calin::simulation::waveform_processor
