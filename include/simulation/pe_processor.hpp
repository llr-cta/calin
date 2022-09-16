/*

   calin/simulation/pe_processor.hpp -- Stephen Fegan -- 2017-01-16

   Multi-purpose PE (weight, scope, pixel & time) processor.

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <vector>
#include <string>
#include <stdexcept>

#include <math/rng.hpp>
#include <math/accumulator.hpp>
#include <math/moments_calc.hpp>
#include <simulation/detector_efficiency.hpp>

namespace calin { namespace simulation { namespace pe_processor {

class PEProcessor
{
public:
  virtual ~PEProcessor();
  virtual void start_processing();
  virtual void process_focal_plane_hit(unsigned scope_id, int pixel_id,
    double x, double y, double ux, double uy, double t, double pe_weight);
  virtual void finish_processing();
};

class RecordingPEProcessor: public PEProcessor
{
public:
  RecordingPEProcessor(unsigned nmax = 0, bool auto_clear = true);
  virtual ~RecordingPEProcessor();
  void start_processing() override;
  void process_focal_plane_hit(unsigned scope_id, int pixel_id,
    double x, double y, double ux, double uy, double t, double pe_weight) override;
  void clear_all_pes();
  unsigned npe() const { return x_.size(); }
  unsigned scope_id(unsigned ipe) { return sid_.at(ipe); }
  unsigned pixel_id(unsigned ipe) { return pid_.at(ipe); }
  double x(unsigned ipe) { return x_.at(ipe); }
  double y(unsigned ipe) { return y_.at(ipe); }
  double ux(unsigned ipe) { return ux_.at(ipe); }
  double uz(unsigned ipe) { return uy_.at(ipe); }
  double t(unsigned ipe) { return t_.at(ipe); }
  double w(unsigned ipe) { return w_.at(ipe); }
private:
  unsigned nmax_;
  bool auto_clear_;
  std::vector<unsigned> sid_;
  std::vector<unsigned> pid_;
  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<double> ux_;
  std::vector<double> uy_;
  std::vector<double> t_;
  std::vector<double> w_;
};

class SimpleImagePEProcessor: public PEProcessor
{
public:
  SimpleImagePEProcessor(unsigned nscope, unsigned npix,
    bool auto_clear = true);
  SimpleImagePEProcessor(const std::vector<unsigned> npix,
    bool auto_clear = true);
  virtual ~SimpleImagePEProcessor();
  void start_processing() override;
  void process_focal_plane_hit(unsigned scope_id, int pixel_id,
    double x, double y, double ux, double uy, double t, double pe_weight) override;
  const std::vector<double> scope_image(unsigned iscope) const;
  void clear_all_images();
private:
  CALIN_TYPEALIAS(Accumulator, calin::math::accumulator::RecommendedAccumulator);
  bool auto_clear_ = false;
  std::vector<std::vector<Accumulator>> images_;
};

class WaveformPEProcessor: public PEProcessor
{
public:
  WaveformPEProcessor(unsigned nscope, unsigned npix, unsigned nsamp, double delta_t,
    calin::math::rng::RNG* rng = nullptr, bool auto_clear = true);
  virtual ~WaveformPEProcessor();
  void start_processing() override;
  void process_focal_plane_hit(unsigned scope_id, int pixel_id,
    double x, double y, double ux, double uy, double t, double pe_weight) override;
  const Eigen::MatrixXd& scope_traces(unsigned iscope) const { check_iscope(iscope); return traces_[iscope]; }
  double scope_t0(unsigned iscope) const { check_iscope(iscope); return t0_(iscope); }
  int scope_nmin(unsigned iscope) const { check_iscope(iscope); return nmin_(iscope); }
  int scope_nmax(unsigned iscope) const { check_iscope(iscope); return nmax_(iscope); }
  Eigen::VectorXd scope_overflow(unsigned iscope) const { check_iscope(iscope); return overflow_.row(iscope); }
  void clear_all_traces();
  void add_nsb(double rate_ghz,
    calin::simulation::detector_efficiency::PEAmplitudeGenerator* pegen = nullptr);
  void add_nsb(const Eigen::VectorXd rate_per_pixel_ghz,
    calin::simulation::detector_efficiency::PEAmplitudeGenerator* pegen = nullptr);
private:
  void check_iscope(unsigned iscope) const {
    if(iscope >= traces_.size()) {
      throw std::out_of_range("iscope out of range : " + std::to_string(iscope)
        + " >= " + std::to_string(traces_.size()));
    }
  }
  unsigned nsamp_;
  unsigned npix_;
  double delta_t_inv_;
  std::vector<Eigen::MatrixXd> traces_;
  Eigen::VectorXd t0_;
  Eigen::VectorXi nmin_;
  Eigen::VectorXi nmax_;
  Eigen::MatrixXd overflow_;
  bool auto_clear_ = false;
  bool warning_sent_ = false;
  calin::math::rng::RNG* rng_ = nullptr;
};

class UnbinnedWaveformPEProcessor: public PEProcessor
{
public:
  UnbinnedWaveformPEProcessor(unsigned nscope, unsigned npix,
    double scope_trace_delta_t = 1.0 /*ns*/, unsigned scope_trace_nsamp = 4096,
    bool auto_clear = true);
  virtual ~UnbinnedWaveformPEProcessor();
  void start_processing() override;
  void process_focal_plane_hit(unsigned scope_id, int pixel_id,
    double x, double y, double ux, double uy, double t, double pe_weight) override;
  Eigen::MatrixXd pixel_traces(double& trace_t0, Eigen::VectorXd& trace_overflow, unsigned iscope,
    double trace_delta_t, unsigned trace_nsamp, double trace_advance_time,
    double nsb_rate_ghz = 0, calin::math::rng::RNG* rng_ = nullptr,
    calin::simulation::detector_efficiency::PEAmplitudeGenerator* nsb_pegen = nullptr,
    bool ac_couple=false) const;
  Eigen::VectorXd scope_trace(unsigned iscope) const { check_iscope(iscope); return scope_trace_.row(iscope); }
  double scope_t0(unsigned iscope) const { check_iscope(iscope); return t0_(iscope); }
  int scope_nmin(unsigned iscope) const { check_iscope(iscope); return nmin_(iscope); }
  int scope_nmax(unsigned iscope) const { check_iscope(iscope); return nmax_(iscope); }
  double scope_trace_overflow(unsigned iscope) const { check_iscope(iscope); return scope_trace_overflow_(iscope); }
  void clear_all_traces();
  static Eigen::MatrixXd convolve_instrument_response(const Eigen::MatrixXd& traces,
    const Eigen::VectorXd& impulse_response_dft, double pedestal=0);
private:
  void check_iscope(unsigned iscope) const {
    if(iscope >= nscope_) {
      throw std::out_of_range("iscope out of range : " + std::to_string(iscope)
        + " >= " + std::to_string(nscope_));
    }
  }
  unsigned nscope_;
  unsigned npix_;

  unsigned scope_trace_nsamp_;
  double scope_trace_delta_t_inv_;

  std::vector<unsigned> pe_iscope_;
  std::vector<unsigned> pe_ipix_;
  std::vector<double> pe_t_;
  std::vector<double> pe_q_;

  Eigen::MatrixXd scope_trace_;
  Eigen::VectorXd t0_;
  Eigen::VectorXi nmin_;
  Eigen::VectorXi nmax_;
  Eigen::VectorXd scope_trace_overflow_;
  bool auto_clear_ = false;
  bool warning_sent_ = false;
};

class TelescopePSFCalcPEProcessor: public PEProcessor
{
public:
  TelescopePSFCalcPEProcessor(unsigned iscope = 0, bool auto_clear = true);
  virtual ~TelescopePSFCalcPEProcessor();
  void start_processing() override;
  void process_focal_plane_hit(unsigned scope_id, int pixel_id,
    double x, double y, double ux, double uy, double t, double pe_weight) override;
  void clear() { mom_.reset(); }
  const calin::math::moments_calc::SecondMomentsCalc2D mom() { return mom_; }
private:
  bool auto_clear_ = false;
  unsigned iscope_ = 0;
  calin::math::moments_calc::SecondMomentsCalc2D mom_;
};

class TelescopePSFCalcThirdMomentPEProcessor: public PEProcessor
{
public:
  TelescopePSFCalcThirdMomentPEProcessor(unsigned iscope = 0, bool auto_clear = true);
  virtual ~TelescopePSFCalcThirdMomentPEProcessor();
  void start_processing() override;
  void process_focal_plane_hit(unsigned scope_id, int pixel_id,
    double x, double y, double ux, double uy, double t, double pe_weight) override;
  void clear() { mom_.reset(); }
  const calin::math::moments_calc::ThirdMomentsCalc2D mom() { return mom_; }
private:
  bool auto_clear_ = false;
  unsigned iscope_ = 0;
  calin::math::moments_calc::ThirdMomentsCalc2D mom_;
};

} } } // namespace calin::simulation::pe_processor
