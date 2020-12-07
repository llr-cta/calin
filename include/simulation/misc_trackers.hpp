/*

   calin/simulation/misc_trackers.hpp -- Stephen Fegan -- 2017-07-10

   Miscellaneous class for all air shower track visitors

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

#include <Eigen/Dense>
#include <simulation/tracker.pb.h>
#include <simulation/tracker.hpp>
#include <simulation/atmosphere.hpp>
#include <simulation/air_cherenkov_tracker.hpp>
#include <simulation/detector_efficiency.hpp>

namespace calin { namespace simulation { namespace misc_trackers {

class TrackOnlyDelegateTrackVisitor: public calin::simulation::tracker::TrackVisitor
{
public:
  TrackOnlyDelegateTrackVisitor(
    calin::simulation::tracker::TrackVisitor* delegate, bool adopt_delegate = false);
  TrackOnlyDelegateTrackVisitor(
    calin::simulation::tracker::TrackVisitor* delegate, double t0, double weight,
    bool adopt_delegate = false);
  virtual ~TrackOnlyDelegateTrackVisitor();
  void visit_track(const calin::simulation::tracker::Track& track, bool& kill_track) override;
private:
  calin::simulation::tracker::TrackVisitor* delegate_ = nullptr;
  bool adopt_delegate_ = false;
  bool has_time_offset_ = false;
  double t0_ = 0;
  double weight_ = 0;
};

class LengthLimitingTrackVisitor: public calin::simulation::tracker::TrackVisitor
{
public:
  LengthLimitingTrackVisitor(calin::simulation::tracker::TrackVisitor* visitor, double dx_max,
    double z_max = std::numeric_limits<double>::infinity(),
    bool adopt_visitor = false);
  virtual ~LengthLimitingTrackVisitor();
  void visit_event(const calin::simulation::tracker::Event& event, bool& kill_event) override;
  void visit_track(const calin::simulation::tracker::Track& track, bool& kill_track) override;
  void leave_event() override;
private:
  TrackVisitor* visitor_ = nullptr;
  bool adopt_visitor_ = false;
  double dx_max_ = 0;
  double z_max_ = std::numeric_limits<double>::infinity();
};

class TimeLimitingTrackVisitor: public calin::simulation::tracker::TrackVisitor
{
public:
  TimeLimitingTrackVisitor(double t_max_sec = 60.0);
  virtual ~TimeLimitingTrackVisitor();
  void visit_event(const calin::simulation::tracker::Event& event, bool& kill_event) override;
  void visit_track(const calin::simulation::tracker::Track& track, bool& kill_track) override;
private:
  static uint64_t get_time();
  uint64_t t_event_max_ = 60000000000;
  uint64_t t_event_start_ = 0;
  bool do_kill_track_ = false;
  unsigned num_track_since_check_ = 0;
  unsigned num_track_before_check_ = 1000000;
};

class DebugStatsTrackVisitor: public calin::simulation::tracker::TrackVisitor
{
public:
  DebugStatsTrackVisitor();
  virtual ~DebugStatsTrackVisitor();
  void visit_event(const calin::simulation::tracker::Event& event, bool& kill_event) override;
  void visit_track(const calin::simulation::tracker::Track& track, bool& kill_track) override;
  std::vector<int> particle_pdg_types() const;
  unsigned track_length_zero(int pdg_type) const;
  const calin::math::histogram::SimpleHist& track_length_log_hist(int pdg_type) const;
private:
  struct TrackStats
  {
    unsigned track_length_zero = { 0 };
    calin::math::histogram::SimpleHist track_length_hist = { 0.1 };
  };
  const TrackStats& find_stats(int pdg_type) const;

  std::map<int, TrackStats> stats_;
};

class GroundInterceptTrackVisitor: public calin::simulation::tracker::TrackVisitor
{
public:
  GroundInterceptTrackVisitor(double ground_level_cm);
  virtual ~GroundInterceptTrackVisitor();
  void visit_event(const calin::simulation::tracker::Event& event, bool& kill_event) override;
  void visit_track(const calin::simulation::tracker::Track& track, bool& kill_track) override;
  void add_particle_type_filter(calin::simulation::tracker::ParticleType pt);
  std::vector<calin::simulation::tracker::Track> ground_tracks() const { return ground_tracks_; }
  void clear_tracks() { ground_tracks_.clear(); }
private:
  std::vector<calin::simulation::tracker::Track> ground_tracks_;
  double ground_level_cm_ = 0;
  std::set<calin::simulation::tracker::ParticleType> type_filter_;
};

class RecordingTrackVisitor: public calin::simulation::tracker::TrackVisitor
{
public:
  RecordingTrackVisitor();
  virtual ~RecordingTrackVisitor();
  void visit_event(const calin::simulation::tracker::Event& event, bool& kill_event) override;
  void visit_track(const calin::simulation::tracker::Track& track, bool& kill_track) override;
  void add_particle_type_filter(calin::simulation::tracker::ParticleType pt);
  calin::simulation::tracker::Event event() const { return event_; }
  const std::vector<calin::simulation::tracker::Track>& tracks() const { return tracks_; }
  void replay_event(calin::simulation::tracker::TrackVisitor* visitor) const;
  void clear_tracks() { tracks_.clear(); }
private:
  calin::simulation::tracker::Event event_;
  std::vector<calin::simulation::tracker::Track> tracks_;
  std::set<calin::simulation::tracker::ParticleType> type_filter_;
};

class SubshowerTrackVisitor: public calin::simulation::tracker::TrackVisitor
{
public:
  SubshowerTrackVisitor(double subshower_energy_mev = 1000.0);
  virtual ~SubshowerTrackVisitor();
  void visit_event(const calin::simulation::tracker::Event& event, bool& kill_event) override;
  void visit_track(const calin::simulation::tracker::Track& track, bool& kill_track) override;
  void leave_event() override;
  void clear();
  void generate_showers(calin::simulation::tracker::ShowerGenerator* generator,
    calin::simulation::tracker::TrackVisitor* visitor, unsigned num_event = 1);

  calin::simulation::tracker::Event event() const;
  const std::vector<calin::simulation::tracker::Track>& trunk_tracks() const {
    return trunk_track_visitor_.tracks();
  }
  const std::vector<calin::simulation::tracker::Event>& subshowers() const {
    return subshowers_;
  }

private:
  RecordingTrackVisitor trunk_track_visitor_;
  double subshower_energy_mev_;
  std::vector<calin::simulation::tracker::Event> subshowers_;
};

class ShowerMovieProducerTrackVisitor:
  public calin::simulation::tracker::TrackVisitor
{
public:
  CALIN_TYPEALIAS(LineSegment, std::pair<Eigen::Vector3d,Eigen::Vector3d>);

  ShowerMovieProducerTrackVisitor(calin::simulation::atmosphere::Atmosphere* atm,
    double zground,
    calin::simulation::detector_efficiency::AtmosphericAbsorption* atm_abs = nullptr,
    const calin::ix::simulation::tracker::ShowerMovieProducerTrackVisitorConfig& config = default_config(),
    bool adopt_atm = false, bool adopt_atm_abs = false);
  virtual ~ShowerMovieProducerTrackVisitor();
  void visit_event(const calin::simulation::tracker::Event& event, bool& kill_event) override;
  void visit_track(const calin::simulation::tracker::Track& track, bool& kill_track) override;
  void leave_event() override;
  static calin::ix::simulation::tracker::ShowerMovieProducerTrackVisitorConfig default_config();

  int first_frame_number() const { return frames_.empty() ? 0 : frames_.begin()->first; }
  int last_frame_number() const { return frames_.empty() ? 0 : frames_.rbegin()->first; }
  bool has_frame(int iframe) { return frames_.find(iframe) != frames_.end(); }

#define CREATE_NSEGMENT_FUNCTION(x)                  \
  unsigned nsegment_ ## x(int iframe) const {        \
    const auto frame = frames_.find(iframe);         \
    if(frame == frames_.end())return 0;              \
    return frame->second.x.size();                   \
  }

  CREATE_NSEGMENT_FUNCTION(gamma);
  CREATE_NSEGMENT_FUNCTION(electron);
  CREATE_NSEGMENT_FUNCTION(muon);
  CREATE_NSEGMENT_FUNCTION(other);
  CREATE_NSEGMENT_FUNCTION(cherenkov);

#define CREATE_SEGMENT_FUNCTION(x)                                    \
  void segment_ ## x(int iframe, unsigned isegment,                   \
      Eigen::Vector3d& x0, Eigen::Vector3d& x1) const {               \
    const auto frame = frames_.find(iframe);                          \
    if(frame == frames_.end() or isegment >= frame->second.x.size())  \
      x0 = x1 = Eigen::Vector3d::Zero();                              \
    else                                                              \
      x0 = frame->second.x[isegment].first,                           \
      x1 = frame->second.x[isegment].second;                          \
  }

  CREATE_SEGMENT_FUNCTION(gamma);
  CREATE_SEGMENT_FUNCTION(electron);
  CREATE_SEGMENT_FUNCTION(muon);
  CREATE_SEGMENT_FUNCTION(other);
  CREATE_SEGMENT_FUNCTION(cherenkov);

#define CREATE_ALL_SEGMENTS_FUNCTION(x)                               \
  void all_segments_ ## x(int iframe0, int iframe1,                   \
      Eigen::MatrixXd& x0, Eigen::MatrixXd& x1) const {               \
    unsigned nseg = 0;                                                \
    for(int iframe=iframe0;iframe<iframe1;iframe++)                   \
      nseg += nsegment_ ## x(iframe);                                 \
    x0.resize(3,nseg);                                                \
    x1.resize(3,nseg);                                                \
    nseg = 0;                                                         \
    for(int iframe=iframe0;iframe<iframe1;iframe++) {                 \
      const auto frame = frames_.find(iframe);                        \
      if(frame != frames_.end()) {                                    \
        for(const auto& seg : frame->second.x) {                      \
          x0.col(nseg) = seg.first;                                   \
          x1.col(nseg) = seg.second;                                  \
          ++nseg;                                                     \
        }                                                             \
      }                                                               \
    }                                                                 \
  }

  CREATE_ALL_SEGMENTS_FUNCTION(gamma);
  CREATE_ALL_SEGMENTS_FUNCTION(electron);
  CREATE_ALL_SEGMENTS_FUNCTION(muon);
  CREATE_ALL_SEGMENTS_FUNCTION(other);
  CREATE_ALL_SEGMENTS_FUNCTION(cherenkov);

#undef CREATE_NSEGMENT_FUNCTION
#undef CREATE_SEGMENT_FUNCTION
#undef CREATE_ALL_SEGMENTS_FUNCTION

private:
#ifndef SWIG
  class ShowerMovieProducerCherenkovPhotonVisitor:
    public calin::simulation::air_cherenkov_tracker::CherenkovPhotonVisitor
  {
  public:
    ShowerMovieProducerCherenkovPhotonVisitor(ShowerMovieProducerTrackVisitor* parent,
      calin::simulation::atmosphere::Atmosphere* atm, calin::math::rng::RNG* rng,
      calin::simulation::detector_efficiency::AtmosphericAbsorption* atm_abs,
      bool adopt_atm_abs);
    virtual ~ShowerMovieProducerCherenkovPhotonVisitor();
    void visit_event(const calin::simulation::tracker::Event& event, bool& kill_event) override;
    void visit_cherenkov_photon(const calin::simulation::air_cherenkov_tracker::CherenkovPhoton& cherenkov_photon) override;
  private:
    ShowerMovieProducerTrackVisitor* parent_;
    calin::simulation::atmosphere::Atmosphere* atm_;
    calin::math::rng::RNG* rng_;
    calin::simulation::detector_efficiency::AtmosphericAbsorption* atm_abs_;
    bool adopt_atm_abs_;
    calin::simulation::detector_efficiency::ACTEffectiveBandwidth bw_;
  };

  struct Frame {
    std::vector<LineSegment> gamma;
    std::vector<LineSegment> electron;
    std::vector<LineSegment> muon;
    std::vector<LineSegment> other;
    std::vector<LineSegment> cherenkov;
  };
  std::map<int, Frame> frames_;
  calin::ix::simulation::tracker::ShowerMovieProducerTrackVisitorConfig config_;
  double zground_;
  calin::simulation::air_cherenkov_tracker::AirCherenkovParameterCalculatorTrackVisitor* cherenkov_ = nullptr;
#endif
};

} } } // namerspace calin::simulation::misc_trackers
