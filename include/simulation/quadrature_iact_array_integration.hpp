/*

   calin/simulation/quadrature_iact_array_integration.hpp
                                              -- Stephen Fegan -- 2016-07-27

   IACT array tracker that does quadrature integration over the cherenkov
   cone.

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

#pragma once

#include <simulation/iact_array_tracker.hpp>
#include <simulation/tracker.pb.h>

namespace calin { namespace simulation { namespace quadrature_iact_array_integration {

class QuadratureIACTArrayPEProcessor
{
public:
  virtual ~QuadratureIACTArrayPEProcessor();
  virtual void visit_event(const calin::simulation::tracker::Event& event,
    bool& kill_event);
  virtual void visit_cherenkov_track(
    const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack& cherenkov_track,
    bool& kill_track);
  virtual void process_pe(unsigned scope_id, unsigned pixel_id,
    double x, double y, double t0, double pe_weight);
  virtual void leave_cherenkov_track();
  virtual void leave_event();
};

class SimpleImagePEProcessor: public QuadratureIACTArrayPEProcessor
{
public:
  SimpleImagePEProcessor(unsigned nscope, unsigned npix);
  SimpleImagePEProcessor(const std::vector<unsigned> npix);
  virtual ~SimpleImagePEProcessor();
  void visit_event(const calin::simulation::tracker::Event& event,
    bool& kill_event) override;
  void process_pe(unsigned scope_id, unsigned pixel_id,
    double x, double y, double t0, double pe_weight) override;
  const std::vector<double> scope_image(unsigned iscope) const;
  void clear_all_images();
private:
  std::vector<std::vector<double>> images_;
};

} } } // namespace calin::simulation::quadrature_iact_array_integration
