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

#include <simulation/vso_array.hpp>
#include <simulation/iact_array_tracker.hpp>

namespace calin { namespace simulation { namespace quadrature_iact_array_integration {

class VSO_IACTDetectorSphereHitProcessor: public IACTDetectorSphereHitProcessor
{
public:
  VSO_IACTDetectorSphereHitProcessor(
    calin::simulation::vs_optics::VSOTelescope* scope);
  virtual ~VSO_IACTDetectorSphereHitProcessor();
  virtual void process_hit(const IACTDetectorSphereHit& hit);
private:
  calin::simulation::vs_optics::VSOTelescope* scope_ = nullptr;
};

class VSO_QuadratureIACTArrayIntegrationHitVisitor:
  public calin::simulation::iact_array_tracker::HitIACTVisitor
{
public:
  VSO_QuadratureIACTArrayIntegrationHitVisitor(
    calin::simulation::vs_optics::VSOArray* array, adopt_array = false);
  ~VSO_QuadratureIACTArrayIntegrationHitVisitor();
  std::vector<calin::simulation::iact_array_tracker::IACTDetectorSphere> spheres() override;
  void visit_event(const calin::simulation::tracker::Event& event,
    bool& kill_event) override;
  void visit_cherenkov_track(
    const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack& cherenkov_track,
    bool& kill_track) override;
  void leave_cherenkov_track() override;
  void leave_event() override;
private:
  calin::simulation::vs_optics::VSOArray* array_ = nullptr;
  bool adopt_array_ = false;
};

} } } // namespace calin::simulation::quadrature_iact_array_integration
