/*

   calin/simulation/quadrature_iact_array_integration.cpp
                                        -- Stephen Fegan -- 2016-07-28

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

#include <limits>

#include <simulation/iact_array_tracker.hpp>
#include <simulation/quadrature_iact_array_integration.hpp>
#include <math/special.hpp>

using namespace calin::simulation::quadrature_iact_array_integration;
using calin::math::special::SQR;
using calin::simulation::iact_array_tracker::IACTDetectorSphere;

VSO_IACTDetectorSphereHitProcessor::VSO_IACTDetectorSphereHitProcessor(
    calin::simulation::vs_optics::VSOTelescope* scope):
  IACTDetectorSphereHitProcessor(), scope_(scope)
{
  // nothing to see here
}

~VSO_IACTDetectorSphereHitProcessor()
{
  // nothing to see here
}

  virtual void process_hit(const IACTDetectorSphereHit& hit);
private:
  calin::simulation::vs_optics::VSOTelescope* scope_ = nullptr;
};


VSO_QuadratureIACTArrayIntegrationHitVisitor::
VSO_QuadratureIACTArrayIntegrationHitVisitor(
    calin::simulation::vs_optics::VSOArray* array, adopt_array):
  calin::simulation::iact_array_tracker::HitIACTVisitor(),
  array_(array), adopt_array_(adopt_array)
{
  // nothing to see here
}

VSO_QuadratureIACTArrayIntegrationHitVisitor::
~VSO_QuadratureIACTArrayIntegrationHitVisitor()
{
  if(adopt_array_)delete array_;
}

std::vector<IACTDetectorSphere>
VSO_QuadratureIACTArrayIntegrationHitVisitor::spheres()
{
  std::vector<IACTDetectorSphere> s;
  for(unsigned iscope=0; iscope<array_->numTelescopes(); iscope++)
  {
    const auto* scope = array_->telescope(iscope);
    s.emplace_back(scope->pos(), scope->curvatureRadius(), nullptr);
  }
  return s;
}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::
visit_event(const calin::simulation::tracker::Event& event, bool& kill_event)
{

}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::visit_cherenkov_track(
  const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack& cherenkov_track,
  bool& kill_track)
{

}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::leave_cherenkov_track()
{

}

void VSO_QuadratureIACTArrayIntegrationHitVisitor::leave_event()
{

}
