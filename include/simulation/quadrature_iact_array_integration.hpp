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
#include <simulation/vso_raytracer.hpp>
#include <simulation/iact_array_tracker.hpp>

namespace calin { namespace simulation { namespace quadrature_iact_array_integration {

class QuadratureIACTArrayPEProcessor
{
public:
  virtual ~QuadratureIACTArrayPEProcessor();
  virtual void process_pe(unsigned scope_id, unsigned pixel_id, double t0,
    double pe_weight) = 0;
};

class SimpleImagePEProcessor: public QuadratureIACTArrayPEProcessor
{
public:
  SimpleImagePEProcessor(unsigned nscope, unsigned npix);
  SimpleImagePEProcessor(const std::vector<unsigned> npix);
  virtual ~SimpleImagePEProcessor();
  void process_pe(unsigned scope_id, unsigned pixel_id, double t0,
    double pe_weight) override;
  const std::vector<double> scope_image(unsigned iscope);
  void clear_all_images();
private:
  std::vector<std::vector<double>> images_;
};

class VSO_QuadratureIACTArrayIntegrationHitVisitor;

class VSO_IACTDetectorSphereHitProcessor:
  public calin::simulation::iact_array_tracker::IACTDetectorSphereHitProcessor
{
public:
  VSO_IACTDetectorSphereHitProcessor(
    VSO_QuadratureIACTArrayIntegrationHitVisitor* quadrature,
    calin::simulation::vs_optics::VSOTelescope* scope);
  virtual ~VSO_IACTDetectorSphereHitProcessor();
  virtual void process_hit(
    const calin::simulation::iact_array_tracker::IACTDetectorSphereHit& hit);
private:
  VSO_QuadratureIACTArrayIntegrationHitVisitor* quadrature_ = nullptr;
  calin::simulation::vs_optics::VSOTelescope* scope_ = nullptr;
};

class VSO_QuadratureIACTArrayIntegrationHitVisitor:
  public calin::simulation::iact_array_tracker::HitIACTVisitor
{
public:
  VSO_QuadratureIACTArrayIntegrationHitVisitor(double test_ray_spacing,
    calin::simulation::vs_optics::VSOArray* array,
    QuadratureIACTArrayPEProcessor* visitor,
    bool adopt_array = false, bool adopt_visitor = false);
  ~VSO_QuadratureIACTArrayIntegrationHitVisitor();
  std::vector<calin::simulation::iact_array_tracker::IACTDetectorSphere> spheres() override;
  void visit_event(const calin::simulation::tracker::Event& event,
    bool& kill_event) override;
  void visit_cherenkov_track(
    const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack& cherenkov_track,
    bool& kill_track) override;
  void leave_cherenkov_track() override;
  void leave_event() override;
protected:
  friend class VSO_IACTDetectorSphereHitProcessor;

  void process_test_ray(
    const calin::simulation::iact_array_tracker::IACTDetectorSphereHit& hit,
    calin::simulation::vs_optics::VSOTelescope* scope,
    double cos_phi, double sin_phi, double dphi);
  void process_hit(
    const calin::simulation::iact_array_tracker::IACTDetectorSphereHit& hit,
    calin::simulation::vs_optics::VSOTelescope* scope);

  double test_ray_spacing_;
  calin::simulation::vs_optics::VSOArray* array_ = nullptr;
  bool adopt_array_ = false;
  QuadratureIACTArrayPEProcessor* visitor_ = nullptr;
  bool adopt_visitor_ = false;
  calin::simulation::vs_optics::VSORayTracer* ray_tracer_ = nullptr;
};

} } } // namespace calin::simulation::quadrature_iact_array_integration
