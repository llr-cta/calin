/*

   calin/simulation/ray_processor.cpp -- Stephen Fegan -- 2017-01-16

   Multi-purpose ray (weight, scope, ray) processor.

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <io/log.hpp>
#include <simulation/ray_processor.hpp>

using namespace calin::io::log;
using namespace calin::simulation::ray_processor;

RayProcessor::~RayProcessor()
{
  // nothing to see here
}

std::vector<RayProcessorDetectorSphere> RayProcessor::detector_spheres()
{
  return {};
}

void RayProcessor::start_processing()
{
  // nothing to see here
}

void RayProcessor::process_ray(unsigned scope_id,
  const calin::math::ray::Ray& ray, double pe_weight)
{
  // nothing to see here
}

void RayProcessor::finish_processing()
{
  // nothing to see here
}

unsigned RayProcessor::process_all_from_ray_generator(
  calin::math::ray_generator::RayGenerator* gen, unsigned scope_id)
{
  unsigned n = 0;
  calin::math::ray::Ray ray;
  double weight;
  this->start_processing();
  while(gen->next(ray, weight)) {
#if 0
    if(n<10)LOG(INFO) << "[ "<< ray.position().transpose() << " ] [ "
                      << ray.direction().transpose() << " ]";
#endif
    ++n;
    this->process_ray(scope_id, ray, weight);
  }
  this->finish_processing();
  return n;
}

FederatedPEProcessor::
FederatedPEProcessor(calin::simulation::pe_processor::PEProcessor* pe_processor,
    unsigned scope_id_base):
  calin::simulation::pe_processor::PEProcessor(),
  scope_id_base_(scope_id_base), pe_processor_(pe_processor)
{
  // nothing to see here
}

FederatedPEProcessor::~FederatedPEProcessor()
{
  // nothing to see here
}

void FederatedPEProcessor::start_processing()
{
  if(scope_id_base_ == 0)pe_processor_->start_processing();
}

void FederatedPEProcessor::process_pe(unsigned scope_id, int pixel_id,
  double x, double y, double t0, double pe_weight)
{
  pe_processor_->
    process_pe(scope_id+scope_id_base_, pixel_id, x, y, t0, pe_weight);
}

void FederatedPEProcessor::finish_processing()
{
  if(scope_id_base_ == 0)pe_processor_->finish_processing();
}

FederatedRayProcessor::FederatedRayProcessor(): RayProcessor()
{

}

FederatedRayProcessor::~FederatedRayProcessor()
{
  for(auto* iprocessor : adopted_processors_)delete iprocessor;
}

std::vector<RayProcessorDetectorSphere>
FederatedRayProcessor::detector_spheres()
{
  std::vector<RayProcessorDetectorSphere> spheres;
  for(auto* iprocessor : processors_)
  {
    std::vector<RayProcessorDetectorSphere> iprocessor_spheres =
      iprocessor->detector_spheres();
    for(auto isphere : iprocessor_spheres)spheres.push_back(isphere);
  }
  return spheres;
}

void FederatedRayProcessor::start_processing()
{
  for(auto* iprocessor : processors_)iprocessor->start_processing();
}

void FederatedRayProcessor::process_ray(unsigned scope_id,
  const calin::math::ray::Ray& ray, double pe_weight)
{
  auto iprocessor =
    std::upper_bound(nsphere_.begin(), nsphere_.end(), scope_id);
  assert(iprocessor != nsphere_.end());
  if(iprocessor == nsphere_.begin())
    processors_.front()->process_ray(scope_id, ray, pe_weight);
  else
    processors_[iprocessor - nsphere_.begin()]->
      process_ray(scope_id - *(iprocessor-1), ray, pe_weight);
}

void FederatedRayProcessor::finish_processing()
{
  for(auto* iprocessor : processors_)iprocessor->finish_processing();
}

void FederatedRayProcessor::add_processor(RayProcessor* processor,
  bool adopt_processor)
{
  processors_.emplace_back(processor);
  if(adopt_processor)adopted_processors_.emplace_back(processor);
  nsphere_.emplace_back(processor->detector_spheres().size());
}

FederatedPEProcessor* FederatedRayProcessor::add_processor_and_pe_visitor(
  RayProcessor* ray_processor,
  calin::simulation::pe_processor::PEProcessor* pe_processor,
  bool adopt_ray_processor)
{
  unsigned scope_id_base = 0;
  if(nsphere_.empty())scope_id_base = nsphere_.back();
  add_processor(ray_processor, adopt_ray_processor);
  return new FederatedPEProcessor(pe_processor, scope_id_base);
}
