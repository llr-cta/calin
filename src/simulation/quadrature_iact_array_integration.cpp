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
#include <io/log.hpp>

using namespace calin::simulation::quadrature_iact_array_integration;
using namespace calin::io::log;
using calin::math::special::SQR;
using calin::simulation::iact_array_tracker::IACTDetectorSphere;

QuadratureIACTArrayPEProcessor::~QuadratureIACTArrayPEProcessor()
{
  // nothing to see here
}

void QuadratureIACTArrayPEProcessor::visit_event(
  const calin::simulation::tracker::Event& event,
  bool& kill_event)
{
  // nothing to see here
}

void QuadratureIACTArrayPEProcessor::visit_cherenkov_track(
  const calin::simulation::air_cherenkov_tracker::AirCherenkovTrack& cherenkov_track,
  bool& kill_track)
{
  // nothing to see here
}

void QuadratureIACTArrayPEProcessor::process_pe(unsigned scope_id,
  unsigned pixel_id, double x, double y, double t0, double pe_weight)
{
  // nothing to see here
}

void QuadratureIACTArrayPEProcessor::leave_cherenkov_track()
{
  // nothing to see here
}

void QuadratureIACTArrayPEProcessor::leave_event()
{
  // nothing to see here
}

SimpleImagePEProcessor::
SimpleImagePEProcessor(unsigned nscope, unsigned npix):
  QuadratureIACTArrayPEProcessor(),
  images_(nscope, std::vector<double>(npix))
{
  // nothing to see here
}

SimpleImagePEProcessor::
SimpleImagePEProcessor(const std::vector<unsigned> npix):
  QuadratureIACTArrayPEProcessor(), images_()
{
  for(auto inpix : npix)images_.emplace_back(inpix);
}

SimpleImagePEProcessor::~SimpleImagePEProcessor()
{
  // nothing to see here
}

void SimpleImagePEProcessor::visit_event(
  const calin::simulation::tracker::Event& event, bool& kill_event)
{
  clear_all_images();
}

void SimpleImagePEProcessor::
process_pe(unsigned scope_id, unsigned pixel_id, double x, double y,
  double t0, double pe_weight)
{
  if(scope_id >= images_.size())
    throw std::out_of_range("SimpleImagePEProcessor::process_pe: scope_id out "
      "of range");
  if(pixel_id >= images_[scope_id].size())
    throw std::out_of_range("SimpleImagePEProcessor::process_pe: pixel_id out "
      "of range");
  images_[scope_id][pixel_id] += pe_weight;
}

const std::vector<double>
SimpleImagePEProcessor::scope_image(unsigned iscope) const
{
  if(iscope >= images_.size())
    throw std::out_of_range("SimpleImagePEProcessor::scope_image: iscope out "
      "of range");
  return images_[iscope];
}

void SimpleImagePEProcessor::clear_all_images()
{
  for(auto& image : images_)std::fill(image.begin(), image.end(), 0.0);
}
