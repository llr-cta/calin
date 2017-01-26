/*

   calin/simulation/pe_processor.hpp -- Stephen Fegan -- 2017-01-16

   Multi-purpose PE (weight, scope, pixel & time) processor.

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

#pragma once

#include <vector>

#include <math/moments_calc.hpp>

namespace calin { namespace simulation { namespace pe_processor {

class PEProcessor
{
public:
  virtual ~PEProcessor();
  virtual void start_processing();
  virtual void process_pe(unsigned scope_id, int pixel_id,
    double x, double y, double t0, double pe_weight);
  virtual void finish_processing();
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
  void process_pe(unsigned scope_id, int pixel_id,
    double x, double y, double t0, double pe_weight) override;
  const std::vector<double> scope_image(unsigned iscope) const;
  void clear_all_images();
private:
  bool auto_clear_ = false;
  std::vector<std::vector<double>> images_;
};

class TelescopePSFCalcPEProcessor: public PEProcessor
{
public:
  TelescopePSFCalcPEProcessor(unsigned iscope = 0, bool auto_clear = true);
  virtual ~TelescopePSFCalcPEProcessor();
  void start_processing() override;
  void process_pe(unsigned scope_id, int pixel_id,
    double x, double y, double t0, double pe_weight) override;
  void clear() { mom_.reset(); }
  const calin::math::moments_calc::SecondMomentsCalc2D mom() { return mom_; }
private:
  bool auto_clear_ = false;
  unsigned iscope_ = 0;
  calin::math::moments_calc::SecondMomentsCalc2D mom_;
};

} } } // namespace calin::simulation::pe_processor
