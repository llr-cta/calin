/*

   calin/proto/diagnostics/run_coherence_mif.cpp -- Stephen Fegan -- 2018-10-26

   Protobufs for run coherence diagnostics

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <algorithm>

namespace {
  struct range_and_value {
    uint64_t begin_index;
    uint64_t end_index;
    int64_t value;
  };

  bool operator<(const range_and_value& a, const range_and_value& b) {
    return a.begin_index<b.begin_index; }
}

#include "diagnostics/run_coherence.pb.h"

void calin::ix::diagnostics::run_coherence::IndexRange::IntegrateFrom(
  const calin::ix::diagnostics::run_coherence::IndexRange& from)
{
  if(this->begin_index_size() != this->end_index_size())
    throw std::runtime_error("IntegrateFrom: index and value arrays must be same size");

  if(from.begin_index_size() != from.end_index_size())
    throw std::runtime_error("IntegrateFrom: index and value arrays must be same size");

  std::vector<range_and_value> values;
  values.reserve(this->begin_index_size() + from.begin_index_size());
  for(unsigned i=0; i<this->begin_index_size(); i++)
    values.push_back({this->begin_index(i),this->end_index(i),0});
  for(unsigned i=0; i<from.begin_index_size(); i++)
    values.push_back({from.begin_index(i),from.end_index(i),0});

  if(values.size() == 0)return;
  std::sort(values.begin(), values.end());

  this->clear_begin_index();
  this->clear_end_index();

  range_and_value v = values.front();
  for(unsigned i=1; i<values.size(); i++)
  {
    if(values[i].begin_index < v.end_index)
      throw std::runtime_error("Cannot integrate IntegrateFrom, ranges overlap: "
        + std::to_string(v.end_index) + " > " + std::to_string(values[i].begin_index));

    if(values[i].begin_index == v.end_index and values[i].value == v.value) {
      v.end_index = values[i].end_index;
    } else {
      this->add_begin_index(v.begin_index);
      this->add_end_index(v.end_index);
      v = values[i];
    }
  }
  this->add_begin_index(v.begin_index);
  this->add_end_index(v.end_index);
}

void calin::ix::diagnostics::run_coherence::CoherentValueRange::IntegrateFrom(
  const calin::ix::diagnostics::run_coherence::CoherentValueRange& from)
{
  if(this->value_size() != this->begin_index_size() or
      this->value_size() != this->end_index_size())
    throw std::runtime_error("CoherentValueRange: index and value arrays must be same size");

  if(from.value_size() != from.begin_index_size() or
      from.value_size() != from.end_index_size())
    throw std::runtime_error("CoherentValueRange: index and value arrays must be same size");

    std::vector<range_and_value> values;
    values.reserve(this->begin_index_size() + from.begin_index_size());
    for(unsigned i=0; i<this->begin_index_size(); i++)
      values.push_back({this->begin_index(i),this->end_index(i),this->value(i)});
    for(unsigned i=0; i<from.begin_index_size(); i++)
      values.push_back({from.begin_index(i),from.end_index(i),from.value(i)});

  if(values.size() == 0)return;
  std::sort(values.begin(), values.end());

  this->clear_value();
  this->clear_begin_index();
  this->clear_end_index();

  range_and_value v = values.front();
  for(unsigned i=1; i<values.size(); i++)
  {
    if(values[i].begin_index < v.end_index)
      throw std::runtime_error("Cannot integrate CoherentValueRange, ranges overlap: "
        + std::to_string(v.end_index) + " > " + std::to_string(values[i].begin_index));

    if(values[i].begin_index == v.end_index and values[i].value == v.value) {
      v.end_index = values[i].end_index;
    } else {
      this->add_begin_index(v.begin_index);
      this->add_end_index(v.end_index);
      this->add_value(v.value);
      v = values[i];
    }
  }
  this->add_begin_index(v.begin_index);
  this->add_end_index(v.end_index);
  this->add_value(v.value);
}
