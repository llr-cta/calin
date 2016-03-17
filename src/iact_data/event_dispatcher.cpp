/*

   calin/iact_data/event_dispatcher.cpp -- Stephen Fegan -- 2016-02-09

   A dispatcher of run and event data

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

#include <chrono>

#include <util/string.hpp>
#include <io/log.hpp>
#include <iact_data/event_dispatcher.hpp>

using namespace calin::util::string;
using namespace calin::io::log;
using namespace calin::iact_data::event_dispatcher;
using namespace calin::iact_data::telescope_data_source;
using namespace calin::ix::iact_data::telescope_event;

TelescopeEventDispatcher::TelescopeEventDispatcher()
{
  // nothing to see here
}

TelescopeEventDispatcher::~TelescopeEventDispatcher()
{
  for(auto ivisitor : adopted_visitors_)delete ivisitor;
}

void TelescopeEventDispatcher::
add_visitor(event_visitor::TelescopeEventVisitor* visitor, bool adopt_visitor)
{
  visitors_.emplace_back(visitor);
  if(adopt_visitor)adopted_visitors_.emplace_back(visitor);
}

void TelescopeEventDispatcher::accept(TelescopeEvent* event)
{
  for(auto ivisitor : visitors_)ivisitor->visit_telescope_event(event);

  Waveforms* hg_wf = nullptr;
  if(event->has_high_gain_image() and
    event->high_gain_image().has_camera_waveforms())
  {
    hg_wf = event->mutable_high_gain_image()->mutable_camera_waveforms();
    assert(hg_wf->channel_id_size() == 0 or
      hg_wf->channel_id_size() == hg_wf->waveform_size());
  }

  Waveforms* lg_wf = nullptr;
  if(event->has_low_gain_image() and
    event->low_gain_image().has_camera_waveforms())
  {
    lg_wf = event->mutable_low_gain_image()->mutable_camera_waveforms();
    assert(lg_wf->channel_id_size() == 0 or
      lg_wf->channel_id_size() == lg_wf->waveform_size());
  }

  int hg_iwf = 0;
  int lg_iwf = 0;

  while((!hg_wf or hg_iwf<hg_wf->waveform_size()) or
        (!lg_wf or lg_iwf<lg_wf->waveform_size()))
  {
    int hg_ichan = hg_iwf;
    if(hg_wf and hg_iwf<hg_wf->channel_id_size())
      hg_ichan = hg_wf->channel_id(hg_iwf);
    int lg_ichan = lg_iwf;
    if(lg_wf and lg_iwf<lg_wf->channel_id_size())
      lg_ichan = lg_wf->channel_id(lg_iwf);
    int ichan = std::min(hg_ichan, lg_ichan);
    ChannelWaveform* hg_wf_ptr = nullptr;
    if(hg_wf and hg_iwf<hg_wf->waveform_size() and ichan==hg_ichan)
      hg_wf_ptr = hg_wf->mutable_waveform(hg_iwf++);
    ChannelWaveform* lg_wf_ptr = nullptr;
    if(lg_wf and lg_iwf<lg_wf->waveform_size() and ichan==lg_ichan)
      lg_wf_ptr = lg_wf->mutable_waveform(lg_iwf++);
    for(auto ivisitor : visitors_)
      ivisitor->visit_waveform(ichan, hg_wf_ptr, lg_wf_ptr);
  }

  for(auto ivisitor : visitors_)ivisitor->leave_telescope_event();
}

void TelescopeEventDispatcher::accept_all_from_src(
  calin::io::data_source::DataSource<
    calin::ix::iact_data::telescope_event::TelescopeEvent>* src,
  unsigned log_frequency, bool use_buffered_reader)
{
  if(!use_buffered_reader)return accept_from_src(src,log_frequency);
  std::unique_ptr<MultiThreadTelescopeDataSourceBuffer> buffer {
    new MultiThreadTelescopeDataSourceBuffer(src) };
  std::unique_ptr<BufferedTelescopeDataSource> bsrc {
    buffer->new_data_source(100) };
  accept_from_src(bsrc.get(),log_frequency);
  bsrc.reset();
}

void TelescopeEventDispatcher::
accept_from_src(TelescopeDataSource* src, unsigned log_frequency,
  unsigned num_event_max)
{
  using namespace std::chrono;
  uint64_t ndispatched = 0;
  auto start_time = system_clock::now();
  while(TelescopeEvent* event = src->get_next())
  {
    accept(event);
    delete event;
    ++ndispatched;
    if(log_frequency and ndispatched % log_frequency == 0)
    {
      auto dt = system_clock::now() - start_time;
      LOG(INFO) << "Dispatched "
        << to_string_with_commas(ndispatched) << " events in "
        << to_string_with_commas(duration_cast<seconds>(dt).count()) << " sec";
    }
    if(num_event_max and not --num_event_max)break;
  }
  if(log_frequency and ndispatched % log_frequency != 0)
  {
    auto dt = system_clock::now() - start_time;
    LOG(INFO) << "Dispatched "
      << to_string_with_commas(ndispatched) << " events in "
      << to_string_with_commas(duration_cast<seconds>(dt).count())
      << " sec (finished)";
  }
}
