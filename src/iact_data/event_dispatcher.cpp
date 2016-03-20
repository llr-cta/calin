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
#include <io/one_to_n_data_sink.hpp>

using namespace calin::util::string;
using namespace calin::io::log;
using namespace calin::iact_data::event_dispatcher;
using namespace calin::iact_data::telescope_data_source;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;

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
  if(visitor->demand_waveforms())wf_visitors_.emplace_back(visitor);
  if(adopt_visitor)adopted_visitors_.emplace_back(visitor);
}

void TelescopeEventDispatcher::
process_run(TelescopeRandomAccessDataSourceWithRunConfig* src,
  unsigned log_frequency, int nthread)
{
  TelescopeRunConfiguration* run_config = src->get_run_configuration();
  accept_run_configuration(run_config);
  if(nthread <= 0 or std::none_of(visitors_.begin(), visitors_.end(),
    [](event_visitor::TelescopeEventVisitor* v){
      return v->is_parallelizable(); }))
  {
    accept_all_from_src(src, log_frequency, nthread==0);
  }
  else
  {
    auto visitors = visitors_;
    visitors_.clear();
    wf_visitors_.clear();
    std::vector<TelescopeEventDispatcher*> sub_dispatchers;
    for(int ithread=0;ithread<nthread;ithread++)
      sub_dispatchers.emplace_back(new TelescopeEventDispatcher);
    for(auto* v : visitors)
    {
      if(v->is_parallelizable())
        for(auto* d : sub_dispatchers)
          d->add_visitor(v->new_sub_visitor(), true);
      else
        add_visitor(v, false);
    }

    auto* sink = new io::data_source::OneToNDataSink<TelescopeEvent>;
    std::vector<std::thread> threads;
    std::atomic<unsigned> threads_active { 0 };
    // Go go gadget threads
    for(auto* d : sub_dispatchers)
    {
      d->accept_run_configuration(run_config);
      threads_active++;
      threads.emplace_back([d,sink,&threads_active](){
        auto* bsrc = sink->new_data_source();
        while(TelescopeEvent* event = bsrc->get_next())
        {
          d->accept_event(event);
          delete event;
        }
        delete bsrc;
        threads_active--;
      });
    }

    accept_all_from_src(src, log_frequency, true, sink);

    while(threads_active)
    {
      sink->put_nullptr();
      std::this_thread::sleep_for(std::chrono::microseconds(100));
    }

    for(auto& i : threads)i.join();
    delete sink;
    for(auto* d : sub_dispatchers)
    {
      d->leave_run();
      d->merge_results();
      delete d;
    }
    visitors_.clear();
    wf_visitors_.clear();
    for(auto* v : visitors)add_visitor(v, false);
  }
  leave_run();
}

void TelescopeEventDispatcher::
accept_run_configuration(TelescopeRunConfiguration* run_config)
{
  for(auto ivisitor : visitors_)ivisitor->visit_telescope_run(run_config);
}

void TelescopeEventDispatcher::accept_event(TelescopeEvent* event)
{
  for(auto ivisitor : visitors_)ivisitor->visit_telescope_event(event);

  if(!wf_visitors_.empty())
  {
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

      for(auto ivisitor : wf_visitors_)
        ivisitor->visit_waveform(ichan, hg_wf_ptr, lg_wf_ptr);
    }
  }

  for(auto ivisitor : visitors_)ivisitor->leave_telescope_event();
}

void TelescopeEventDispatcher::accept_all_from_src(
  calin::io::data_source::DataSource<
    calin::ix::iact_data::telescope_event::TelescopeEvent>* src,
  unsigned log_frequency, bool use_buffered_reader,
  calin::io::data_source::DataSink<
    calin::ix::iact_data::telescope_event::TelescopeEvent>* sink)
{
  if(!use_buffered_reader)return do_accept_from_src(src,log_frequency,sink);
  std::unique_ptr<MultiThreadTelescopeDataSourceBuffer> buffer {
    new MultiThreadTelescopeDataSourceBuffer(src) };
  std::unique_ptr<BufferedTelescopeDataSource> bsrc {
    buffer->new_data_source() };
  do_accept_from_src(bsrc.get(),log_frequency, sink);
  bsrc.reset();
}

void TelescopeEventDispatcher::
do_accept_from_src(TelescopeDataSource* src, unsigned log_frequency,
  calin::io::data_source::DataSink<
    calin::ix::iact_data::telescope_event::TelescopeEvent>* sink)
{
  using namespace std::chrono;
  uint64_t ndispatched = 0;
  auto start_time = system_clock::now();
  while(TelescopeEvent* event = src->get_next())
  {
    accept_event(event);
    if(sink)sink->put_next(event, true);
    else delete event;
    ++ndispatched;
    if(log_frequency and ndispatched % log_frequency == 0)
    {
      auto dt = system_clock::now() - start_time;
      LOG(INFO) << "Dispatched "
        << to_string_with_commas(ndispatched) << " events in "
        << to_string_with_commas(duration_cast<seconds>(dt).count()) << " sec";
    }
    //if(num_event_max and not --num_event_max)break;
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

void TelescopeEventDispatcher::leave_run()
{
  for(auto ivisitor : visitors_)ivisitor->leave_telescope_run();
}

void TelescopeEventDispatcher::merge_results()
{
  for(auto ivisitor : visitors_)ivisitor->merge_results();
}
