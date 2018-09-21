/*

   calin/iact_data/parallel_event_dispatcher.cpp -- Stephen Fegan -- 2018-02-08

   A parallel-only dispatcher of run and event data

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

#include <chrono>
#include <type_traits>

#include <util/string.hpp>
#include <util/log.hpp>
#include <iact_data/event_dispatcher.hpp>
#include <iact_data/parallel_event_dispatcher.hpp>
#include <io/data_source.hpp>
#include <io/buffered_data_source.hpp>
#include <util/file.hpp>
#include <iact_data/nectarcam_actl_event_decoder.hpp>

using namespace calin::util::string;
using namespace calin::util::log;
using namespace calin::iact_data::event_dispatcher;
using namespace calin::io::data_source;
using namespace calin::iact_data::telescope_data_source;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using calin::iact_data::event_visitor::ParallelEventVisitor;

ParallelEventDispatcher::ParallelEventDispatcher()
{
  // nothing to see here
}

ParallelEventDispatcher::~ParallelEventDispatcher()
{
  for(auto ivisitor : adopted_visitors_)delete ivisitor;

}

void ParallelEventDispatcher::add_event_to_keep(
  calin::ix::iact_data::telescope_event::TelescopeEvent* event,
  uint64_t seq_index, google::protobuf::Arena* arena)
{
  if(event_keep_.find(event) != event_keep_.end())
    throw std::logic_error("add_event_to_keep: event already on managed event list");
  managed_event me { 0, seq_index, arena };
  event_keep_.emplace(event, me);
}

void ParallelEventDispatcher::keep_event(
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  auto ifind = event_keep_.find(event);
  if(ifind == event_keep_.end())
    throw std::logic_error("keep_event: event not on managed event list");
  ifind->second.usage_count++;
  // std::cerr << "Keep: " << event << " count " << ifind->second.usage_count << '\n';
}

void ParallelEventDispatcher::release_event(
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  auto ifind = event_keep_.find(event);
  if(ifind == event_keep_.end())
    throw std::logic_error("release_event: event not on managed event list");
  ifind->second.usage_count--;
  // std::cerr << "Release: " << event << " count " << ifind->second.usage_count << '\n';
  if(ifind->second.usage_count == 0) {
    // std::cerr << "Freeing: " << event << " arena: " << ifind->second.arena
    //   << " seq: " << ifind->second.seq_index << '\n';
    if(ifind->second.arena)delete ifind->second.arena;
    else delete event;
    event_keep_.erase(ifind);
  }
}

void ParallelEventDispatcher::
add_visitor(ParallelEventVisitor* visitor, bool adopt_visitor)
{
  visitors_.emplace_back(visitor);
  if(adopt_visitor)adopted_visitors_.emplace_back(visitor);
}

void ParallelEventDispatcher::
process_run(TelescopeRandomAccessDataSourceWithRunConfig* src,
  unsigned log_frequency, int nthread)
{
  TelescopeRunConfiguration* run_config = src->get_run_configuration();
  process_run(src, run_config, log_frequency, nthread);
  delete run_config;
}

void ParallelEventDispatcher::process_run(calin::io::data_source::DataSource<
    calin::ix::iact_data::telescope_event::TelescopeEvent>* src,
  calin::ix::iact_data::
    telescope_run_configuration::TelescopeRunConfiguration* run_config,
  unsigned log_frequency, int nthread)
{
  auto start_time = std::chrono::system_clock::now();
  std::atomic<uint_fast64_t> ndispatched { 0 };

  dispatch_run_configuration(run_config);
  if(nthread <= 0)
  {
    do_dispatcher_loop(src, log_frequency, start_time, ndispatched);
  }
  else
  {
    io::data_source::UnidirectionalBufferedDataSourcePump<TelescopeEvent> pump(src);
    do_parallel_dispatcher_loops(run_config, &pump, nthread, log_frequency,
      start_time, ndispatched);
  }
  dispatch_leave_run();
  write_final_log_message(log_frequency, start_time, ndispatched);
}

#ifdef CALIN_HAVE_CTA_CAMERASTOACTL
void ParallelEventDispatcher::process_nectarcam_zfits_run(
  const std::string& filename,
  unsigned log_frequency, int nthread,
  const calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig& decoder_config,
  const calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig& zfits_config)
{
  auto start_time = std::chrono::system_clock::now();
  std::atomic<uint_fast64_t> ndispatched { 0 };

  calin::iact_data::zfits_actl_data_source::
    ZFITSACTL_L0_CameraEventDataSource zfits_actl_src(filename, zfits_config);
  calin::iact_data::nectarcam_actl_event_decoder::
    NectarCAM_ACTL_L0_CameraEventDecoder decoder(filename,
      calin::util::file::extract_first_number_from_filename(filename),
      decoder_config);

  const DataModel::CameraEvent* actl_sample_event = nullptr;
  const DataModel::CameraRunHeader* actl_run_header = nullptr;
  try {
    zfits_actl_src.set_next_index(0);
    uint64_t unused_seq_index = 0;
    actl_sample_event = zfits_actl_src.borrow_next_event(unused_seq_index);
  } catch(...) {
    // ignore errors that occur reading sample event;
  }
  try {
    actl_run_header = zfits_actl_src.get_run_header();
  } catch(...) {
    // ignore errors that occur reading run header
  }
  calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration run_config;
  decoder.decode_run_config(&run_config, actl_run_header, actl_sample_event);
  delete actl_run_header;
  zfits_actl_src.release_borrowed_event(actl_sample_event);
  zfits_actl_src.set_next_index(0);
  run_config.clear_fragment_filename();
  for(const auto& ifilename : zfits_actl_src.source_names())
    run_config.add_fragment_filename(ifilename);

  dispatch_run_configuration(&run_config);

  zfits_actl_data_source::ZFITSConstACTL_L0_CameraEventDataSourceBorrowAdapter zfits_actl_borrow_src(&zfits_actl_src);
  zfits_actl_data_source::ZFITSConstACTL_L0_CameraEventDataSourceReleaseAdapter zfits_actl_release_sink(&zfits_actl_src);

  if(nthread <= 0)
  {
    calin::iact_data::actl_event_decoder::DecodedConstACTL_L0_CameraEventDataSource src(
      &zfits_actl_borrow_src, &zfits_actl_release_sink, &decoder);
    do_dispatcher_loop(&src, log_frequency, start_time, ndispatched);
  }
  else
  {
    calin::io::data_source::BidirectionalBufferedDataSourcePump<
      const DataModel::CameraEvent> pump(&zfits_actl_borrow_src, &zfits_actl_release_sink,
        /* buffer_size = */ 100, /* sink_unsent_data = */ true);
    calin::iact_data::actl_event_decoder::DecodedConstACTL_L0_CameraEventDataSourceFactory factory(
      &pump, &decoder);

    do_parallel_dispatcher_loops(&run_config, &factory, nthread, log_frequency,
      start_time, ndispatched);
  }

  dispatch_leave_run();
  write_final_log_message(log_frequency, start_time, ndispatched);
}
#endif // defined(CALIN_HAVE_CTA_CAMERASTOACTL)

void ParallelEventDispatcher::
dispatch_event(uint64_t seq_index, TelescopeEvent* event)
{
  for(auto iv : visitors_)iv->visit_telescope_event(seq_index, event);
}

void ParallelEventDispatcher::do_parallel_dispatcher_loops(
  calin::ix::iact_data::
    telescope_run_configuration::TelescopeRunConfiguration* run_config,
  calin::io::data_source::DataSourceFactory<
    calin::ix::iact_data::telescope_event::TelescopeEvent>* src_factory,
  unsigned nthread, unsigned log_frequency,
  const std::chrono::system_clock::time_point& start_time,
  std::atomic<uint_fast64_t>& ndispatched)
{
  std::vector<ParallelEventDispatcher*> sub_dispatchers;
  for(unsigned ithread=0;ithread<nthread;ithread++)
  {
    auto* d = new ParallelEventDispatcher;
    sub_dispatchers.emplace_back(d);

    std::map<ParallelEventVisitor*,ParallelEventVisitor*>
      antecedent_visitors;
    for(auto* v : visitors_)
    {
      ParallelEventVisitor* sv = v->new_sub_visitor(antecedent_visitors);
      if(sv) {
        d->add_visitor(sv, true);
        antecedent_visitors[v] = sv;
      } else {
        d->add_visitor(v, false);
      }
    }
    d->dispatch_run_configuration(run_config, /* dispatch_only_to_adopted_visitors = */true);
  }

  std::vector<std::thread> threads;
  std::atomic<unsigned> threads_active { 0 };

  // Go go gadget threads
  for(auto* d : sub_dispatchers)
  {
    threads_active++;
    threads.emplace_back([d,src_factory,&threads_active,log_frequency,start_time,&ndispatched](){
      try {
        auto* bsrc = src_factory->new_data_source();
        d->do_dispatcher_loop(bsrc, log_frequency, start_time, ndispatched);
        delete bsrc;
        threads_active--;
      } catch(const std::exception& x) {
        util::log::LOG(util::log::FATAL) << x.what();
        throw;
      }
    });
  }

  for(auto& i : threads)i.join();

  for(auto* d : sub_dispatchers)
  {
    d->dispatch_leave_run(/* dispatch_only_to_adopted_visitors = */ true);
    d->dispatch_merge_results(/* dispatch_only_to_adopted_visitors = */ true);
    delete d;
  }
}

void ParallelEventDispatcher::do_dispatcher_loop(
  calin::io::data_source::DataSource<
    calin::ix::iact_data::telescope_event::TelescopeEvent>* src,
  unsigned log_frequency,
  const std::chrono::system_clock::time_point& start_time,
  std::atomic<uint_fast64_t>& ndispatched)
{
  using namespace std::chrono;
  google::protobuf::Arena* arena = nullptr;
  uint64_t seq_index;
  while(TelescopeEvent* event = src->get_next(seq_index, &arena))
  {
    unsigned ndispatched_val = ndispatched.fetch_add(1);
    if(log_frequency and ndispatched_val and ndispatched_val % log_frequency == 0)
    {
      auto dt = std::chrono::system_clock::now() - start_time;
      LOG(INFO) << "Dispatched "
        << to_string_with_commas(ndispatched_val) << " events in "
        << to_string_with_commas(duration_cast<seconds>(dt).count()) << " sec";
    }

    add_event_to_keep(event, seq_index, arena);
    keep_event(event);
    dispatch_event(seq_index, event);
    release_event(event);
    arena = nullptr;
  }
}

void ParallelEventDispatcher::write_final_log_message(
  unsigned log_frequency, const std::chrono::system_clock::time_point& start_time,
  std::atomic<uint_fast64_t>& ndispatched)
{
  using namespace std::chrono;
  if(log_frequency)
  {
    auto dt = system_clock::now() - start_time;
    LOG(INFO) << "Dispatched "
      << to_string_with_commas(uint64_t(ndispatched)) << " events in "
      << to_string_with_commas(duration_cast<seconds>(dt).count()) << " sec, "
      << to_string_with_commas(duration_cast<microseconds>(dt).count()/ndispatched)
      << " us/event (finished)";
  }
}

void ParallelEventDispatcher::
dispatch_run_configuration(TelescopeRunConfiguration* run_config,
  bool dispatch_only_to_adopted_visitors)
{
  if(dispatch_only_to_adopted_visitors) {
    for(auto iv : adopted_visitors_)iv->visit_telescope_run(run_config, this);
  } else {
    for(auto iv : visitors_)iv->visit_telescope_run(run_config, this);
  }
}

void ParallelEventDispatcher::
dispatch_leave_run(bool dispatch_only_to_adopted_visitors)
{
  if(dispatch_only_to_adopted_visitors) {
    for(auto iv : adopted_visitors_)iv->leave_telescope_run();
  } else {
    for(auto iv : visitors_)iv->leave_telescope_run();
  }
}

void ParallelEventDispatcher::
dispatch_merge_results(bool dispatch_only_to_adopted_visitors)
{
  if(dispatch_only_to_adopted_visitors) {
    for(auto iv : adopted_visitors_)iv->merge_results();
  } else {
    for(auto iv : visitors_)iv->merge_results();
  }
}
