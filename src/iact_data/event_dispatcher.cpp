/*

   calin/iact_data/event_dispatcher.cpp -- Stephen Fegan -- 2016-02-09

   A dispatcher of run and event data

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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
#include <io/one_to_n_data_sink.hpp>
#include <util/file.hpp>

using namespace calin::util::string;
using namespace calin::util::log;
using namespace calin::iact_data::event_dispatcher;
using namespace calin::iact_data::telescope_data_source;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using calin::iact_data::event_visitor::TelescopeEventVisitor;

TelescopeEventDispatcher::TelescopeEventDispatcher()
{
  // nothing to see here
}

TelescopeEventDispatcher::~TelescopeEventDispatcher()
{
  for(auto ivisitor : adopted_visitors_)delete ivisitor;
}

void TelescopeEventDispatcher::
add_visitor(TelescopeEventVisitor* visitor,
  VisitorExecutionMode execution_mode, bool adopt_visitor)
{
  if((execution_mode == EXECUTE_SEQUENTIAL_AND_PARALLEL or
      execution_mode == EXECUTE_PARALLEL) and !visitor->is_parallelizable())
    throw std::runtime_error(
      "add_visitor: parallel execution mode requested on incompatible visitor");

  visitors_.emplace_back(execution_mode, visitor);
  if(visitor->demand_waveforms())wf_visitors_.emplace_back(visitor);
  if(adopt_visitor)adopted_visitors_.emplace_back(visitor);
}

void TelescopeEventDispatcher::
process_run(TelescopeRandomAccessDataSourceWithRunConfig* src,
  unsigned log_frequency, int nthread)
{
  TelescopeRunConfiguration* run_config = src->get_run_configuration();
  process_run(src, run_config, log_frequency, nthread);
  delete run_config;
}

void TelescopeEventDispatcher::process_run(calin::io::data_source::DataSource<
    calin::ix::iact_data::telescope_event::TelescopeEvent>* src,
  calin::ix::iact_data::
    telescope_run_configuration::TelescopeRunConfiguration* run_config,
  unsigned log_frequency, int nthread, bool use_buffered_reader)
{
  accept_run_configuration(run_config);
  if(nthread <= 0 or std::none_of(visitors_.begin(), visitors_.end(),
    [](const std::pair<VisitorExecutionMode,TelescopeEventVisitor*>& iv)
      { return (iv.first == EXECUTE_PARALLEL_IF_POSSIBLE or
          iv.first == EXECUTE_SEQUENTIAL_AND_PARALLEL or
          iv.first == EXECUTE_SEQUENTIAL_AND_PARALLEL_IF_POSSIBLE or
          iv.first == EXECUTE_PARALLEL) and iv.second->is_parallelizable(); }))
  {
    accept_all_from_src(src, log_frequency, nthread==0 and use_buffered_reader);
  }
  else
  {
    auto visitors = visitors_;
    visitors_.clear();
    wf_visitors_.clear();
    std::vector<TelescopeEventDispatcher*> sub_dispatchers;
    for(int ithread=0;ithread<nthread;ithread++)
      sub_dispatchers.emplace_back(new TelescopeEventDispatcher);

    std::vector<TelescopeEventVisitor*> parallelizable_visitors;
    for(auto iv : visitors)
    {
      bool add_sequential = false;
      bool add_parallel = false;
      auto* v = iv.second;
      switch(iv.first)
      {
      case EXECUTE_PARALLEL_IF_POSSIBLE:
        if(v->is_parallelizable())add_parallel=true;
        else add_sequential=true;
        break;
      case EXECUTE_SEQUENTIAL:
        add_sequential=true;
        break;
      case EXECUTE_SEQUENTIAL_AND_PARALLEL:
        add_sequential = true;
        add_parallel = true;
        break;
      case EXECUTE_SEQUENTIAL_AND_PARALLEL_IF_POSSIBLE:
        if(v->is_parallelizable())add_parallel=true;
        add_sequential=true;
        break;
      case EXECUTE_PARALLEL:
        add_parallel = true;
        break;
      }
      if(add_parallel)parallelizable_visitors.emplace_back(v);
      if(add_sequential)add_visitor(v, EXECUTE_SEQUENTIAL, false);
    }

    for(auto* d : sub_dispatchers)
    {
      std::map<TelescopeEventVisitor*,TelescopeEventVisitor*>
        antecedent_visitors;
      for(auto* v : parallelizable_visitors)
      {
        TelescopeEventVisitor* sv = v->new_sub_visitor(antecedent_visitors);
        d->add_visitor(sv, EXECUTE_SEQUENTIAL, true);
        antecedent_visitors[v] = sv;
      }
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
        google::protobuf::Arena* arena = nullptr;
        uint64_t seq_index;
        while(TelescopeEvent* event = bsrc->get_next(seq_index, &arena))
        {
          d->accept_event(seq_index, event);
          if(arena)delete arena;
          else delete event;
          arena = nullptr;
        }
        delete bsrc;
        threads_active--;
      });
    }

    accept_all_from_src(src, log_frequency, use_buffered_reader, sink);

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
    for(auto iv : visitors)add_visitor(iv.second, iv.first, false);
  }
  leave_run();
}

#ifdef CALIN_HAVE_CTA_CAMERASTOACTL
void TelescopeEventDispatcher::process_nectarcam_zfits_run(
  const std::string& filename,
  unsigned log_frequency, int nthread,
  const calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig& decoder_config,
  const calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig& zfits_config)
{
  calin::iact_data::zfits_actl_data_source::
    ZFITSACTLDataSource zfits_actl_src(filename, zfits_config);
  calin::iact_data::nectarcam_data_source::
    NectarCamCameraEventDecoder decoder(filename,
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

  zfits_actl_data_source::ZFITSConstACTLDataSourceBorrowAdapter zfits_actl_borrow_src(&zfits_actl_src);
  zfits_actl_data_source::ZFITSConstACTLDataSourceReleaseAdapter zfits_actl_release_sink(&zfits_actl_src);

  calin::io::data_source::BidirectionalBufferedDataSourcePump<
    const DataModel::CameraEvent>* pump_actl_src = nullptr;

  calin::iact_data::zfits_data_source::DecodedConstACTLDataSource* src;

  if(nthread == -1) {
    src = new calin::iact_data::zfits_data_source::DecodedConstACTLDataSource(
      &zfits_actl_borrow_src, &zfits_actl_release_sink, &decoder);
  } else {
    pump_actl_src = new calin::io::data_source::BidirectionalBufferedDataSourcePump<
      const DataModel::CameraEvent>(&zfits_actl_borrow_src, &zfits_actl_release_sink,
        /* buffer_size = */ 100, /* sink_unsent_data = */ true);
    src = new calin::iact_data::zfits_data_source::DecodedConstACTLDataSource(
      pump_actl_src->new_data_source(), pump_actl_src->new_data_sink(),
      &decoder, /* adopt_actl_src = */ true, /* adopt_actl_sink = */ true);
  }

  process_run(src, &run_config, log_frequency, nthread, false);

  delete src;
  delete pump_actl_src;
}
#endif // defined(CALIN_HAVE_CTA_CAMERASTOACTL)

void TelescopeEventDispatcher::
accept_run_configuration(TelescopeRunConfiguration* run_config)
{
  for(auto iv : visitors_)iv.second->visit_telescope_run(run_config);
}

void TelescopeEventDispatcher::
accept_event(uint64_t seq_index, TelescopeEvent* event)
{
  for(auto iv : visitors_)iv.second->visit_telescope_event(seq_index, event);

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

  for(auto iv : visitors_)iv.second->leave_telescope_event();
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
  google::protobuf::Arena* arena = nullptr;
  uint64_t seq_index;
  while(TelescopeEvent* event = src->get_next(seq_index, &arena))
  {
    if(log_frequency and ndispatched and ndispatched % log_frequency == 0)
    {
      auto dt = system_clock::now() - start_time;
      LOG(INFO) << "Dispatched "
        << to_string_with_commas(ndispatched) << " events in "
        << to_string_with_commas(duration_cast<seconds>(dt).count()) << " sec";
    }

    accept_event(seq_index, event);
    if(sink)sink->put_next(event, seq_index, arena, true);
    else if(arena)delete arena;
    else delete event;
    ++ndispatched;
    arena = nullptr;
    //if(num_event_max and not --num_event_max)break;
  }
  if(log_frequency)
  {
    auto dt = system_clock::now() - start_time;
    LOG(INFO) << "Dispatched "
      << to_string_with_commas(ndispatched) << " events in "
      << to_string_with_commas(duration_cast<seconds>(dt).count()) << " sec, "
      << to_string_with_commas(duration_cast<microseconds>(dt).count()/ndispatched)
      << " us/event (finished)";
  }
}

void TelescopeEventDispatcher::leave_run()
{
  for(auto iv : visitors_)iv.second->leave_telescope_run();
}

void TelescopeEventDispatcher::merge_results()
{
  for(auto iv : visitors_)iv.second->merge_results();
}
