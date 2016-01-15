/*

   calin/io/nectarcam_data_source.cpp -- Stephen Fegan -- 2016-01-11

   A supplier of single telescope data from NectarCam DAQ data files

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

#include <stdexcept>

#include <io/nectarcam_data_source.hpp>

using namespace calin::io::nectarcam_data_source;
using namespace calin::ix::iact::telescope_event;

#ifdef CALIN_HAVE_CTA_CAMERASTOACTL

#include <ProtobufIFits.h>
#include <L0.pb.h>

NectarCamZFITSDataSource::
NectarCamZFITSDataSource(const std::string& filename):
  calin::io::telescope_data_source::TelescopeDataSource(),
  filename_(filename)
{
  zfits_ = new ACTL::IO::ProtobufIFits(filename_.c_str());
}

NectarCamZFITSDataSource::~NectarCamZFITSDataSource()
{
  delete zfits_;
}

TelescopeEvent* NectarCamZFITSDataSource::getNextEvent()
{
  if(zfits_ == nullptr)return nullptr;
  if(next_event_index_ >= zfits_->getNumMessagesInTable() or
    (config_.num_events_max() and next_event_index_>=config_.num_events_max()))
  {
    delete zfits_;
    zfits_ = nullptr;
    return nullptr;
  }

  const auto* cta_event =
    zfits_->readTypedMessage<DataModel::CameraEvent>(next_event_index_+1);
  assert(cta_event);
  ++next_event_index_;

  auto* calin_event = new TelescopeEvent;
  calin_event->set_telescope_id(cta_event->telescopeid());
  calin_event->set_local_event_number(cta_event->eventnumber());
  calin_event->set_trigger_type(TRIGGER_SCIENCE);
  calin_event->set_array_trigger_received(false);
  calin_event->set_array_event_number(-1);
  //calin_event->local_clock_time
  calin_event->set_image_treatment(TREATMENT_SCIENCE);

  if(cta_event->has_higain())
    copy_single_gain_image(cta_event->higain(),
      calin_event->mutable_high_gain_image());

  if(cta_event->has_logain())
    copy_single_gain_image(cta_event->logain(),
      calin_event->mutable_low_gain_image());

  if(cta_event->has_cameracounters())
  {
    struct NectarCounters {
      uint32_t global_event_counter;
      uint16_t bunch_counter;
      uint16_t event_counter;
      uint32_t ts1;
      uint8_t  ts2_bunch;
      uint8_t  ts2_event;
      uint16_t ts2_empty;
    }__attribute__((packed));

    auto& cta_counters = cta_event->cameracounters();
    assert(cta_counters.has_counters());
    assert(cta_counters.counters().type() == DataModel::AnyArray::U16);
    assert(cta_counters.counters().data().size()%sizeof(NectarCounters)
      == 0);
    unsigned nmod =
      cta_counters.counters().data().size()/sizeof(NectarCounters);
    const auto* mod_counter = reinterpret_cast<const NectarCounters*>(
      &cta_counters.counters().data().front());
    for(unsigned imod=0;imod<nmod;imod++)
    {
#define ts2_decode(x) ((x)&0xF0?((x)&0xC0?((x)&0x80?0:1):((x)&0x20?2:3)):\
                                ((x)&0x0C?((x)&0x08?4:5):((x)&0x02?6:7)))
      int ts2_bunch = ts2_decode(mod_counter->ts2_bunch);
      int ts2_event = ts2_decode(mod_counter->ts2_event);
      int64_t time_ns = mod_counter->bunch_counter*1000000000LL
        + mod_counter->ts1*8LL + ts2_event - ts2_bunch;
      auto* clock = calin_event->add_local_clock_time();
      clock->set_clock_id(imod);
      clock->mutable_clock_time()->set_time_ns(time_ns);
      mod_counter++;
    }
  }

  delete cta_event;
  return calin_event;
}

void NectarCamZFITSDataSource::
copy_single_gain_image(const DataModel::PixelsChannel& cta_image,
  calin::ix::iact::telescope_event::DigitizedSkyImage* calin_image)
{
  if(cta_image.has_waveforms())
  {
    const auto& cta_wf = cta_image.waveforms();
    auto* calin_wf_image = calin_image->mutable_camera_waveforms();
    assert(cta_wf.has_samples());
    assert(cta_wf.samples().type() == DataModel::AnyArray::U16);
    unsigned nsample = config_.demand_nsample();
    if(nsample == 0)nsample = cta_wf.num_samples();
    assert(cta_wf.samples().data().size() % (sizeof(uint16_t)*nsample) == 0);
    unsigned npix = cta_wf.samples().data().size()/(sizeof(uint16_t)*nsample);
    const uint16_t* cta_wf_data =
      reinterpret_cast<const uint16_t*>(&cta_wf.samples().data().front());
    for(unsigned ipix=0;ipix<npix;ipix++)
    {
      auto* calin_wf = calin_wf_image->add_waveform();
      for(unsigned isample=0;isample<nsample;isample++)
        calin_wf->add_samples(*cta_wf_data++);
    }
  }

  if(cta_image.has_integrals())
  {
    const auto& cta_q = cta_image.integrals();
    auto* calin_q_image = calin_image->mutable_camera_charges();
    assert(cta_q.has_gains());
    assert(cta_q.gains().type() == DataModel::AnyArray::U16);
    assert(cta_q.gains().data().size() % sizeof(uint16_t) == 0);
    unsigned npix = cta_q.gains().data().size()/sizeof(uint16_t);
    const uint16_t* cta_q_data =
      reinterpret_cast<const uint16_t*>(&cta_q.gains().data().front());
    for(unsigned ipix=0;ipix<npix;ipix++)
      calin_q_image->add_charge(*cta_q_data++);

  }
}

#else // #ifdef CALIN_HAVE_CTA_CAMERASTOACTL

NectarCamZFITSDataSource::
NectarCamZFITSDataSource(const std::string& filename):
  calin::io::telescope_data_source::TelescopeDataSource()
{
  // nothing to see here
}

NectarCamZFITSDataSource::~NectarCamZFITSDataSource()
{
  // nothing to see here
}

TelescopeEvent* NectarCamZFITSDataSource::getNextEvent()
{
  throw(std::logic_error("NectarCamZFITSDataSource::getNextEvent(): calin "
    "compiled without ZFits support."))
  return nullptr;
}

#endif // #ifdef CALIN_HAVE_CTA_CAMERASTOACTL
