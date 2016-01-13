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

#include <io/nectarcam_data_source.hpp>

using namespace calin::io::nectarcam_data_source;
using namespace calin::ix::iact::telescope_event;

#ifdef CALIN_HAVE_CTA_CAMERASTOACTL

NectarCamZFITSDataSource::NectarCamZFITSDataSource(const std::string& filename):
  calin::io::telescope_data_source::TelescopeDataSource(),
  filename_(filename)
{

}

NectarCamZFITSDataSource::~NectarCamZFITSDataSource()
{
  // need to close ZFITS files
}

TelescopeEvent* NectarCamZFITSDataSource::getNextEvent()
{
  if(next_event_index_ >= zfits_->getNumMessagesInTable())return nullptr;

  auto* event =
    zfits_->readTypedMessage<DataModel::CameraEvent>(next_event_index_+1);
  assert(event);
  ++next_event_index_;

  auto* calin_event = new TelescopeEvent;
  calin_event.set_telescope_id(event->telescopeid());
  calin_event.set_local_event_number(event->eventnumber());
  calin_event.set_trigger_type(
    calin::ix::iact::telescope_event::TRIGGER_SCIENCE);
  calin_event.set_array_trigger_received(false);
  calin_event.set_array_event_number(-1);
  //c_event.local_clock_time
  calin_event.set_image_treatment(
    calin::ix::iact::telescope_event::TREATMENT_SCIENCE);

  if(event->has_higain())
  {
    auto& image = event->higain();
    auto* c_image = c_event.mutable_high_gain_image();
    if(image.has_waveforms())
    {
      auto& wf = image.waveforms();
      auto* c_wf_image = c_image->mutable_camera_waveforms();
      assert(wf.samples().type() == DataModel::AnyArray::U16);
      assert(wf.pixelsindices().type() == DataModel::AnyArray::U16);
      unsigned nsample = demand_nsample;
      if(nsample == 0)nsample = wf.num_samples();
      assert(wf.samples().data().size() % (sizeof(uint16_t)*nsample) == 0);
      unsigned npix = wf.samples().data().size()/(sizeof(uint16_t)*nsample);
      const uint16_t* wf_data =
        reinterpret_cast<const uint16_t*>(&wf.samples().data().front());
      for(unsigned ipix=0;ipix<npix;ipix++)
      {
        auto* c_wf = c_wf_image->add_waveform();
        for(unsigned isample=0;isample<nsample;isample++)
          c_wf->add_samples(*wf_data++);
      }
    }
  }

  if(event->has_logain())
  {

  }

}

#endif // #ifdef CALIN_HAVE_CTA_CAMERASTOACTL
