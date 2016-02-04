/*

   calin/iact_data/nectarcam_data_source.cpp -- Stephen Fegan -- 2016-01-11

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
#include <string>

#include <io/log.hpp>
#include <util/file.hpp>
#include <iact_data/nectarcam_data_source.hpp>
#include <iact_data/zfits_data_source.hpp>

using namespace calin::iact_data::nectarcam_data_source;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::io::log;
using calin::iact_data::zfits_data_source::ZFITSSingleFileDataSource;

#include <ProtobufIFits.h>
#include <L0.pb.h>

#define TEST_ANYARRAY_TYPES 0

NectarCamCameraEventDecoder::~NectarCamCameraEventDecoder()
{
  // nothing to see here
}

calin::ix::iact_data::telescope_event::TelescopeEvent*
NectarCamCameraEventDecoder::decode(const DataModel::CameraEvent* cta_event)
{
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
      calin_event->mutable_high_gain_image(), "high");

  if(cta_event->has_logain())
    copy_single_gain_image(cta_event->logain(),
      calin_event->mutable_low_gain_image(), "low");

  if(cta_event->has_cameracounters() and
    cta_event->cameracounters().has_counters())
  {
    struct NectarCounters {
      uint32_t global_event_counter;
      uint16_t bunch_counter;
      uint16_t event_counter;
      uint32_t ts1;
      uint16_t ts2_empty;
      uint8_t  ts2_bunch;
      uint8_t  ts2_event;
    }__attribute__((packed));

    const auto& cta_counters = cta_event->cameracounters().counters();
#if TEST_ANYARRAY_TYPES
    if(cta_counters.type() != DataModel::AnyArray::U16)
      throw std::runtime_error("Camera counters type not U16");
#endif
    if(cta_counters.data().size()%sizeof(NectarCounters) != 0)
      throw std::runtime_error("Camera counters data array not integral "
        "multiple of expected structure size.");
    unsigned nmod =
      cta_counters.data().size()/sizeof(NectarCounters);
    const auto* mod_counter =
      reinterpret_cast<const NectarCounters*>(&cta_counters.data().front());
    for(unsigned imod=0;imod<nmod;imod++)
    {
      auto* module_counters = calin_event->add_module_counter();
      module_counters->set_module_id(imod);
#define add_counter(id,val) \
      { \
        auto* counter = module_counters->add_counter(); \
        counter->set_counter_id(id); \
        counter->set_value(val); \
      }
      add_counter(0, mod_counter->global_event_counter);
      add_counter(1, mod_counter->bunch_counter);
      add_counter(2, mod_counter->event_counter);
      add_counter(3, mod_counter->ts1);
      add_counter(4, mod_counter->ts2_bunch);
      add_counter(5, mod_counter->ts2_event);

#define ts2_decode(x) ((x)&0xF0?((x)&0xC0?((x)&0x80?0:1):((x)&0x20?2:3)):\
                                ((x)&0x0C?((x)&0x08?4:5):((x)&0x02?6:7)))
      int ts2_bunch = ts2_decode(mod_counter->ts2_bunch);
      int ts2_event = ts2_decode(mod_counter->ts2_event);
      int64_t time_ns = mod_counter->bunch_counter*1000000000LL
        + mod_counter->ts1*8LL + ts2_event - ts2_bunch;
      auto* module_clocks = calin_event->add_module_clock();
      module_clocks->set_module_id(imod);
      auto* clock = module_clocks->add_clock();
      clock->set_clock_id(0);
      clock->mutable_time()->set_time_ns(time_ns);
      mod_counter++;
    }
  }

  return calin_event;
}

void NectarCamCameraEventDecoder::
copy_single_gain_image(const DataModel::PixelsChannel& cta_image,
  calin::ix::iact_data::telescope_event::DigitizedSkyImage* calin_image,
  const std::string& which_gain)
{
  if(cta_image.has_waveforms() and cta_image.waveforms().has_samples())
  {
    const auto& cta_wf = cta_image.waveforms().samples();
    auto* calin_wf_image = calin_image->mutable_camera_waveforms();
#if TEST_ANYARRAY_TYPES
    if(cta_wf.type() != DataModel::AnyArray::U16)
      throw std::runtime_error("Waveform data type not uint16 in " +
        which_gain + " channel.");
#endif
    unsigned nsample = config_.demand_nsample();
    if(nsample == 0)nsample = cta_image.waveforms().num_samples();
    if(nsample == 0)throw std::runtime_error("Number of samples is zero in "
      + which_gain + " channel.");
    if(cta_wf.data().size() % (sizeof(uint16_t)*nsample) != 0)
      throw std::runtime_error("Waveform data array for " + which_gain +
        " gain channel not integral multiple of nsample uint16.");
    calin_wf_image->set_num_samples_per_channel(nsample);
    unsigned npix = cta_wf.data().size()/(sizeof(uint16_t)*nsample);
    const uint16_t* cta_wf_data =
      reinterpret_cast<const uint16_t*>(&cta_wf.data().front());
    for(unsigned ipix=0;ipix<npix;ipix++)
    {
      auto* calin_wf = calin_wf_image->add_waveform();
      for(unsigned isample=0;isample<nsample;isample++)
        calin_wf->add_samples(*cta_wf_data++);
    }
  }

  if(cta_image.has_integrals() and cta_image.integrals().has_gains())
  {
    const auto& cta_q = cta_image.integrals().gains();
    auto* calin_q_image = calin_image->mutable_camera_charges();
#if TEST_ANYARRAY_TYPES
    if(cta_q.type() != DataModel::AnyArray::U16)
      throw std::runtime_error("Integral data type not uint16 in " +
        which_gain + " channel.");
#endif
    if(cta_q.data().size() % sizeof(uint16_t) != 0)
      throw std::runtime_error("Charge data array for " + which_gain +
        " gain channel not integral multiple of uint16.");
    unsigned npix = cta_q.data().size()/sizeof(uint16_t);
    const uint16_t* cta_q_data =
      reinterpret_cast<const uint16_t*>(&cta_q.data().front());
    for(unsigned ipix=0;ipix<npix;ipix++)
      calin_q_image->add_charge(*cta_q_data++);
  }
}

NectarCamZFITSDataSource::
NectarCamZFITSDataSource(const std::string& filename,
  const decoder_config_type& decoder_config, const config_type& config):
  calin::iact_data::zfits_data_source::ZFITSDataSource(filename,
    decoder_ = new NectarCamCameraEventDecoder(decoder_config), false,
    config)
{
  // nothing to see here
}

NectarCamZFITSDataSource::~NectarCamZFITSDataSource()
{
  delete decoder_;
}
