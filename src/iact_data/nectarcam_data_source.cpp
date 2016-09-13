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
#include <iact_data/nectarcam_layout.hpp>

using namespace calin::iact_data::nectarcam_data_source;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::io::log;

#include <ProtobufIFits.h>
#include <L0.pb.h>

#define TEST_ANYARRAY_TYPES 0

NectarCamCameraEventDecoder::~NectarCamCameraEventDecoder()
{
  // nothing to see here
}

bool NectarCamCameraEventDecoder::decode(
  calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
  const DataModel::CameraEvent* cta_event)
{
  calin_event->set_telescope_id(cta_event->telescopeid());
  calin_event->set_local_event_number(cta_event->eventnumber());
  calin_event->set_trigger_type(TRIGGER_SCIENCE);
  calin_event->set_array_trigger_received(false);
  calin_event->set_array_event_number(-1);
  //calin_event->local_clock_time
  calin_event->set_image_treatment(TREATMENT_SCIENCE);

  bool all_modules_present = true;
  if(cta_event->has_drawerstatus() and
    cta_event->drawerstatus().has_status())
  {
    const auto& cta_status = cta_event->drawerstatus().status();
#if TEST_ANYARRAY_TYPES
    if(cta_status.type() != DataModel::AnyArray::U8)
      throw std::runtime_error("Camera status type not U8");
#endif
    unsigned nmod =
      cta_status.data().size();
    const auto* mod_status =
      reinterpret_cast<const uint8_t*>(&cta_status.data().front());
    for(unsigned imod=0, mod_index=0;imod<nmod;imod++)
    {
      if(*(mod_status++)&0x01)
      {
        calin_event->add_module_index(mod_index);
        calin_event->add_module_id(imod);
      }
      else
      {
        calin_event->add_module_index(-1);
        all_modules_present = false;
      }
    }
  }
  else
  {
    unsigned nmod = 0;
    if(cta_event->has_higain() and
      cta_event->higain().has_integrals() and
      cta_event->higain().integrals().has_gains())
    {
      const auto& cta_q = cta_event->higain().integrals().gains();
      nmod = cta_q.data().size()/sizeof(uint16_t)/7;
    }
    else if(cta_event->has_logain() and
      cta_event->logain().has_integrals() and
      cta_event->logain().integrals().has_gains())
    {
      const auto& cta_q = cta_event->logain().integrals().gains();
      nmod = cta_q.data().size()/sizeof(uint16_t)/7;
    }
  }
  calin_event->set_all_modules_present(all_modules_present);

  // ==========================================================================
  //
  // TRANSFER IMAGE DATA
  //
  // ==========================================================================

  if(config_.exchange_gain_channels())
  {
    if(cta_event->has_higain())
      copy_single_gain_image(cta_event, calin_event, cta_event->higain(),
        calin_event->mutable_low_gain_image(), "high");

    if(cta_event->has_logain())
      copy_single_gain_image(cta_event, calin_event, cta_event->logain(),
        calin_event->mutable_high_gain_image(), "low");
  }
  else
  {
    if(cta_event->has_higain())
      copy_single_gain_image(cta_event, calin_event, cta_event->higain(),
        calin_event->mutable_high_gain_image(), "high");

    if(cta_event->has_logain())
      copy_single_gain_image(cta_event, calin_event, cta_event->logain(),
        calin_event->mutable_low_gain_image(), "low");
  }

  // ==========================================================================
  //
  // DECODE NECTARCAM COUNTERS
  //
  // ==========================================================================

  if(cta_event->has_cameracounters() and
    cta_event->cameracounters().has_counters())
  {
    struct NectarCounters {
      uint32_t global_event_counter;
      uint16_t bunch_counter;
      uint16_t event_counter;
      uint32_t ts1;
      uint8_t  ts2_event;
      uint8_t  ts2_bunch;
      uint16_t ts2_empty;
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
    for(unsigned imod=0;imod<nmod;imod++, mod_counter++)
    {
      if(imod < calin_event->module_index_size() and
        calin_event->module_index(imod) == -1)continue;
      auto* module_counters = calin_event->add_module_counter();
      module_counters->set_module_id(imod);
#define add_mod_counter(id,val) \
      { \
        module_counters->add_counter_id(id); \
        module_counters->add_counter_value(val); \
      }
#if 0
        auto* counter = module_counters->add_counter(); \
        counter->set_counter_id(id); \
        counter->set_value(val); \
      }
#endif
      add_mod_counter(0, mod_counter->global_event_counter);
      add_mod_counter(1, mod_counter->bunch_counter);
      add_mod_counter(2, mod_counter->event_counter);
      add_mod_counter(3, mod_counter->ts1);
      add_mod_counter(4, mod_counter->ts2_bunch);
      add_mod_counter(5, mod_counter->ts2_event);
      add_mod_counter(6, mod_counter->ts2_empty);

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
    }
  }

  // ==========================================================================
  //
  // DECODE NECTARCAM CDTS DATA MESSAGE
  //
  // ==========================================================================

  if(cta_event->uctsdatapresence() and cta_event->has_uctsdata() and
    cta_event->uctsdata().has_data())
  {
    struct CDTSMessageData {
      uint32_t event_counter;
      uint32_t pps_counter;
      uint32_t clock_counter;
      uint64_t ucts_timestamp;
      uint64_t camera_timestamp;
      uint8_t trigger_type;
      uint8_t white_rabbit_status;
      uint8_t arbitrary_information;
    } __attribute__((__packed__));

    const auto& cta_cdts_data = cta_event->uctsdata().data();
#if TEST_ANYARRAY_TYPES
    if(cta_cdts_data.type() != DataModel::AnyArray::U32)
      throw std::runtime_error("CDTS counters type not U32");
#endif
    if(cta_cdts_data.data().size() != sizeof(CDTSMessageData))
      throw std::runtime_error("CDTS data array not expected size");
    const auto* cdts_data =
      reinterpret_cast<const CDTSMessageData*>(&cta_cdts_data.data().front());

    auto* calin_cdts_data = calin_event->mutable_cdts_data();
    calin_cdts_data->set_event_counter(cdts_data->event_counter);
    calin_cdts_data->set_pps_counter(cdts_data->pps_counter);
    calin_cdts_data->set_clock_counter(cdts_data->clock_counter);
    calin_cdts_data->set_ucts_timestamp(cdts_data->ucts_timestamp);
    calin_cdts_data->set_camera_timestamp(cdts_data->camera_timestamp);
    calin_cdts_data->set_trigger_type(cdts_data->trigger_type);
    calin_cdts_data->set_white_rabbit_status(cdts_data->white_rabbit_status);
    calin_cdts_data->set_arbitrary_information(cdts_data->arbitrary_information);
  }

  return true;
}

bool NectarCamCameraEventDecoder::decode_run_config(
  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* calin_run_config,
  const DataModel::CameraRunHeader* cta_run_header,
  const DataModel::CameraEvent* cta_event)
{
  nectarcam_layout::nectarcam_19module_layout(
    calin_run_config->mutable_camera_layout());

  if(cta_run_header)
  {
#if 0
    calin_run_config->set_run_number(cta_run_header->runnumber());
    calin_run_config->set_run_start_time(
      make_time_mjd_ns(cta_run_header->datemjd(), cta_run_header->timenanosec());
    calin_run_config->set_num_samples(cta_run_header->numtraces());
#endif
  }

  unsigned nmod = 0;
  if(cta_run_header)
  {

  }
  if(nmod==0 and cta_event)
  {
    nmod = get_nmod_from_event(cta_event);
  }

  for(unsigned imod=0;imod<nmod;imod++)
  {
    calin_run_config->add_configured_module_id(imod);
    calin_run_config->add_configured_module_index(imod);
    for(unsigned ipix=0; ipix<7; ipix++)
    {
      calin_run_config->add_configured_channel_id(imod*7+ipix);
      calin_run_config->add_configured_channel_index(imod*7+ipix);
    }
  }

  unsigned nsample = config_.demand_nsample();
  if(nsample == 0 and cta_run_header)
    nsample = cta_run_header->numtraces();
  if(nsample == 0 and cta_event and cta_event->has_logain() and
      cta_event->logain().has_waveforms())
    nsample = cta_event->logain().waveforms().num_samples();
  if(nsample == 0 and cta_event and cta_event->has_higain() and
      cta_event->higain().has_waveforms())
    nsample = cta_event->higain().waveforms().num_samples();
  calin_run_config->set_num_samples(nsample);

  return true;
}

void NectarCamCameraEventDecoder::
copy_single_gain_image(const DataModel::CameraEvent* cta_event,
  const calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
  const DataModel::PixelsChannel& cta_image,
  calin::ix::iact_data::telescope_event::DigitizedSkyImage* calin_image,
  const std::string& which_gain) const
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
    bool all_channels_present = true;
    for(unsigned ipix=0;ipix<npix;ipix++)
    {
      if(calin_event->module_index(ipix/7) == -1)
      {
        all_channels_present = false;
        calin_wf_image->add_channel_index(-1);
        cta_wf_data += nsample;
      }
      else
      {
        calin_wf_image->add_channel_index(calin_wf_image->channel_id_size());
        calin_wf_image->add_channel_id(ipix);
        auto* calin_samp = calin_wf_image->add_waveform()->mutable_samples();
        calin_samp->Reserve(nsample);
        for(unsigned isample=0;isample<nsample;isample++)
          calin_samp->Add(*cta_wf_data++);
      }
    }
    calin_wf_image->set_all_channels_present(all_channels_present);
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
    bool all_channels_present = true;
    for(unsigned ipix=0;ipix<npix;ipix++, cta_q_data++)
    {
      if(calin_event->module_index(ipix/7) == -1)
      {
        calin_q_image->add_channel_index(-1);
        all_channels_present = false;
      }
      else
      {
        calin_q_image->add_channel_index(calin_q_image->channel_id_size());
        calin_q_image->add_channel_id(ipix);
        calin_q_image->add_charge(*cta_q_data);
      }
    }
    calin_q_image->set_all_channels_present(all_channels_present);
  }
}

unsigned NectarCamCameraEventDecoder::
get_nmod_from_event(const DataModel::CameraEvent* cta_event) const
{
  unsigned nmod = 0;
  if(cta_event->has_drawerstatus() and
    cta_event->drawerstatus().has_status())
  {
    const auto& cta_status = cta_event->drawerstatus().status();
#if TEST_ANYARRAY_TYPES
    if(cta_status.type() != DataModel::AnyArray::U8)
      throw std::runtime_error("Camera status type not U8");
#endif
    nmod = cta_status.data().size();
  }
  else if(cta_event->has_higain() and cta_event->higain().has_integrals()
    and cta_event->higain().integrals().has_gains())
  {
    const auto& cta_q = cta_event->higain().integrals().gains();
#if TEST_ANYARRAY_TYPES
    if(cta_q.type() != DataModel::AnyArray::U16)
      throw std::runtime_error("Integral data type not uint16 in "
        "high gain channel.");
#endif
    nmod = cta_q.data().size()/sizeof(uint16_t)/7;
  }
  else if(cta_event->has_logain() and cta_event->logain().has_integrals()
   and cta_event->logain().integrals().has_gains())
  {
    const auto& cta_q = cta_event->logain().integrals().gains();
#if TEST_ANYARRAY_TYPES
    if(cta_q.type() != DataModel::AnyArray::U16)
      throw std::runtime_error("Integral data type not uint16 in "
        "low gain channel.");
#endif
    nmod = cta_q.data().size()/sizeof(uint16_t)/7;
  }
  else if(cta_event->has_higain() and cta_event->higain().has_waveforms()
    and cta_event->higain().waveforms().has_samples())
  {
    const auto& cta_wf = cta_event->higain().waveforms().samples();
#if TEST_ANYARRAY_TYPES
    if(cta_wf.type() != DataModel::AnyArray::U16)
      throw std::runtime_error("Waveform data type not uint16 in " +
        "high gain channel.");
#endif
    unsigned nsample = config_.demand_nsample();
    if(nsample == 0)nsample = cta_event->higain().waveforms().num_samples();
    if(nsample == 0)throw std::runtime_error("Number of samples is zero in "
      "high gain channel.");
    if(cta_wf.data().size() % (sizeof(uint16_t)*nsample) != 0)
      throw std::runtime_error("Waveform data array for high gain "
        "channel not integral multiple of nsample uint16.");
    nmod = cta_wf.data().size()/(sizeof(uint16_t)*nsample*7);
  }
  else if(cta_event->has_logain() and cta_event->logain().has_waveforms()
    and cta_event->logain().waveforms().has_samples())
  {
    const auto& cta_wf = cta_event->logain().waveforms().samples();
#if TEST_ANYARRAY_TYPES
    if(cta_wf.type() != DataModel::AnyArray::U16)
      throw std::runtime_error("Waveform data type not uint16 in " +
        "low gain channel.");
#endif
    unsigned nsample = config_.demand_nsample();
    if(nsample == 0)nsample = cta_event->logain().waveforms().num_samples();
    if(nsample == 0)throw std::runtime_error("Number of samples is zero in "
      "low gain channel.");
    if(cta_wf.data().size() % (sizeof(uint16_t)*nsample) != 0)
      throw std::runtime_error("Waveform data array for low gain "
        "channel not integral multiple of nsample uint16.");
    nmod = cta_wf.data().size()/(sizeof(uint16_t)*nsample*7);
  }
  return nmod;
}

NectarCamZFITSDataSource::
NectarCamZFITSDataSource(const std::string& filename,
  const config_type& config, const decoder_config_type& decoder_config):
  calin::iact_data::zfits_data_source::ZFITSDataSource(filename,
    decoder_ = new NectarCamCameraEventDecoder(decoder_config), false,
    config)
{
  // nothing to see here
}

NectarCamZFITSDataSource::
NectarCamZFITSDataSource(const std::string& filename,
  const decoder_config_type& decoder_config, const config_type& config):
    NectarCamZFITSDataSource(filename, config, decoder_config)
{
  // nothing to see here
}

NectarCamZFITSDataSource::~NectarCamZFITSDataSource()
{
  delete decoder_;
}
