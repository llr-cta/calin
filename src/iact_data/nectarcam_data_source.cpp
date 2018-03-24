/*

   calin/iact_data/nectarcam_data_source.cpp -- Stephen Fegan -- 2016-01-11

   A supplier of single telescope data from NectarCam DAQ data files

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

#include <stdexcept>
#include <string>
#include <memory>
#include <numeric>

#if !defined(__clang__) and defined(__GNUC__) && __GNUC__<5
namespace std {

inline void *align( std::size_t alignment, std::size_t size,
                void *&ptr, std::size_t &space ) {
  std::uintptr_t pn = reinterpret_cast< std::uintptr_t >( ptr );
  std::uintptr_t aligned = ( pn + alignment - 1 ) & - alignment;
  std::size_t padding = aligned - pn;
  if ( space < size + padding ) return nullptr;
  space -= padding;
  return ptr = reinterpret_cast< void * >( aligned );
}

}
#endif

#include <util/log.hpp>
#include <util/file.hpp>
#include <iact_data/nectarcam_data_source.hpp>
#include <iact_data/zfits_data_source.hpp>
#include <iact_data/nectarcam_layout.hpp>
#include <iact_data/nectarcam_module_configuration.hpp>
#include <math/simd.hpp>
#include <provenance/system_info.hpp>

using namespace calin::iact_data::nectarcam_data_source;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::ix::iact_data::nectarcam_data_source;
using namespace calin::util::log;

#include <ProtobufIFits.h>
#include <L0.pb.h>

#define TEST_ANYARRAY_TYPES 0

NectarCamCameraEventDecoder::NectarCamCameraEventDecoder(
  const std::string& filename,
  unsigned run_number, const config_type& config):
  zfits_data_source::CTACameraEventDecoder(), config_(config),
  filename_(filename), run_number_(run_number)
{
  switch(config.exchange_gain_channels()) {
    case ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig::EXCHANGE_GAIN_MODE_NONE:
      exchange_gain_channels_ = false;
      break;
    case ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig::EXCHANGE_GAIN_MODE_FORCED:
      exchange_gain_channels_ = true;
      break;
    case ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig::EXCHANGE_GAIN_MODE_AUTOMATIC:
    default:
      exchange_gain_channels_ = run_number>=32 and run_number<621;
      if(exchange_gain_channels_)
        LOG(WARNING) << "High/Low gain exchange automatically configured.";
      break;
  }
}

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
  if(cta_event->has_modulestatus() and
    cta_event->modulestatus().has_status())
  {
    const auto& cta_status = cta_event->modulestatus().status();
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
#if 0 // Is this needed?
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
#endif
  }
  calin_event->set_all_modules_present(all_modules_present);

  // ==========================================================================
  //
  // TRANSFER IMAGE DATA
  //
  // ==========================================================================

  if(exchange_gain_channels_)
  {
    if(cta_event->has_higain()) {
      copy_single_gain_integrals(cta_event, calin_event, cta_event->higain(),
        calin_event->mutable_low_gain_image(), "high");
      copy_single_gain_waveforms(cta_event, calin_event, cta_event->higain(),
        calin_event->mutable_low_gain_image(), "high");
    }

    if(cta_event->has_logain()) {
      copy_single_gain_integrals(cta_event, calin_event, cta_event->logain(),
        calin_event->mutable_high_gain_image(), "low");
      copy_single_gain_waveforms(cta_event, calin_event, cta_event->logain(),
        calin_event->mutable_high_gain_image(), "low");
    }
  }
  else
  {
    if(cta_event->has_higain()) {
      copy_single_gain_integrals(cta_event, calin_event, cta_event->higain(),
        calin_event->mutable_high_gain_image(), "high");
      copy_single_gain_waveforms(cta_event, calin_event, cta_event->higain(),
        calin_event->mutable_high_gain_image(), "high");
    }

    if(cta_event->has_logain()) {
      copy_single_gain_integrals(cta_event, calin_event, cta_event->logain(),
        calin_event->mutable_low_gain_image(), "low");
      copy_single_gain_waveforms(cta_event, calin_event, cta_event->logain(),
        calin_event->mutable_low_gain_image(), "low");
    }
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
      int8_t   ts2_event;
      int8_t   ts2_bunch;
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
      if(imod < static_cast<unsigned>(calin_event->module_index_size()) and
        calin_event->module_index(imod) == -1)continue;

      auto* module_counters = calin_event->add_module_counter();
      module_counters->set_module_id(imod);
#define add_mod_counter(id,val) \
      { \
        module_counters->add_counter_id(id); \
        module_counters->add_counter_value(val); \
      }
      add_mod_counter(0, mod_counter->global_event_counter);
      add_mod_counter(1, mod_counter->bunch_counter);
      add_mod_counter(2, mod_counter->event_counter);
      add_mod_counter(3, mod_counter->ts1);
      add_mod_counter(4, mod_counter->ts2_bunch);
      add_mod_counter(5, mod_counter->ts2_event);
      add_mod_counter(6, mod_counter->ts2_empty);

      auto* module_data = calin_event->add_module_data()->mutable_nectarcam();
      module_data->set_module_id(imod);
      module_data->set_global_event_counter(mod_counter->global_event_counter);
      module_data->set_bunch_counter(mod_counter->bunch_counter);
      module_data->set_event_counter(mod_counter->event_counter);
      module_data->set_ts1(mod_counter->ts1);
      module_data->set_ts2_bunch(mod_counter->ts2_bunch);
      module_data->set_ts2_event(mod_counter->ts2_event);
      module_data->set_ts2_empty(mod_counter->ts2_empty);

#if 0 // OBSOLETE version of TS2 definition
#define ts2_decode(x) ((x)&0xF0?((x)&0xC0?((x)&0x80?0:1):((x)&0x20?2:3)):\
                                ((x)&0x0C?((x)&0x08?4:5):((x)&0x02?6:7)))
#else
#define ts2_decode(x) int32_t(x)
#endif
      int32_t ts2_bunch = ts2_decode(mod_counter->ts2_bunch);
      int32_t ts2_event = ts2_decode(mod_counter->ts2_event);
      int32_t ts = mod_counter->ts1*8 + ts2_event - ts2_bunch;

      module_data->set_bunch_event_time(ts);

      int64_t time_ns = mod_counter->bunch_counter*1000000000LL + ts;
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
    calin::iact_data::zfits_data_source::decode_cdts_data(
      calin_event->mutable_cdts_data(), cta_event->uctsdata().data());

    if(calin_event->cdts_data().white_rabbit_status() == 1) {
      auto* calin_clock = calin_event->add_camera_clock();
      calin_clock->set_clock_id(0);
      calin_clock->mutable_time()->set_time_ns(calin_event->cdts_data().ucts_timestamp());
    }
  }

  // ==========================================================================
  //
  // DECODE NECTARCAM TIB DATA MESSAGE
  //
  // ==========================================================================

  if(cta_event->tibdatapresence() and cta_event->has_tibdata() and
    cta_event->tibdata().has_data())
  {
    calin::iact_data::zfits_data_source::decode_tib_data(
      calin_event->mutable_tib_data(), cta_event->tibdata().data());
  }

  // ==========================================================================
  //
  // FIGURE OUT EVENT TIME
  //
  // ==========================================================================

  if(calin_event->has_cdts_data()) {
    calin_event->mutable_absolute_event_time()->set_time_ns(
      calin_event->cdts_data().ucts_timestamp());
  } else {
    // Now what cat? Now what?
  }

  if(calin_event->has_absolute_event_time() and run_start_time_!=0) {
    calin_event->mutable_elapsed_event_time()->set_time_ns(
      calin_event->absolute_event_time().time_ns() - run_start_time_);
  }

  // ==========================================================================
  //
  // SERIALIZE RAW DATA
  //
  // ==========================================================================

  if(config_.include_serialized_raw_data())
  {
    calin_event->set_serialized_raw_event_type(
      SerializedRawEventType::SERIALIZED_RAW_EVENT_ACTL_PROTOBUF);
    cta_event->SerializeToString(calin_event->mutable_serialized_raw_event());
  } else {
    calin_event->set_serialized_raw_event_type(
      SerializedRawEventType::SERIALIZED_RAW_EVENT_NONE);
  }

  return true;
}

bool NectarCamCameraEventDecoder::decode_run_config(
  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* calin_run_config,
  const DataModel::CameraRunHeader* cta_run_header,
  const DataModel::CameraEvent* cta_event)
{
  calin_run_config->set_filename(filename_);
  calin_run_config->add_fragment_filename(filename_);
  calin_run_config->set_run_number(run_number_);

  switch(config_.camera_type())
  {
  case NectarCamCameraEventDecoderConfig::NECTARCAM:
    nectarcam_layout::nectarcam_layout(
      calin_run_config->mutable_camera_layout());
    break;
  case NectarCamCameraEventDecoderConfig::AUTOMATIC:
  case NectarCamCameraEventDecoderConfig::NECTARCAM_TESTBENCH_19CHANNEL:
  default:
    nectarcam_layout::nectarcam_19module_layout(
      calin_run_config->mutable_camera_layout());
    break;
  }

  // ---------------------------------------------------------------------------
  //
  // Try to read the NectarCam module configuration XML file
  //
  // ---------------------------------------------------------------------------

  std::vector<std::string> nmc_file_tried;
  std::string nmc_file;

  if(not config_.demand_nmc_xml_file().empty()) {
    if(calin::util::file::is_readable(config_.demand_nmc_xml_file())) {
      nmc_file = config_.demand_nmc_xml_file();
    } else {
      nmc_file_tried.emplace_back(config_.demand_nmc_xml_file());
    }
  } else {
    std::string nmc_dirname = calin::util::file::dirname(filename_);
    if(nmc_dirname == ".") {
      nmc_dirname = "";
    } else {
      nmc_dirname += '/';
    }
    std::string nmc_basename = calin::util::file::basename(filename_);
    while(not nmc_basename.empty()) {
      std::string test_file = nmc_dirname + nmc_basename + config_.nmc_xml_suffix();
      if(calin::util::file::is_readable(test_file)) {
        nmc_file = test_file;
        break;
      } else {
        nmc_file_tried.emplace_back(test_file);
      }
      nmc_basename = calin::util::file::strip_extension(nmc_basename);
    }
  }

  if(not nmc_file.empty()) {
    calin::ix::iact_data::nectarcam_module_configuration::NectarCamCameraConfiguration* nccc =
      calin::iact_data::nectarcam_module_configuration::decode_nmc_xml_file(nmc_file);
    if(nccc) {
      calin_run_config->mutable_nectarcam()->CopyFrom(*nccc);
    } else {
      LOG(WARNING) << "Could not parse NectarCAM module configuration XML file "
        << nmc_file;
    }
  } else {
    auto logger = LOG(WARNING);
    logger << "Could not find NectarCAM module configuration XML file, tried:\n";
    for(auto try_fn : nmc_file_tried) {
      logger << "- " << try_fn << '\n';
    }
    logger << "Set the \"demand_nmc_xml_file\" decoder option if you wish to "
      "specify a different file.";
  }

  // ---------------------------------------------------------------------------
  //
  // Get the list of configured modules
  //
  // ---------------------------------------------------------------------------

  if(cta_run_header)
  {
#if 0
    calin_run_config->set_run_number(cta_run_header->runnumber());
    calin_run_config->set_run_start_time(
      make_time_mjd_ns(cta_run_header->datemjd(), cta_run_header->timenanosec());
    calin_run_config->set_num_samples(cta_run_header->numtraces());
#endif
  }

  int nmod = 0;
  std::set<unsigned> config_mod_id;
  if(cta_run_header)
  {
    // what to do here
  }
  if(nmod==0 and cta_event)
  {
    nmod = get_nmod_from_event(cta_event);
    for(int imod=0;imod<nmod;imod++)config_mod_id.insert(imod);
  }
  unsigned nmod_camera = calin_run_config->camera_layout().module_size();
  if(config_.demand_configured_module_id_size() != 0)
  {
    if(config_.demand_configured_module_id_size() != nmod)
      throw std::runtime_error("NectarCamCameraEventDecoder::decode_run_config: "
        "Demand module list size must equal number of modules in data.");
    config_mod_id.clear();
    for(int imod=0;imod<nmod;imod++) {
      unsigned mod_id = config_.demand_configured_module_id(imod);
      if(mod_id >= nmod_camera)
        throw std::runtime_error("NectarCamCameraEventDecoder::decode_run_config: "
          "Demand module id out of range: " + std::to_string(mod_id) + " >= " +
          std::to_string(nmod_camera));
      config_mod_id.insert(mod_id);
    }
  } else if(calin_run_config->has_nectarcam() and
      calin_run_config->nectarcam().module_size() == nmod) {
    config_mod_id.clear();
    for(int imod=0; imod<nmod; imod++) {
      unsigned mod_id = calin_run_config->nectarcam().module(imod).module_id();
      if(mod_id >= nmod_camera)
        throw std::runtime_error("NectarCamCameraEventDecoder::decode_run_config: "
          "Demand module id out of range: " + std::to_string(mod_id) + " >= " +
          std::to_string(nmod_camera));
      config_mod_id.insert(mod_id);
    }
  }

  for(unsigned mod_id=0, mod_index=0; mod_id<nmod_camera; mod_id++)
  {
    if(config_mod_id.find(mod_id) == config_mod_id.end()) {
      calin_run_config->add_configured_module_index(-1);
      for(unsigned ipix=0; ipix<7; ipix++)
        calin_run_config->add_configured_channel_index(-1);
    } else {
      calin_run_config->add_configured_module_id(mod_id);
      calin_run_config->add_configured_module_index(mod_index);
      for(unsigned ipix=0; ipix<7; ipix++) {
        calin_run_config->add_configured_channel_id(mod_id*7+ipix);
        calin_run_config->add_configured_channel_index(mod_index*7+ipix);
      }
      ++mod_index;
    }
  }

  unsigned nsample = config_.demand_nsample();
  if(nsample == 0 and cta_run_header)
    nsample = cta_run_header->numtraces();
  if(nsample == 0 and calin_run_config->has_nectarcam() and
      calin_run_config->nectarcam().module_size()>0) {
    nsample = calin_run_config->nectarcam().module(0).daq_num_samples();
    for(int imod=1; imod<calin_run_config->nectarcam().module_size(); imod++)
      if(calin_run_config->nectarcam().module(imod).daq_num_samples() != nsample)
        nsample = 0;
  }
  if(nsample == 0 and cta_event and cta_event->has_logain() and
      cta_event->logain().has_waveforms())
    nsample = cta_event->logain().waveforms().num_samples();
  if(nsample == 0 and cta_event and cta_event->has_higain() and
      cta_event->higain().has_waveforms())
    nsample = cta_event->higain().waveforms().num_samples();
  calin_run_config->set_num_samples(nsample);

  // ==========================================================================
  //
  // RUN START TIME
  //
  // ==========================================================================

  if(cta_event->uctsdatapresence() and cta_event->has_uctsdata() and
    cta_event->uctsdata().has_data())
  {
    calin::ix::iact_data::telescope_event::CDTSData calin_cdts_data;
    calin::iact_data::zfits_data_source::decode_cdts_data(
      &calin_cdts_data, cta_event->uctsdata().data());

    if(calin_cdts_data.white_rabbit_status() == 1) {
      run_start_time_ = calin_cdts_data.ucts_timestamp();
      calin_run_config->mutable_run_start_time()->set_time_ns(run_start_time_);
    }
  }

  // ==========================================================================
  //
  // SERIALIZE RAW DATA
  //
  // ==========================================================================

  if(cta_run_header and config_.include_serialized_raw_data())
  {
    calin_run_config->set_serialized_raw_header_type(
      SerializedRawHeaderType::SERIALIZED_RAW_HEADER_ACTL_PROTOBUF);
    cta_run_header->SerializeToString(calin_run_config->mutable_serialized_raw_header());
  } else {
    calin_run_config->set_serialized_raw_header_type(
      SerializedRawHeaderType::SERIALIZED_RAW_HEADER_NONE);
  }

  return true;
}

void NectarCamCameraEventDecoder::
copy_single_gain_integrals(const DataModel::CameraEvent* cta_event,
  const calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
  const DataModel::PixelsChannel& cta_image,
  calin::ix::iact_data::telescope_event::DigitizedSkyImage* calin_image,
  const std::string& which_gain) const
{
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

void NectarCamCameraEventDecoder::
copy_single_gain_waveforms(const DataModel::CameraEvent* cta_event,
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

    std::string* calin_wf_raw_data_string = calin_wf_image->mutable_raw_samples_array();
    unsigned simd_vec_size = 0;
#if defined(__SSE2__)
    if(calin::provenance::system_info::the_host_info()->cpu_has_sse2())simd_vec_size=128/8;
#endif
#if defined(__AVX__)
    if(calin::provenance::system_info::the_host_info()->cpu_has_avx())simd_vec_size=256/8;
#endif
#if defined(__AVX512F__)
    if(calin::provenance::system_info::the_host_info()->cpu_has_avx512f())simd_vec_size=512/8;
#endif
    calin_wf_raw_data_string->resize(cta_wf.data().size() + 2*simd_vec_size);
    char* cp = &calin_wf_raw_data_string->front();
    void* vp = cp;
    std::size_t space = calin_wf_raw_data_string->size();
    if(simd_vec_size>0) {
      if(std::align(simd_vec_size, cta_wf.data().size(), vp, space)==nullptr)
        throw std::runtime_error("NectarCamCameraEventDecoder: cannot align data.");
      std::fill(cp, (char*)vp, uint8_t(0));
      std::fill((char*)vp+cta_wf.data().size(), cp+calin_wf_raw_data_string->size(), uint8_t(0));
    }
    calin_wf_image->set_raw_samples_array_start((char*)vp - cp);
    uint16_t* calin_wf_raw_data = reinterpret_cast<uint16_t*>(vp);

    unsigned npix = cta_wf.data().size()/(sizeof(uint16_t)*nsample);
    const uint16_t* cta_wf_data =
      reinterpret_cast<const uint16_t*>(&cta_wf.data().front());
    bool all_channels_present = true;

    calin_wf_image->mutable_channel_index()->Reserve(npix);
    calin_wf_image->mutable_channel_id()->Reserve(npix);

    for(unsigned ipix=0;ipix<npix;ipix++)
    {
      constexpr unsigned nblock_copy = 16;
      if((ipix&(nblock_copy-1)) == 0) {
        unsigned icount = (std::min(ipix+nblock_copy,npix) - ipix)*nsample;
        std::copy(cta_wf_data, cta_wf_data+icount, calin_wf_raw_data);
        // std::iota(calin_wf_raw_data, calin_wf_raw_data+icount, uint16_t(ipix*nsample));
        calin_wf_raw_data += icount;
      }

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
        if(config_.separate_channel_waveforms())
        {
          auto* calin_samp = calin_wf_image->add_waveform()->mutable_samples();
          calin_samp->Reserve(nsample);
          for(unsigned isample=0;isample<nsample;isample++)
            calin_samp->Add(*cta_wf_data++);
        }
      }
    }
    calin_wf_image->set_all_channels_present(all_channels_present);
  }
}

unsigned NectarCamCameraEventDecoder::
get_nmod_from_event(const DataModel::CameraEvent* cta_event) const
{
  unsigned nmod = 0;
  if(cta_event->has_modulestatus() and
    cta_event->modulestatus().has_status())
  {
    const auto& cta_status = cta_event->modulestatus().status();
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

// =============================================================================
// =============================================================================
//
// NectarCamZFITSDataSource
//
// =============================================================================
// =============================================================================

NectarCamZFITSDataSource::
NectarCamZFITSDataSource(const std::string& filename,
  const config_type& config, const decoder_config_type& decoder_config):
  calin::iact_data::zfits_data_source::ZFITSDataSource(filename,
    decoder_ = new NectarCamCameraEventDecoder(filename,
      calin::util::file::extract_first_number_from_filename(filename),
      decoder_config), false /* we delete it! */, config)
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
