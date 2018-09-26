/*

   calin/iact_data/nectarcam_actl_event_decoder.cpp -- Stephen Fegan -- 2018-09-21

   A decoder of NectarCAM ACTL data

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

#include <util/log.hpp>
#include <util/file.hpp>
#include <iact_data/nectarcam_actl_event_decoder.hpp>
#include <iact_data/nectarcam_layout.hpp>
#include <iact_data/nectarcam_module_configuration.hpp>
#include <math/simd.hpp>
#include <provenance/system_info.hpp>

using namespace calin::iact_data::nectarcam_actl_event_decoder;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::ix::iact_data::nectarcam_data_source;
using namespace calin::util::log;

#include <ProtobufIFits.h>
#include <R1.pb.h>

#define TEST_ANYARRAY_TYPES 0

/*

                      RRRRRRRRRRRRRRRRR          1111111
                      R::::::::::::::::R        1::::::1
                      R::::::RRRRRR:::::R      1:::::::1
                      RR:::::R     R:::::R     111:::::1
                        R::::R     R:::::R        1::::1
                        R::::R     R:::::R        1::::1
                        R::::RRRRRR:::::R         1::::1
                        R:::::::::::::RR          1::::l
                        R::::RRRRRR:::::R         1::::l
                        R::::R     R:::::R        1::::l
                        R::::R     R:::::R        1::::l
                        R::::R     R:::::R        1::::l
                      RR:::::R     R:::::R     111::::::111
                      R::::::R     R:::::R     1::::::::::1
                      R::::::R     R:::::R     1::::::::::1
                      RRRRRRRR     RRRRRRR     111111111111

*/

NectarCAM_ACTL_R1_CameraEventDecoder::NectarCAM_ACTL_R1_CameraEventDecoder(
    const std::string& filename, unsigned run_number, const config_type& config):
  actl_event_decoder::ACTL_R1_CameraEventDecoder(), config_(config),
  filename_(filename), run_number_(run_number)
{
  if(config_.demand_configured_module_id_size() != 0)
    LOG(WARNING) << "Decoder option \"demand_configured_module_id_size\" not supported in R1 data at this time.";

  if(config_.exchange_gain_channels() ==
      ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig::EXCHANGE_GAIN_MODE_FORCED)
    LOG(WARNING) << "Decoder option \"exchange_gain_channels\" not supported in R1 data at this time.";
}

NectarCAM_ACTL_R1_CameraEventDecoder::~NectarCAM_ACTL_R1_CameraEventDecoder()
{
  // nothing to see here
}

bool NectarCAM_ACTL_R1_CameraEventDecoder::decode(
  calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
  const R1::CameraEvent* cta_event)
{
  if(!cta_event->has_nectarcam())
    throw(std::runtime_error("NectarCAM_ACTL_R1_CameraEventDecoder::decode: "
      "ACTL event does not have NectarCAM extension"));

  calin_event->set_telescope_id(telescope_id_);
  calin_event->set_local_event_number(cta_event->tel_event_id());
  calin_event->set_trigger_type(TRIGGER_SCIENCE);
  calin_event->set_array_trigger_received(false);
  calin_event->set_array_event_number(cta_event->event_id());
  //calin_event->local_clock_time
  calin_event->set_image_treatment(TREATMENT_SCIENCE);
  calin_event->set_configuration_id(cta_event->configuration_id());
  calin_event->set_pedestal_dataset_id(cta_event->ped_id());

  calin_event->mutable_absolute_event_time()->set_time_ns(
    uint64_t(cta_event->trigger_time_s())*uint64_t(1000000000)
    + uint64_t(cta_event->trigger_time_qns()>>2));
  calin_event->mutable_elapsed_event_time()->set_time_ns(
    calin_event->absolute_event_time().time_ns() - run_start_time_);

  bool all_modules_present = true;
  if(cta_event->nectarcam().has_module_status())
  {
    const auto& cta_status = cta_event->nectarcam().module_status();
#if TEST_ANYARRAY_TYPES
    if(cta_status.type() != DataModel::AnyArray::U8)
      throw std::runtime_error("NectarCAM_ACTL_R1_CameraEventDecoder::decode: "
        "Camera status type not U8");
#endif
    if(nmod_ != cta_status.data().size())
      throw std::runtime_error("NectarCAM_ACTL_R1_CameraEventDecoder::decode: "
        "Module status array size does not match number of nodukes.");
    const auto* mod_status =
      reinterpret_cast<const uint8_t*>(&cta_status.data().front());
    for(unsigned imod=0, mod_index=0;imod<nmod_;imod++)
    {
      if(*(mod_status++)&0x01)
      {
        calin_event->add_module_index(mod_index);
        calin_event->add_module_id(imod);
        mod_index++;
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
    throw(std::runtime_error("NectarCAM_ACTL_R1_CameraEventDecoder::decode: "
      "ACTL event does not have NectarCAM module_status"));
  }
  calin_event->set_all_modules_present(all_modules_present);

  // ==========================================================================
  //
  // TRANSFER IMAGE DATA
  //
  // ==========================================================================

  if(cta_event->has_waveform() and cta_event->has_pixel_status())
  {
    unsigned single_gain_dataset_size = 7*nmod_*nsample_*sizeof(int16_t);
    if(cta_event->waveform().data().size() != 2*single_gain_dataset_size)
      throw(std::runtime_error("NectarCAM_ACTL_R1_CameraEventDecoder::decode: "
        "Samples array incorrect size: "
        + std::to_string(cta_event->waveform().data().size())
        + ", expected: " + std::to_string(2*single_gain_dataset_size)));
    if(cta_event->pixel_status().data().size() != nmod_*7)
      throw(std::runtime_error("NectarCAM_ACTL_R1_CameraEventDecoder::decode: "
        "Pixel status array incorrect size: "
        + std::to_string(cta_event->pixel_status().data().size())
        + ", expected: " + std::to_string(nmod_*7)));

    const uint8_t* pix_status =
      reinterpret_cast<const uint8_t*>(cta_event->pixel_status().data().data());
    const int16_t* waveforms =
      reinterpret_cast<const int16_t*>(cta_event->waveform().data().data());

    copy_single_gain_waveforms(calin_event, waveforms, pix_status,
      calin_event->mutable_high_gain_image()->mutable_camera_waveforms(),
      0x08, "high");

    copy_single_gain_waveforms(calin_event, waveforms+nmod_*nsample_, pix_status,
      calin_event->mutable_low_gain_image()->mutable_camera_waveforms(),
      0x04, "low");
  }
#if 0

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
#endif

  // ==========================================================================
  //
  // DECODE NECTARCAM COUNTERS
  //
  // ==========================================================================

  if(cta_event->nectarcam().has_counters())
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

    const auto& cta_counters = cta_event->nectarcam().counters();
#if TEST_ANYARRAY_TYPES
    if(cta_counters.type() != DataModel::AnyArray::U16)
      throw std::runtime_error("Camera counters type not U16");
#endif
    if(cta_counters.data().size()%sizeof(NectarCounters) != 0)
      throw std::runtime_error("Camera counters data array not integral "
        "multiple of expected structure size.");
    unsigned nmod = cta_counters.data().size()/sizeof(NectarCounters);
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

#define ts2_decode(x) int32_t(x)

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

  if(cta_event->nectarcam().has_cdts_data()
    and cta_event->nectarcam().cdts_data().has_data())
  {
    calin::iact_data::actl_event_decoder::decode_cdts_data(
      calin_event->mutable_cdts_data(), cta_event->nectarcam().cdts_data());

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

  if(cta_event->nectarcam().has_tib_data()
    and cta_event->nectarcam().tib_data().has_data())
  {
    calin::iact_data::actl_event_decoder::decode_tib_data(
      calin_event->mutable_tib_data(), cta_event->nectarcam().tib_data());
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

bool NectarCAM_ACTL_R1_CameraEventDecoder::decode_run_config(
  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* calin_run_config,
  const R1::CameraConfiguration* cta_run_header,
  const R1::CameraEvent* cta_event)
{
  calin_run_config->set_data_transcoder(
    "calin::iact_data::nectarcam_actl_event_decoder::NectarCAM_ACTL_R1_CameraEventDecoder");
  calin_run_config->set_filename(filename_);
  calin_run_config->add_fragment_filename(filename_);
  calin_run_config->set_run_number(run_number_);

  switch(config_.camera_type())
  {
  case NectarCamCameraEventDecoderConfig::AUTOMATIC:
  case NectarCamCameraEventDecoderConfig::NECTARCAM:
    nectarcam_layout::nectarcam_layout(
      calin_run_config->mutable_camera_layout());
    break;
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
  // Extract CTA run header data
  //
  // ---------------------------------------------------------------------------

  if(cta_run_header)
  {
    telescope_id_ = cta_run_header->telescope_id();
    calin_run_config->set_telescope_id(telescope_id_);
    calin_run_config->set_camera_server_serial_number(cta_run_header->cs_serial());
    calin_run_config->set_configuration_id(cta_run_header->configuration_id());
    calin_run_config->set_data_model_version(cta_run_header->data_model_version());
    run_start_time_ = uint64_t(cta_run_header->date()) * uint64_t(1000000000);
    calin_run_config->mutable_run_start_time()->set_time_ns(run_start_time_);
    if(cta_run_header->has_nectarcam()) {
      const auto* cta_nectarcam_header = &cta_run_header->nectarcam();
      auto* ncc = calin_run_config->mutable_nectarcam();
      if(ncc->daq_mode() == "") {
        switch(cta_nectarcam_header->acquisition_mode()) {
          // See MST-CAM-TN-0008-LPNHE
          case 0: ncc->set_daq_mode("DAQCHARGE_C"); break;
          case 1: ncc->set_daq_mode("DAQCHARGE_T"); break;
          case 2: ncc->set_daq_mode("DAQSAMPLE"); break;
          case 3: ncc->set_daq_mode("DAQSAMPLE_C"); break;
          case 4: ncc->set_daq_mode("DAQSAMPLE_D"); break;
          case 5: ncc->set_daq_mode("DAQCHARGESAMPLE"); break;
          case 6: ncc->set_daq_mode("DAQCHARGESAMPLE_D"); break;
          default:
            ncc->set_daq_mode("Unrecognised DAQ mode: "
              + std::to_string(cta_nectarcam_header->acquisition_mode()));
            break;
        }
      }
      ncc->set_daq_processing_algorithms(cta_nectarcam_header->algorithms());
      ncc->set_idaq_version(cta_nectarcam_header->idaq_version());
      ncc->set_cdhs_version(cta_nectarcam_header->cdhs_version());
    }
  }

  // ---------------------------------------------------------------------------
  //
  // Get the list of configured modules
  //
  // ---------------------------------------------------------------------------

  // First try to determine number of modules - require this to come from data
  nmod_ = 0;
  if(nmod_==0 and cta_run_header and cta_run_header->has_nectarcam()) {
    nmod_ = cta_run_header->nectarcam().num_modules();
  }
  if(nmod_==0 and cta_event and cta_event->has_nectarcam()
      and cta_event->nectarcam().has_module_status()) {
    nmod_ = cta_event->nectarcam().module_status().data().size();
  }
  if(nmod_ == 0) {
    throw std::runtime_error("NectarCAM_ACTL_R1_CameraEventDecoder::decode_run_config: "
      "Could not determine number of modules");
  }

  // Next get list of IDs
  unsigned nmod_camera = calin_run_config->camera_layout().module_size();
  std::vector<unsigned> config_mod_id;
  config_mod_id.reserve(nmod_camera);
  if(config_.demand_configured_module_id_size() != 0)
  {
    if(config_.demand_configured_module_id_size() != nmod_)
      throw std::runtime_error("NectarCAM_ACTL_R1_CameraEventDecoder::decode_run_config: "
        "Demand module list size must equal number of modules in data.");
    for(int imod=0;imod<nmod_;imod++) {
      unsigned mod_id = config_.demand_configured_module_id(imod);
      if(mod_id >= nmod_camera)
        throw std::runtime_error("NectarCAM_ACTL_R1_CameraEventDecoder::decode_run_config: "
          "Demand module id out of range: " + std::to_string(mod_id) + " >= " +
          std::to_string(nmod_camera));
      config_mod_id.push_back(mod_id);
    }
  }
  else if(calin_run_config->has_nectarcam() and
    calin_run_config->nectarcam().module_size() == nmod_)
  {
    for(int imod=0; imod<nmod_; imod++) {
      unsigned mod_id = calin_run_config->nectarcam().module(imod).module_id();
      if(mod_id >= nmod_camera)
        throw std::runtime_error("NectarCAM_ACTL_R1_CameraEventDecoder::decode_run_config: "
          "NMC module id out of range: " + std::to_string(mod_id) + " >= " +
          std::to_string(nmod_camera));
      config_mod_id.push_back(mod_id);
    }
  }
  else if(cta_run_header
    and cta_run_header->has_nectarcam()
    and cta_run_header->nectarcam().has_expected_modules_id()
    and cta_run_header->nectarcam().expected_modules_id().data().size() == nmod_*sizeof(uint16_t))
  {
    const uint16_t* mod_id =
      reinterpret_cast<const uint16_t*>(&
        cta_run_header->nectarcam().expected_modules_id().data().front());
    for(int imod=0;imod<nmod_;imod++)config_mod_id.push_back(mod_id[imod]);
  }
  else
  {
    for(int imod=0;imod<nmod_;imod++)config_mod_id.push_back(imod);
  }

  for(unsigned mod_id=0; mod_id<nmod_camera; mod_id++)
  {
    calin_run_config->add_configured_module_index(-1);
    for(unsigned ipix=0; ipix<7; ipix++)
      calin_run_config->add_configured_channel_index(-1);
  }
  for(unsigned imod=0; imod<config_mod_id.size();imod++)
  {
    unsigned mod_id = config_mod_id[imod];
    calin_run_config->add_configured_module_id(mod_id);
    calin_run_config->set_configured_module_index(mod_id, imod);
    for(unsigned ipix=0; ipix<7; ipix++) {
      calin_run_config->add_configured_channel_id(mod_id*7+ipix);
      calin_run_config->set_configured_channel_index(mod_id*7+ipix, imod*7+ipix);
    }
  }

  // ---------------------------------------------------------------------------
  //
  // Determine nsample
  //
  // ---------------------------------------------------------------------------

  nsample_ = config_.demand_nsample();
  if(nsample_ == 0 and cta_run_header)
    nsample_ = cta_run_header->num_samples();
  if(nsample_ == 0 and calin_run_config->has_nectarcam() and
      calin_run_config->nectarcam().module_size()>0) {
    nsample_ = calin_run_config->nectarcam().module(0).num_samples();
    for(int imod=1; imod<calin_run_config->nectarcam().module_size(); imod++)
      if(calin_run_config->nectarcam().module(imod).num_samples() != nsample_)
        nsample_ = 0;
  }
#if 0
  if(nsample == 0 and cta_event and cta_event->has_logain() and
      cta_event->logain().has_waveforms())
    nsample = cta_event->logain().waveforms().num_samples();
  if(nsample == 0 and cta_event and cta_event->has_higain() and
      cta_event->higain().has_waveforms())
    nsample = cta_event->higain().waveforms().num_samples();
#endif
  if(nsample_ == 0) {
    throw std::runtime_error("NectarCAM_ACTL_R1_CameraEventDecoder::decode_run_config: "
      "Could not determine number of samples");
  }
  calin_run_config->set_num_samples(nsample_);

  // ==========================================================================
  //
  // RUN SAMPLING FREQUENCY
  //
  // ==========================================================================

  double nominal_sampling_frequency = config_.demand_sampling_frequency();
  if(nominal_sampling_frequency == 0.0 and calin_run_config->has_nectarcam() and
      calin_run_config->nectarcam().module_size()>0) {
    nominal_sampling_frequency =
      calin_run_config->nectarcam().module(0).nominal_sampling_frequency();
    for(int imod=1; imod<calin_run_config->nectarcam().module_size(); imod++)
      if(calin_run_config->nectarcam().module(imod).nominal_sampling_frequency()
        != nominal_sampling_frequency)nominal_sampling_frequency = 0;
  }
  calin_run_config->set_nominal_sampling_frequency(nominal_sampling_frequency);

  // ==========================================================================
  //
  // RUN START TIME
  //
  // ==========================================================================

  if(cta_event->nectarcam().has_cdts_data()
    and cta_event->nectarcam().cdts_data().has_data())
  {
    calin::ix::iact_data::telescope_event::CDTSData calin_cdts_data;
    calin::iact_data::actl_event_decoder::decode_cdts_data(
      &calin_cdts_data, cta_event->nectarcam().cdts_data());

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
      SerializedRawHeaderType::SERIALIZED_RAW_HEADER_ACTL_R1_PROTOBUF);
    cta_run_header->SerializeToString(calin_run_config->mutable_serialized_raw_header());
  } else {
    calin_run_config->set_serialized_raw_header_type(
      SerializedRawHeaderType::SERIALIZED_RAW_HEADER_NONE);
  }
  return true;
}

void NectarCAM_ACTL_R1_CameraEventDecoder::
copy_single_gain_integrals(const DataModel::CameraEvent* cta_event,
  const calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
  const DataModel::PixelsChannel& cta_image,
  calin::ix::iact_data::telescope_event::DigitizedSkyImage* calin_image,
  const std::string& which_gain) const
{
#if 0
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
#endif
}

void NectarCAM_ACTL_R1_CameraEventDecoder::
copy_single_gain_waveforms(
  const calin::ix::iact_data::telescope_event::TelescopeEvent* calin_event,
  const int16_t* cta_waveforms, const uint8_t* cta_pixel_mask,
  calin::ix::iact_data::telescope_event::Waveforms* calin_waveforms,
  uint8 has_gain_mask, const std::string& which_gain) const
{
  unsigned npix = nmod_*7;
  bool all_channels_present = true;

  calin_waveforms->mutable_channel_index()->Reserve(npix);
  calin_waveforms->mutable_channel_id()->Reserve(npix);

  std::string* calin_wf_raw_data_string = calin_waveforms->mutable_raw_samples_array();
  unsigned simd_vec_size = calin::provenance::system_info::the_host_info()->simd_vec_size();
  calin_wf_raw_data_string->resize(nsample_*npix*sizeof(int16_t) + simd_vec_size);
  char* cp = &calin_wf_raw_data_string->front();
  std::fill(cp+nsample_*npix*sizeof(int16_t), cp+calin_wf_raw_data_string->size(), int8_t(0));
  int16_t* calin_wf_raw_data = reinterpret_cast<int16_t*>(cp);

  for(unsigned ipix=0;ipix<npix;ipix++)
  {
    if(cta_pixel_mask[ipix] & has_gain_mask) {
      std::copy(cta_waveforms, cta_waveforms+nsample_, calin_wf_raw_data);
      calin_waveforms->add_channel_index(calin_waveforms->channel_id_size());
      calin_waveforms->add_channel_id(ipix);
      if(config_.separate_channel_waveforms()) {
        auto* calin_samp = calin_waveforms->add_waveform()->mutable_samples();
        calin_samp->Reserve(nsample_);
        for(unsigned isample=0;isample<nsample_;isample++)
          calin_samp->Add(*cta_waveforms++);
      }
    } else {
      std::fill(calin_wf_raw_data, calin_wf_raw_data+nsample_, 0);
      all_channels_present = false;
      calin_waveforms->add_channel_index(-1);
      cta_waveforms += nsample_;
    }
  }

  calin_waveforms->set_all_channels_present(all_channels_present);
}
