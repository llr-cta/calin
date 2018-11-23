/*

   calin/iact_data/cta_actl_r1_event_decoder.cpp -- Stephen Fegan -- 2018-11-23

   A supplier of single telescope data from CTA DAQ data files

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

#include <string>

#include <util/log.hpp>
#include <iact_data/cta_actl_event_decoder.hpp>
#include <iact_data/nectarcam_actl_event_decoder.hpp>
#include <iact_data/lstcam_actl_event_decoder.hpp>

using namespace calin::ix::iact_data::cta_data_source;
using namespace calin::iact_data::cta_actl_event_decoder;
using namespace calin::iact_data::nectarcam_actl_event_decoder;
using namespace calin::iact_data::lstcam_actl_event_decoder;
using namespace calin::pattern::delegation;

#include <ProtobufIFits.h>
#include <IFits.h>
#include <L0.pb.h>
#include <R1.pb.h>

CTA_ACTL_R1_CameraEventDecoder::
CTA_ACTL_R1_CameraEventDecoder(const std::string& filename, unsigned run_number,
    const calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig& config):
  ACTL_R1_CameraEventDecoder(), Delegator<ACTL_R1_CameraEventDecoder>(nullptr),
  filename_(filename), run_number_(run_number), config_(config)
{
  // nothing to see here
}

CTA_ACTL_R1_CameraEventDecoder::~CTA_ACTL_R1_CameraEventDecoder()
{
  // nothing to see here
}

bool CTA_ACTL_R1_CameraEventDecoder::decode(
  calin::ix::iact_data::telescope_event::TelescopeEvent* event,
  const R1::CameraEvent* cta_event)
{
  if(this->delegate()==nullptr) {
    if((config_.camera_type() == CTACameraEventDecoderConfig::NECTARCAM) or
      ((config_.camera_type() == CTACameraEventDecoderConfig::AUTO_DETECT) and
        (cta_event!=nullptr and cta_event->has_nectarcam())))
      this->set_delegate(new NectarCam_ACTL_R1_CameraEventDecoder(filename_, run_number_, config_.nectarcam()),true);
    else if((config_.camera_type() == CTACameraEventDecoderConfig::LSTCAM) or
      ((config_.camera_type() == CTACameraEventDecoderConfig::AUTO_DETECT) and
        (cta_event!=nullptr and cta_event->has_lstcam())))
      this->set_delegate(new LSTCam_ACTL_R1_CameraEventDecoder(filename_, run_number_, config_.lstcam()),true);
    else
      throw std::runtime_error("CTA_ACTL_R1_CameraEventDecoder: event does not "
        "have NectarCAM or LSTCam extensions");
  } else if(config_.camera_type() == CTACameraEventDecoderConfig::AUTO_DETECT
      and dynamic_cast<LSTCam_ACTL_R1_CameraEventDecoder*>(this->delegate())
      and cta_event!=nullptr and cta_event->has_nectarcam()) {
    this->set_delegate(new NectarCam_ACTL_R1_CameraEventDecoder(filename_, run_number_, config_.nectarcam()),true);
  } else if(config_.camera_type() == CTACameraEventDecoderConfig::AUTO_DETECT
      and dynamic_cast<NectarCam_ACTL_R1_CameraEventDecoder*>(this->delegate())
      and cta_event!=nullptr and cta_event->has_lstcam()) {
    this->set_delegate(new LSTCam_ACTL_R1_CameraEventDecoder(filename_, run_number_, config_.lstcam()),true);
  }
  return this->delegate()->decode(event, cta_event);
}

bool CTA_ACTL_R1_CameraEventDecoder::decode_run_config(
  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* run_config,
  const R1::CameraConfiguration* cta_run_header,
  const R1::CameraEvent* cta_event)
{
  if(this->delegate()==nullptr) {
    if((config_.camera_type() == CTACameraEventDecoderConfig::NECTARCAM) or
      ((config_.camera_type() == CTACameraEventDecoderConfig::AUTO_DETECT) and
        ((cta_event!=nullptr and cta_event->has_nectarcam()) or
        (cta_run_header!=nullptr and cta_run_header->has_nectarcam()))))
      this->set_delegate(new NectarCam_ACTL_R1_CameraEventDecoder(filename_, run_number_, config_.nectarcam()),true);
    else if((config_.camera_type() == CTACameraEventDecoderConfig::LSTCAM) or
      ((config_.camera_type() == CTACameraEventDecoderConfig::AUTO_DETECT) and
        ((cta_event!=nullptr and cta_event->has_lstcam()) or
        (cta_run_header!=nullptr and cta_run_header->has_lstcam()))))
      this->set_delegate(new LSTCam_ACTL_R1_CameraEventDecoder(filename_, run_number_, config_.lstcam()),true);
    else
      throw std::runtime_error("CTA_ACTL_R1_CameraEventDecoder: event does not "
        "have NectarCAM or LSTCam extensions");
  } else if(config_.camera_type() == CTACameraEventDecoderConfig::AUTO_DETECT
      and dynamic_cast<LSTCam_ACTL_R1_CameraEventDecoder*>(this->delegate())
      and ((cta_event!=nullptr and cta_event->has_nectarcam()) or
        (cta_run_header!=nullptr and cta_run_header->has_nectarcam()))) {
    this->set_delegate(new NectarCam_ACTL_R1_CameraEventDecoder(filename_, run_number_, config_.nectarcam()),true);
  } else if(config_.camera_type() == CTACameraEventDecoderConfig::AUTO_DETECT
      and dynamic_cast<NectarCam_ACTL_R1_CameraEventDecoder*>(this->delegate())
      and ((cta_event!=nullptr and cta_event->has_lstcam()) or
        (cta_run_header!=nullptr and cta_run_header->has_lstcam()))) {
    this->set_delegate(new LSTCam_ACTL_R1_CameraEventDecoder(filename_, run_number_, config_.lstcam()),true);
  }
  return this->delegate()->decode_run_config(run_config, cta_run_header, cta_event);
}

calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig
CTA_ACTL_R1_CameraEventDecoder::default_config()
{
  calin::ix::iact_data::cta_data_source::CTACameraEventDecoderConfig config;
  config.mutable_nectarcam()->CopyFrom(
    NectarCam_ACTL_R1_CameraEventDecoder::default_config());
  config.mutable_lstcam()->CopyFrom(
    LSTCam_ACTL_R1_CameraEventDecoder::default_config());
  return config;
}
