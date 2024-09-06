/*

   calin/iact_data/acada_event_decoder.hpp -- Stephen Fegan -- 2024-09-05

   Base for all decoders of ACADA event types

   Copyright 2024, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#pragma once

#include <string>

#include <calin_global_definitions.hpp>
#include <calin_global_config.hpp>
#include <io/chained_data_source.hpp>
#include <iact_data/acada_data_source.hpp>
#include <iact_data/telescope_data_source.hpp>

namespace calin { namespace iact_data { namespace acada_event_decoder {

#ifdef CALIN_HAVE_CTA_CAMERASTOACTL

void decode_cdts_data(calin::ix::iact_data::telescope_event::CDTSData* calin_cdts_data,
  const DataModel::AnyArray& cta_array);

void decode_tib_data(calin::ix::iact_data::telescope_event::TIBData* calin_tib_data,
  const DataModel::AnyArray& cta_array);

calin::ix::iact_data::telescope_event::TriggerType determine_trigger_type(
    const calin::ix::iact_data::telescope_event::TIBData* calin_tib_data,
    const calin::ix::iact_data::telescope_event::CDTSData* calin_cdts_data);

template<typename EventMessage, typename HeaderMessage>
class ACADACameraEventDecoder
{
public:
  CALIN_TYPEALIAS(event_type, EventMessage);
  CALIN_TYPEALIAS(header_type, HeaderMessage);

  virtual ~ACADACameraEventDecoder();
  virtual bool decode(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    const EventMessage* cta_event) = 0;
  virtual bool decode_run_config(
    calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config,
    const HeaderMessage* cta_run_header,
    const EventMessage* cta_event) = 0;
  virtual ACADACameraEventDecoder* clone() const = 0;
};

template<typename EventMessage, typename HeaderMessage>
class DecodedACADACameraEventDataSource:
  public calin::iact_data::telescope_data_source::TelescopeDataSource
{
public:
  CALIN_TYPEALIAS(event_type, EventMessage);
  CALIN_TYPEALIAS(header_type, HeaderMessage);

  CALIN_TYPEALIAS(ACADACameraEventDataSource, 
    calin::iact_data::acada_data_source::ACADACameraEventDataSource<EventMessage>);
  CALIN_TYPEALIAS(ACADACameraEventDecoder,
    ACADACameraEventDecoder<EventMessage,HeaderMessage>);

  DecodedACADACameraEventDataSource(
    ACADACameraEventDataSource* acada_src, ACADACameraEventDecoder* decoder,
    bool adopt_acada_src = false, bool adopt_decoder = false);

  virtual ~DecodedACADACameraEventDataSource();

  calin::ix::iact_data::telescope_event::TelescopeEvent* get_next(
    uint64_t& seq_index_out, google::protobuf::Arena** arena = nullptr) override;

protected:
  ACADACameraEventDecoder* decoder_;
  virtual const EventMessage* borrow_next_acada_event(uint64_t& seq_index_out);
  virtual void release_borrowed_acada_event(const EventMessage* event);

private:
  bool adopt_decoder_ = false;
  ACADACameraEventDataSource* acada_src_ = nullptr;
  bool adopt_acada_src_ = false;
};

template<typename EventMessage, typename HeaderMessage>
class DecodedACADACameraEventDataSourceWithRunConfig:
  public virtual calin::iact_data::telescope_data_source::TelescopeDataSourceWithRunConfig,
  public virtual DecodedACADACameraEventDataSource<EventMessage,HeaderMessage>
{
public:
  CALIN_TYPEALIAS(event_type, EventMessage);
  CALIN_TYPEALIAS(header_type, HeaderMessage);

  CALIN_TYPEALIAS(ACADACameraEventDataSourceWithRunHeader, 
    calin::iact_data::acada_data_source::ACADACameraEventDataSourceWithRunHeader<EventMessage,HeaderMessage>);
  CALIN_TYPEALIAS(ACADACameraEventDecoder,
    ACADACameraEventDecoder<EventMessage,HeaderMessage>);

  DecodedACADACameraEventDataSourceWithRunConfig(
    ACADACameraEventDataSourceWithRunHeader* acada_src, ACADACameraEventDecoder* decoder,
    bool adopt_acada_src = false, bool adopt_decoder = false);

  virtual ~DecodedACADACameraEventDataSourceWithRunConfig();

  calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* get_run_configuration() override;

protected:
  const EventMessage* borrow_next_acada_event(uint64_t& seq_index_out) override;
  void ensure_run_config();

private:
  ACADACameraEventDataSourceWithRunHeader* acada_src_ = nullptr;
  const EventMessage* saved_event_ = nullptr;
  uint64_t saved_seq_index_ = 0;
  calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config_ = nullptr;
};

#endif

} } } // namespace calin::iact_data::acada_event_decoder
