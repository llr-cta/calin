/*

   calin/iact_data/acada_data_source.hpp -- Stephen Fegan -- 2022-09-05

   Definitions of sources of "raw" ACADA data types

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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
#include <io/data_source.hpp>
#include <io/chained_data_source.hpp>

#include <ProtoDataModel.pb.h>
#include <ProtoR1.pb.h>
#include <R1v1.pb.h>

namespace calin { namespace iact_data { namespace acada_data_source {

template<typename EventMessage>
class ACADACameraEventDataSource:
  public virtual calin::io::data_source::DataSource<EventMessage>
{
public:
  CALIN_TYPEALIAS(event_type, EventMessage);

  ACADACameraEventDataSource(): 
      calin::io::data_source::DataSource<EventMessage>() {
    /* nothing to see here */ }
  virtual ~ACADACameraEventDataSource();

  virtual const EventMessage* borrow_next_event(uint64_t& seq_index_out) = 0;
  virtual void release_borrowed_event(const EventMessage* event) = 0;
};

template<typename MessageSet>
class ACADACameraEventDataSourceWithRunHeader:
  public virtual ACADACameraEventDataSource<typename MessageSet::event_type>
{
public:
  CALIN_TYPEALIAS(message_set_type, MessageSet);
  CALIN_TYPEALIAS(event_type, typename MessageSet::event_type);
  CALIN_TYPEALIAS(header_type, typename MessageSet::header_type);
  CALIN_TYPEALIAS(data_stream_type, typename MessageSet::data_stream_type);

  ACADACameraEventDataSourceWithRunHeader(): 
      ACADACameraEventDataSource<event_type>() {
    /* nothing to see here */ }
  virtual ~ACADACameraEventDataSourceWithRunHeader();
  virtual header_type* get_run_header() = 0;
  virtual data_stream_type* get_data_stream() = 0;
};

template<typename MessageSet>
class ACADACameraEventRandomAccessDataSourceWithRunHeader:
  public virtual calin::io::data_source::RandomAccessDataSource<typename MessageSet::event_type>,
  public virtual ACADACameraEventDataSourceWithRunHeader<MessageSet>
{
public:
  CALIN_TYPEALIAS(message_set_type, MessageSet);
  CALIN_TYPEALIAS(event_type, typename MessageSet::event_type);
  CALIN_TYPEALIAS(header_type, typename MessageSet::header_type);
  CALIN_TYPEALIAS(data_stream_type, typename MessageSet::data_stream_type);

  ACADACameraEventRandomAccessDataSourceWithRunHeader(): 
      calin::io::data_source::RandomAccessDataSource<event_type>(), 
      ACADACameraEventDataSourceWithRunHeader<MessageSet>() {
    /* nothing to see here */ }
  virtual ~ACADACameraEventRandomAccessDataSourceWithRunHeader();
};

/*

              LLLLLLLLLLL                       000000000
              L:::::::::L                     00:::::::::00
              L:::::::::L                   00:::::::::::::00
              LL:::::::LL                  0:::::::000:::::::0
                L:::::L                    0::::::0   0::::::0
                L:::::L                    0:::::0     0:::::0
                L:::::L                    0:::::0     0:::::0
                L:::::L                    0:::::0 000 0:::::0
                L:::::L                    0:::::0 000 0:::::0
                L:::::L                    0:::::0     0:::::0
                L:::::L                    0:::::0     0:::::0
                L:::::L         LLLLLL     0::::::0   0::::::0
              LL:::::::LLLLLLLLL:::::L     0:::::::000:::::::0
              L::::::::::::::::::::::L      00:::::::::::::00
              L::::::::::::::::::::::L        00:::::::::00
              LLLLLLLLLLLLLLLLLLLLLLLL          000000000

*/

CALIN_TYPEALIAS(ACADA_EventMessage_L0, ProtoDataModel::CameraEvent);
CALIN_TYPEALIAS(ACADA_HeaderMessage_L0, ProtoDataModel::CameraRunHeader);

struct ACADA_MessageSet_L0 {
  CALIN_TYPEALIAS(event_type, ACADA_EventMessage_L0);
  CALIN_TYPEALIAS(header_type, ACADA_HeaderMessage_L0);
  CALIN_TYPEALIAS(data_stream_type, void);

  const event_type* event = nullptr;
  const header_type* header = nullptr;
  const data_stream_type* data_stream = nullptr;
};

/*

    RRRRRRRRRRRRRRRRR     1111111                              000000000     
    R::::::::::::::::R   1::::::1                            00:::::::::00   
    R::::::RRRRRR:::::R 1:::::::1                          00:::::::::::::00 
    RR:::::R     R:::::R111:::::1                         0:::::::000:::::::0
      R::::R     R:::::R   1::::1vvvvvvv           vvvvvvv0::::::0   0::::::0
      R::::R     R:::::R   1::::1 v:::::v         v:::::v 0:::::0     0:::::0
      R::::RRRRRR:::::R    1::::1  v:::::v       v:::::v  0:::::0     0:::::0
      R:::::::::::::RR     1::::l   v:::::v     v:::::v   0:::::0 000 0:::::0
      R::::RRRRRR:::::R    1::::l    v:::::v   v:::::v    0:::::0 000 0:::::0
      R::::R     R:::::R   1::::l     v:::::v v:::::v     0:::::0     0:::::0
      R::::R     R:::::R   1::::l      v:::::v:::::v      0:::::0     0:::::0
      R::::R     R:::::R   1::::l       v:::::::::v       0::::::0   0::::::0
    RR:::::R     R:::::R111::::::111     v:::::::v        0:::::::000:::::::0
    R::::::R     R:::::R1::::::::::1      v:::::v          00:::::::::::::00 
    R::::::R     R:::::R1::::::::::1       v:::v             00:::::::::00   
    RRRRRRRR     RRRRRRR111111111111        vvv                000000000     

*/

CALIN_TYPEALIAS(ACADA_EventMessage_R1v0, ProtoR1::CameraEvent);
CALIN_TYPEALIAS(ACADA_HeaderMessage_R1v0, ProtoR1::CameraConfiguration);

struct ACADA_MessageSet_R1v0 {
  CALIN_TYPEALIAS(event_type, ACADA_EventMessage_R1v0);
  CALIN_TYPEALIAS(header_type, ACADA_HeaderMessage_R1v0);
  CALIN_TYPEALIAS(data_stream_type, void);

  const event_type* event = nullptr;
  const header_type* header = nullptr;
  const data_stream_type* data_stream = nullptr;
};

/*

        RRRRRRRRRRRRRRRRR     1111111                        1111111   
        R::::::::::::::::R   1::::::1                       1::::::1   
        R::::::RRRRRR:::::R 1:::::::1                      1:::::::1   
        RR:::::R     R:::::R111:::::1                      111:::::1   
          R::::R     R:::::R   1::::1vvvvvvv           vvvvvvv1::::1   
          R::::R     R:::::R   1::::1 v:::::v         v:::::v 1::::1   
          R::::RRRRRR:::::R    1::::1  v:::::v       v:::::v  1::::1   
          R:::::::::::::RR     1::::l   v:::::v     v:::::v   1::::l   
          R::::RRRRRR:::::R    1::::l    v:::::v   v:::::v    1::::l   
          R::::R     R:::::R   1::::l     v:::::v v:::::v     1::::l   
          R::::R     R:::::R   1::::l      v:::::v:::::v      1::::l   
          R::::R     R:::::R   1::::l       v:::::::::v       1::::l   
        RR:::::R     R:::::R111::::::111     v:::::::v     111::::::111
        R::::::R     R:::::R1::::::::::1      v:::::v      1::::::::::1
        R::::::R     R:::::R1::::::::::1       v:::v       1::::::::::1
        RRRRRRRR     RRRRRRR111111111111        vvv        111111111111
                                                               

*/

CALIN_TYPEALIAS(ACADA_EventMessage_R1v1, R1v1::Event);
CALIN_TYPEALIAS(ACADA_HeaderMessage_R1v1, R1v1::CameraConfiguration);
CALIN_TYPEALIAS(ACADA_DataStreamMessage_R1v1, R1v1::TelescopeDataStream);

struct ACADA_MessageSet_R1v1 {
  CALIN_TYPEALIAS(event_type, ACADA_EventMessage_R1v1);
  CALIN_TYPEALIAS(header_type, ACADA_HeaderMessage_R1v1);
  CALIN_TYPEALIAS(data_stream_type, ACADA_DataStreamMessage_R1v1);

  const event_type* event = nullptr;
  const header_type* header = nullptr;
  const data_stream_type* data_stream = nullptr;
};

} } } // namespace calin::iact_data::acada_data_source

