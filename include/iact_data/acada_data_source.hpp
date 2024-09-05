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

using std::streamoff; // needed for IFits.h

#include <L0.pb.h>
#include <R1.pb.h>
#include <ProtobufIFits.h>

namespace calin { namespace iact_data { namespace acada_data_source {

template<typename EventMessage, typename HeaderMessage>
class ACADACameraEventDataSourceWithRunHeader:
  public calin::io::data_source::DataSource<EventMessage>
{
public:
  CALIN_TYPEALIAS(event_type, EventMessage);
  CALIN_TYPEALIAS(header_type, HeaderMessage);

  ACADACameraEventDataSourceWithRunHeader():
      calin::io::data_source::DataSource<event_type>() {
    /* nothing to see here */ }
  virtual ~ACADACameraEventDataSourceWithRunHeader();
  virtual header_type* get_run_header() = 0;
};

template<typename EventMessage, typename HeaderMessage>
class ACADACameraEventRandomAccessDataSourceWithRunHeader:
  public calin::io::data_source::RandomAccessDataSource<EventMessage>
{
public:
  CALIN_TYPEALIAS(event_type, EventMessage);
  CALIN_TYPEALIAS(header_type, HeaderMessage);

  ACADACameraEventRandomAccessDataSourceWithRunHeader();
  virtual ~ACADACameraEventRandomAccessDataSourceWithRunHeader();
  virtual header_type* get_run_header() = 0;

  virtual const event_type* borrow_next_event(uint64_t& seq_index_out) = 0;
  virtual void release_borrowed_event(const event_type* event) = 0;
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

CALIN_TYPEALIAS(ACADA_L0_EventMessage, DataModel::CameraEvent);
CALIN_TYPEALIAS(ACADA_L0_HeaderMessage, DataModel::CameraRunHeader);

// CALIN_TYPEALIAS(ACADA_L0_CameraEventDataSource,
//   calin::io::data_source::DataSource<ACADA_L0_EventMessage>);
// CALIN_TYPEALIAS(ACADA_L0_CameraEventRandomAccessDataSource,
//   calin::io::data_source::RandomAccessDataSource<ACADA_L0_EventMessage>);

// CALIN_TYPEALIAS(ConstACADA_L0_CameraEventDataSource,
//   calin::io::data_source::DataSource<const ACADA_L0_EventMessage>);
// CALIN_TYPEALIAS(ConstACADA_L0_CameraEventDataSink,
//   calin::io::data_source::DataSink<const ACADA_L0_EventMessage>);

// CALIN_TYPEALIAS(ACADA_L0_CameraEventDataSourceWithRunHeader,
//     ACADA_CameraEventDataSourceWithRunHeader<ACADA_L0_EventMessage,ACADA_L0_HeaderMessage>);
// CALIN_TYPEALIAS(ACADA_L0_CameraEventRandomAccessDataSourceWithRunHeader,
//     ACADA_CameraEventRandomAccessDataSourceWithRunHeader<ACADA_L0_EventMessage,ACADA_L0_HeaderMessage>);

// CALIN_TYPEALIAS(ACADA_L0_CameraEventChainedRandomAccessDataSourceWithRunHeader,
//   calin::io::data_source::BasicChainedRandomAccessDataSource<
//     ACADA_L0_CameraEventRandomAccessDataSourceWithRunHeader>);
// CALIN_TYPEALIAS(ACADA_L0_CameraEventRandomAccessDataSourceOpener,
//   calin::io::data_source::DataSourceOpener<
//     ACADA_L0_CameraEventRandomAccessDataSourceWithRunHeader>);

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

CALIN_TYPEALIAS(ACADA_R1v0_EventMessage, R1::CameraEvent);
CALIN_TYPEALIAS(ACADA_R1v0_HeaderMessage, R1::CameraConfiguration);

// CALIN_TYPEALIAS(ACADA_R1v0_CameraEventDataSource,
//   calin::io::data_source::DataSource<ACADA_R1v0_EventMessage>);
// CALIN_TYPEALIAS(ACADA_R1v0_CameraEventRandomAccessDataSource,
//   calin::io::data_source::RandomAccessDataSource<ACADA_R1v0_EventMessage>);

// CALIN_TYPEALIAS(ConstACADA_R1v0_CameraEventDataSource,
//   calin::io::data_source::DataSource<const ACADA_R1v0_EventMessage>);
// CALIN_TYPEALIAS(ConstACADA_R1v0_CameraEventDataSink,
//   calin::io::data_source::DataSink<const ACADA_R1v0_EventMessage>);

// CALIN_TYPEALIAS(ACADA_R1v0_CameraEventDataSourceWithRunHeader,
//     ACADA_CameraEventDataSourceWithRunHeader<ACADA_R1v0_EventMessage,ACADA_R1v0_HeaderMessage>);
// CALIN_TYPEALIAS(ACADA_R1v0_CameraEventRandomAccessDataSourceWithRunHeader,
//     ACADA_CameraEventRandomAccessDataSourceWithRunHeader<ACADA_R1v0_EventMessage,ACADA_R1v0_HeaderMessage>);

// CALIN_TYPEALIAS(ACADA_R1v0_CameraEventChainedRandomAccessDataSourceWithRunHeader,
//   calin::io::data_source::BasicChainedRandomAccessDataSource<
//     ACADA_R1v0_CameraEventRandomAccessDataSourceWithRunHeader>);
// CALIN_TYPEALIAS(ACADA_R1v0_CameraEventRandomAccessDataSourceOpener,
//   calin::io::data_source::DataSourceOpener<
//     ACADA_R1v0_CameraEventRandomAccessDataSourceWithRunHeader>);

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


} } } // namespace calin::iact_data::zfits_actl_data_source

