/*

   calin/io/data_source.hpp -- Stephen Fegan -- 2016-01-15

   A supplier of generic data, probably Protobufs

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

#pragma once

namespace calin { namespace io { namespace data_source {

template<typename T> class DataSource
{
public:
  virtual ~DataSource();
  virtual T* getNext() = 0;
};

template<typename T> class DataSink
{
public:
  virtual ~DataSink();
  virtual bool putNext(T* d) = 0;
};

template<typename T> class ProtobufPacketDataSource
{
public:
  virtual ~ProtobufStreamDataSource();
  virtual T* getNext() = 0;
};

template<typename T> class ProtobufPacketDataSink
{
public:
  virtual ~ProtobufStreamDataSink();
  virtual bool getNext(T* d) = 0;
};

} } } // namespace calin::io::data_source
