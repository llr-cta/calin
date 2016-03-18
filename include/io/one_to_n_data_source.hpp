/*

   calin/io/one_to_n_data_source.hpp -- Stephen Fegan -- 2016-03-18

   A data sink that passes events to one of N data sources

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

#include <atomic>
#include <thread>
#include <cassert>

#include <util/spinlock.hpp>
#include <io/log.hpp>
#include <io/data_source.hpp>
#include <io/zmq_inproc_push_pull.hpp>
#include <io/buffered_data_source.hpp>

namespace calin { namespace io { namespace data_source {

template<typename T> class OneToNDataSink: public DataSink<T>
{
public:
  CALIN_TYPEALIAS(data_type, typename DataSource<T>::data_type);

  OneToNDataSink(const OneToNDataSink&) = delete;
  OneToNDataSink& operator=(const OneToNDataSink&) = delete;

  OneToNDataSink(): DataSink<T>(),
    zmq_(new zmq_inproc::ZMQInprocPushPull), pusher_(zmq_->new_pusher())
  {
    // nothing to see here
  }

  ~OneToNDataSink()
  {
    delete pusher_;
    delete zmq_;
  }

  BufferedDataSource<T>* new_data_source() {
    return new BufferedDataSource<T>(zmq_->new_puller()); }

  bool put_next(T* data, bool adopt_data = false)
  {
    assert(adopt_data);
    bool good = pusher_->push(&data, sizeof(data));
    if(!good)delete data;
    return good;
  }

protected:
  zmq_inproc::ZMQInprocPushPull* zmq_ = nullptr;
  zmq_inproc::ZMQPusher* pusher_ = nullptr;
};

} } } // namespace calin::io::data_source
