/*

   calin/io/one_to_n_data_sink.hpp -- Stephen Fegan -- 2016-03-18

   A data sink that passes events to one of N data sources

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

#include <atomic>
#include <thread>
#include <cassert>

#include <util/spinlock.hpp>
#include <util/log.hpp>
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

  bool put_next(T* data, uint64_t seq_index,
    google::protobuf::Arena* arena = nullptr, bool adopt_data = false)
  {
    assert(adopt_data);
    Payload<T> payload { data, arena, seq_index };
    bool good = pusher_->push(&payload, sizeof(payload));
    if(!good) {
      if(arena)delete arena;
      else delete data;
    }
    return good;
  }

  bool put_nullptr(bool no_wait = true)
  {
    Payload<T> payload = { nullptr, nullptr, 0 };
    return pusher_->push(&payload, sizeof(payload), no_wait);
  }

protected:
  zmq_inproc::ZMQInprocPushPull* zmq_ = nullptr;
  zmq_inproc::ZMQPusher* pusher_ = nullptr;
};

} } } // namespace calin::io::data_source
