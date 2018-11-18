/*

   calin/io/zmq_protobuf_data_source.hpp -- Stephen Fegan -- 2018-11-16

   A supplier of Protobufs, pulled and decoded from a ZMQ stream

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

#pragma once

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>

#include <io/data_source.hpp>
#include <io/zmq_inproc_push_pull.hpp>
#include <provenance/chronicle.hpp>

namespace calin { namespace io { namespace data_source {

template<typename T> class ZMQProtobufDataSource: public DataSource<T>
{
public:
  CALIN_TYPEALIAS(data_type, typename DataSource<T>::data_type);

  ZMQProtobufDataSource(const ZMQProtobufDataSource&) = delete;
  ZMQProtobufDataSource& operator=(const ZMQProtobufDataSource&) = delete;

  ZMQProtobufDataSource(zmq_inproc::ZMQPuller* puller,
    long timeout_ms = -1, long timeout_ms_zero = -1, bool adopt_puller = false):
    DataSource<T>(), puller_(puller), adopt_puller_(adopt_puller),
    timeout_ms_(timeout_ms), timeout_ms_zero_(timeout_ms_zero)
  {
    // nothing to see here
  }

  ZMQProtobufDataSource(void* zmq_ctx, const std::string& endpoint,
      long timeout_ms = -1, long timeout_ms_zero = -1, int buffer_size = 100,
      calin::io::zmq_inproc::ZMQBindOrConnect bind_or_connect = calin::io::zmq_inproc::ZMQBindOrConnect::CONNECT):
    DataSource<T>(), puller_(
      new zmq_inproc::ZMQPuller(zmq_ctx, endpoint, buffer_size, bind_or_connect)),
    adopt_puller_(true), timeout_ms_(timeout_ms), timeout_ms_zero_(timeout_ms_zero)
  {
    calin::provenance::chronicle::register_network_open(endpoint, __PRETTY_FUNCTION__);
  }

  virtual ~ZMQProtobufDataSource()
  {
    if(adopt_puller_)delete puller_;
  }

  T* get_next(uint64_t& seq_index_out, google::protobuf::Arena** arena = nullptr) override
  {
    if(arena != nullptr)
      throw std::runtime_error(
        "ZMQProtobufDataSource::get_next : arenas not supported");

    long timeout = (seq_index_ == 0) ? timeout_ms_zero_ : timeout_ms_;
    if(puller_->wait_for_data(timeout) == false) {
      return nullptr;
    }

    zmq_msg_t msg;
    zmq_msg_init(&msg);
    try {
      if(puller_->pull(&msg, true) == 0) {
        zmq_msg_close (&msg);
        return nullptr;
      }
    } catch(...) {
      zmq_msg_close (&msg);
      throw;
    }
    T* next = new T();
    if(!next->ParseFromArray(zmq_msg_data(&msg), zmq_msg_size(&msg))) {
      delete next;
      next = nullptr;
    } else {
      seq_index_out = seq_index_++;
    }
    zmq_msg_close (&msg);
    return next;
  }

private:
  zmq_inproc::ZMQPuller* puller_ = nullptr;
  bool adopt_puller_ = false;
  uint64_t seq_index_ = 0;
  long timeout_ms_ = -1;
  long timeout_ms_zero_ = -1;
};


} } }
