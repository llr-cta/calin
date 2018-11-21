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
#include <io/zmq_data_source.pb.h>
#include <provenance/chronicle.hpp>

namespace calin { namespace io { namespace data_source {

template<typename T> class ZMQProtobufDataSource: public DataSource<T>
{
public:
  CALIN_TYPEALIAS(config_type, calin::ix::io::zmq_data_source::ZMQDataSourceConfig);
  CALIN_TYPEALIAS(data_type, typename DataSource<T>::data_type);

  ZMQProtobufDataSource(const ZMQProtobufDataSource&) = delete;
  ZMQProtobufDataSource& operator=(const ZMQProtobufDataSource&) = delete;

  ZMQProtobufDataSource(zmq_inproc::ZMQPuller* puller,
      calin::ix::io::zmq_data_source::ZMQDataSourceConfig config = default_config(),
      bool adopt_puller = false):
    DataSource<T>(), puller_(puller), adopt_puller_(adopt_puller),
     config_(config)
  {
    // nothing to see here
  }

  ZMQProtobufDataSource(const std::string& endpoint, void* zmq_ctx = nullptr,
      calin::ix::io::zmq_data_source::ZMQDataSourceConfig config = default_config(),
      calin::io::zmq_inproc::ZMQBindOrConnect bind_or_connect = calin::io::zmq_inproc::ZMQBindOrConnect::CONNECT):
    DataSource<T>(), adopt_puller_(true), config_(config)
  {
    if(zmq_ctx == nullptr) {
      my_zmq_ctx_ = zmq_ctx = calin::io::zmq_inproc::new_zmq_ctx();
      zmq_ctx_set(my_zmq_ctx_, ZMQ_IO_THREADS, config_.num_io_threads());
    }

    puller_ = new calin::io::zmq_inproc::ZMQPuller(zmq_ctx, endpoint,
      config_.receive_buffer_size(), bind_or_connect);

    network_io_record_ =
      calin::provenance::chronicle::register_network_open(endpoint, __PRETTY_FUNCTION__);
    network_io_record_->set_nbytes_sent(0);
  }

  ZMQProtobufDataSource(const std::string& endpoint,
      calin::ix::io::zmq_data_source::ZMQDataSourceConfig config,
      void* zmq_ctx = nullptr,
      calin::io::zmq_inproc::ZMQBindOrConnect bind_or_connect = calin::io::zmq_inproc::ZMQBindOrConnect::CONNECT):
    ZMQProtobufDataSource(endpoint, zmq_ctx, config, bind_or_connect)
  {
    // nothing to see here
  }

  virtual ~ZMQProtobufDataSource()
  {
    if(network_io_record_)
      calin::provenance::chronicle::register_network_close(network_io_record_,
        puller_->nbytes_received(), 0);
    if(adopt_puller_)delete puller_;
    if(my_zmq_ctx_)calin::io::zmq_inproc::destroy_zmq_ctx(my_zmq_ctx_);
  }

  T* get_next(uint64_t& seq_index_out, google::protobuf::Arena** arena = nullptr) override
  {
    if(arena != nullptr)
      throw std::runtime_error(
        "ZMQProtobufDataSource::get_next : arenas not supported");

    timed_out_ = false;
    long timeout = (seq_index_ == 0) ? config_.initial_receive_timeout_ms() : config_.receive_timeout_ms();
    if(puller_->wait_for_data(timeout) == false) {
      timed_out_ = true;
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
    if(network_io_record_) {
      network_io_record_->set_nbytes_received(puller_->nbytes_received());
    }
    return next;
  }

  bool timed_out() const { return timed_out_; }

  void* my_zmq_ctx() { return my_zmq_ctx_; }

  static calin::ix::io::zmq_data_source::ZMQDataSourceConfig default_config() {
    calin::ix::io::zmq_data_source::ZMQDataSourceConfig config;
    config.set_receive_timeout_ms(-1);
    config.set_initial_receive_timeout_ms(-1);
    config.set_receive_buffer_size(1000);
    config.set_num_io_threads(1);
    return config;
  }

private:
  void* my_zmq_ctx_ = nullptr;
  calin::io::zmq_inproc::ZMQPuller* puller_ = nullptr;
  bool adopt_puller_ = false;
  uint64_t seq_index_ = 0;
  calin::ix::io::zmq_data_source::ZMQDataSourceConfig config_ = default_config();
  bool timed_out_ = false;
  calin::ix::provenance::chronicle::NetworkIORecord* network_io_record_ = nullptr;
};


} } }
