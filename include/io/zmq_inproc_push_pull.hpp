/*

   calin/io/zmq_inproc_push_pull.hpp -- Stephen Fegan -- 2016-03-01

   A class to implement ZMQ inprox push/pull sockets

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

#include <string>
#include <memory>

#include <zmq.h>

namespace calin { namespace io { namespace zmq_inproc {

class ZMQPusher
{
public:
  ZMQPusher(void* zmq_ctx, const std::string& endpoint, int buffer_size = 100);
  bool push(void* data, unsigned size, bool dont_wait = false);
private:
  std::unique_ptr<void,int(*)(void*)> socket_;
};

class ZMQPuller
{
public:
  ZMQPuller(void* zmq_ctx, const std::string& endpoint, int buffer_size = 100);
  bool pull(void* data, unsigned buffer_size, unsigned& bytes_received,
     bool dont_wait = false);
  bool pull_assert_size(void* data, unsigned buffer_size,
      bool dont_wait = false);
private:
  std::unique_ptr<void,int(*)(void*)> socket_;
};

class ZMQInprocPushPull
{
public:
  ZMQInprocPushPull(unsigned buffer_size = 100);
  ~ZMQInprocPushPull();

  ZMQPuller* new_puller();
  ZMQPusher* new_pusher();

  void* zmq_ctx() { return zmq_ctx_; }
  unsigned address_index() { return address_index_; }
  std::string address();

private:
  unsigned buffer_size_ = 100;
  void* my_zmq_ctx_ = nullptr;
  void* zmq_ctx_ = nullptr;
  unsigned address_index_ = 0;
};

} } } // namespace calin::io::zmq_inproc
