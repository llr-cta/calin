/*

   calin/io/zmq_inproc_push_pull.cpp -- Stephen Fegan -- 2016-03-01

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

#include <stdexcept>

#include <io/zmq_inproc_push_pull.hpp>

using namespace calin::io::zmq_inproc;

ZMQPusher::ZMQPusher(void* zmq_ctx, const std::string& endpoint):
  socket_(zmq_socket(zmq_ctx, ZMQ_PUSH), &zmq_close)
{
  // Open the socket
  if(!socket_)
    throw std::runtime_error(std::string("ZMQPusher: error creating socket: ")
      + zmq_strerror(errno));

  // Connect socket to endpoint
  if(zmq_connect(socket_.get(), endpoint.c_str()) < 0)
    throw std::runtime_error(
      std::string("ZMQPusher: error connecting socket: ")
      + endpoint + "\n"
      + "ZMQPusher: " + zmq_strerror(errno));
}

bool ZMQPusher::push(void* data, unsigned size, bool dont_wait)
{
  if(zmq_send(socket_.get(), data, size, dont_wait ? ZMQ_DONTWAIT : 0) < 0)
  {
    if(dont_wait and errno == EAGAIN)return false;
    throw std::runtime_error(std::string("ZMQPusher: error sending data: ")
      + zmq_strerror(errno));
  }
  return true;
}

ZMQPuller::ZMQPuller(void* zmq_ctx, const std::string& endpoint,
  int buffer_size):
  socket_(zmq_socket(zmq_ctx, ZMQ_PULL), &zmq_close)
{
  // Open the socket
  if(!socket_)
    throw std::runtime_error(std::string("ZMQPuller: error creating socket: ")
      + zmq_strerror(errno));

  // Set the buffer size
  if(buffer_size > 0)
  {
    if(zmq_setsockopt(socket_.get(), ZMQ_RCVHWM,
        &buffer_size, sizeof(buffer_size)) < 0)
      throw std::runtime_error(std::string("ZMQPuller: error setting buffer size: ")
        + zmq_strerror(errno));
  }

  // Connect socket to endpoint
  if(zmq_bind(socket_.get(), endpoint.c_str()) < 0)
    throw std::runtime_error(
      std::string("ZMQPuller: error binding socket endpoint : ")
        + endpoint + "\n" +
        + "ZMQPuller: " + zmq_strerror(errno));


}

bool ZMQPuller::pull(void* data, unsigned buffer_size, unsigned& bytes_received,
   bool dont_wait)
{
  int nbytes = zmq_recv(socket_.get(), data, buffer_size,
    dont_wait?ZMQ_DONTWAIT:0);
  if(nbytes < 0)
  {
    bytes_received = 0;
    if(dont_wait and errno == EAGAIN)return false;
    throw std::runtime_error(std::string("ZMQPuller: error receiving data: ")
      + zmq_strerror(errno));
  }
  bytes_received = unsigned(nbytes);
  return true;
}

ZMQInprocPushPull::ZMQInprocPushPull()
{
  my_zmq_ctx_ = zmq_ctx_ = zmq_ctx_new();
}

ZMQInprocPushPull::~ZMQInprocPushPull()
{
  if(my_zmq_ctx_)zmq_ctx_destroy(my_zmq_ctx_);
}

ZMQPuller* ZMQInprocPushPull::new_puller()
{
  return new ZMQPuller(zmq_ctx(), address());
}

ZMQPusher* ZMQInprocPushPull::new_pusher()
{
  return new ZMQPusher(zmq_ctx(), address());
}

std::string ZMQInprocPushPull::address()
{
  return std::string("inproc://#")+std::to_string(address_index_);
}
