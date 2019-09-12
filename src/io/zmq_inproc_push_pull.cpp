/*

   calin/io/zmq_inproc_push_pull.cpp -- Stephen Fegan -- 2016-03-01

   A class to implement ZMQ inprox push/pull sockets

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

#include <stdexcept>

#include <io/zmq_inproc_push_pull.hpp>

using namespace calin::io::zmq_inproc;

void* calin::io::zmq_inproc::new_zmq_ctx()
{
  return zmq_ctx_new();
}

void calin::io::zmq_inproc::destroy_zmq_ctx(void* zmq_cxt)
{
  zmq_ctx_destroy(zmq_cxt);
}

ZMQPusher::ZMQPusher(void* zmq_ctx, const std::string& endpoint,
  int buffer_size, ZMQBindOrConnect bind_or_connect, ZMQProtocol protocol):
  socket_(zmq_socket(zmq_ctx, (protocol==ZMQProtocol::PUSH_PULL)?ZMQ_PUSH:ZMQ_PUB), &zmq_close)
{
  // Open the socket
  if(!socket_)
    throw std::runtime_error(std::string("ZMQPusher: error creating socket: ")
      + zmq_strerror(errno));

  // Set the buffer size
  if(buffer_size > 0)
  {
    if(zmq_setsockopt(socket_.get(), ZMQ_SNDHWM,
        &buffer_size, sizeof(buffer_size)) < 0)
      throw std::runtime_error(std::string("ZMQPusher: error setting buffer size: ")
        + zmq_strerror(errno));
  }

  // Bind oe connect socket to endpoint
  if(bind_or_connect == ZMQBindOrConnect::BIND)
  {
    if(zmq_bind(socket_.get(), endpoint.c_str()) < 0)
      throw std::runtime_error(
        std::string("ZMQPusher: error binding socket endpoint : ")
          + endpoint + "\n" +
          + "ZMQPusher: " + zmq_strerror(errno));
  }
  else
  {
    if(zmq_connect(socket_.get(), endpoint.c_str()) < 0)
      throw std::runtime_error(
        std::string("ZMQPusher: error connecting socket: ")
        + endpoint + "\n"
        + "ZMQPusher: " + zmq_strerror(errno));
  }
}

bool ZMQPusher::push(const std::string& data_push, bool dont_wait)
{
  return this->push(data_push.data(), data_push.size(), dont_wait);
}

bool ZMQPusher::push(const void* data, unsigned size, bool dont_wait)
{
  if(zmq_send(socket_.get(), data, size, dont_wait ? ZMQ_DONTWAIT : 0) < 0)
  {
    if((dont_wait and errno == EAGAIN) or errno == ETERM)return false;
    throw std::runtime_error(std::string("ZMQPusher: error sending data: ")
      + zmq_strerror(errno));
  }
  return true;
}

ZMQPuller::ZMQPuller(void* zmq_ctx, const std::string& endpoint,
  int buffer_size, ZMQBindOrConnect bind_or_connect, ZMQProtocol protocol):
  socket_(zmq_socket(zmq_ctx, (protocol==ZMQProtocol::PUSH_PULL)?ZMQ_PULL:ZMQ_SUB), &zmq_close)
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

  // Bind oe connect socket to endpoint
  if(bind_or_connect == ZMQBindOrConnect::BIND)
  {
    if(zmq_bind(socket_.get(), endpoint.c_str()) < 0)
      throw std::runtime_error(
        std::string("ZMQPuller: error binding socket endpoint : ")
          + endpoint + "\n" +
          + "ZMQPuller: " + zmq_strerror(errno));
  }
  else
  {
    if(zmq_connect(socket_.get(), endpoint.c_str()) < 0)
      throw std::runtime_error(
        std::string("ZMQPuller: error connecting socket: ")
        + endpoint + "\n"
        + "ZMQPuller: " + zmq_strerror(errno));
  }

  // If we are PUB/SUB then we must subscribe to all messages
  if(protocol==ZMQProtocol::PUB_SUB) {
    if(zmq_setsockopt(socket_.get(), ZMQ_SUBSCRIBE, "", 0) < 0)
      throw std::runtime_error("ZMQPuller: error setting PUB/SUB subscription: " + endpoint + "\n"
        + "ZMQPuller: " + zmq_strerror(errno));
  }
}

bool ZMQPuller::pull(zmq_msg_t* msg, bool dont_wait)
{
  int nbytes = zmq_recvmsg(socket_.get(), msg, dont_wait?ZMQ_DONTWAIT:0);
  if(nbytes < 0)
  {
    if((dont_wait and errno == EAGAIN) or errno == ETERM)return false;
    throw std::runtime_error(std::string("ZMQPuller: error receiving data: ")
      + zmq_strerror(errno));
  }
  nbytes_pulled_ += nbytes;
  return true;
}

bool ZMQPuller::pull(std::string& data_pull, bool dont_wait)
{
  /* Create an empty ØMQ message */
  zmq_msg_t msg;
  zmq_msg_init(&msg);
  pull(&msg, dont_wait);
  data_pull.assign(static_cast<char*>(zmq_msg_data(&msg)), zmq_msg_size(&msg));
  zmq_msg_close (&msg);
  return true;
}

bool ZMQPuller::pull(void* data, unsigned buffer_size, unsigned& bytes_received,
   bool dont_wait)
{
  int nbytes = zmq_recv(socket_.get(), data, buffer_size,
    dont_wait?ZMQ_DONTWAIT:0);
  if(nbytes < 0)
  {
    bytes_received = 0;
    if((dont_wait and errno == EAGAIN) or errno == ETERM)return false;
    throw std::runtime_error(std::string("ZMQPuller: error receiving data: ")
      + zmq_strerror(errno));
  }
  bytes_received = unsigned(nbytes);
  nbytes_pulled_ += nbytes;
  return true;
}

bool ZMQPuller::pull_assert_size(void* data, unsigned buffer_size,
  bool dont_wait)
{
  unsigned bytes_received = 0;
  bool good = pull(data, buffer_size, bytes_received, dont_wait);
  if(good and bytes_received != buffer_size)
    throw std::runtime_error(std::string("ZMQPuller: received unexpected "
      "number of bytes: ") + std::to_string(bytes_received) + " != "
      + std::to_string(buffer_size));
  return good;
}

bool ZMQPuller::wait_for_data(long timeout_ms)
{
  zmq_pollitem_t item = this->pollitem();
  int rc = zmq_poll(&item, 1, timeout_ms);
  return rc>0;
}

int ZMQPuller::wait_for_data_multi_source(ZMQPuller* puller2, long timeout_ms)
{
  zmq_pollitem_t items[2];
  items[0] = this->pollitem();
  if(puller2) {
    items[1] = puller2->pollitem();
  }
  int rc = zmq_poll(items, (puller2)?2:1, timeout_ms);
  if(rc<=0)return rc;
  if(items[0].revents & ZMQ_POLLIN)return 1;
  if((puller2) and (items[1].revents & ZMQ_POLLIN))return 2;
  throw std::logic_error("Impossible case in ZMQPuller::wait_for_data_multi_source");
  return -1;
}

ZMQInprocPushPull::ZMQInprocPushPull(void* extern_ctx, unsigned address_index,
    ZMQProtocol protocol, unsigned buffer_size):
  buffer_size_(buffer_size),
  zmq_ctx_(extern_ctx), address_index_(address_index), protocol_(protocol)
{
  // nothing to see here
}

ZMQInprocPushPull::
ZMQInprocPushPull(unsigned buffer_size, ZMQInprocPushPull* shared_ctx, ZMQProtocol protocol):
  buffer_size_(buffer_size),
  my_zmq_ctx_((shared_ctx!=nullptr and shared_ctx->my_zmq_ctx_!=nullptr) ? nullptr : zmq_ctx_new()),
  zmq_ctx_((shared_ctx!=nullptr and shared_ctx->my_zmq_ctx_!=nullptr) ? shared_ctx->my_zmq_ctx_ : my_zmq_ctx_),
  address_index_((shared_ctx!=nullptr and shared_ctx->my_zmq_ctx_!=nullptr) ? shared_ctx->zmq_ctx_address_.fetch_add(1) : zmq_ctx_address_.fetch_add(1)),
  protocol_(protocol)
{
  if(my_zmq_ctx_)zmq_ctx_set(my_zmq_ctx_, ZMQ_IO_THREADS, 0); // inproc only
}

ZMQInprocPushPull::~ZMQInprocPushPull()
{
  if(my_zmq_ctx_)zmq_ctx_destroy(my_zmq_ctx_);
}

ZMQPuller* ZMQInprocPushPull::new_puller(ZMQBindOrConnect bind_or_connect)
{
  return new ZMQPuller(zmq_ctx(), address(), buffer_size_, bind_or_connect, protocol_);
}

ZMQPusher* ZMQInprocPushPull::new_pusher(ZMQBindOrConnect bind_or_connect)
{
  return new ZMQPusher(zmq_ctx(), address(), buffer_size_, bind_or_connect, protocol_);
}

std::string ZMQInprocPushPull::address()
{
  return std::string("inproc://#")+std::to_string(address_index_);
}
