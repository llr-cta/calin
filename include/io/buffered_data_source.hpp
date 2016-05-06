/*

   calin/io/buffered_data_source.hpp -- Stephen Fegan -- 2016-02-29

   A multi-threaded data source which buffers data from a parent source
   making it available

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

namespace calin { namespace io { namespace data_source {

template<typename T> struct Payload
{
public:
  Payload(T* ptr_ = nullptr, google::protobuf::Arena* arena_ = nullptr ,
    uint64_t seq_index_ = 0):
    ptr(ptr_), arena(arena_), seq_index(seq_index_) { }
  T* ptr = nullptr;
  google::protobuf::Arena* arena = nullptr;
  uint64_t seq_index = 0;
};

template<typename T> inline T* safe_downcast(void* p)
{
  // Not really safe at all!
  return static_cast<T*>(p);
}

template<typename T> class BufferedDataSource: public DataSource<T>
{
public:
  CALIN_TYPEALIAS(data_type, typename DataSource<T>::data_type);

  BufferedDataSource(const BufferedDataSource&) = delete;
  BufferedDataSource& operator=(const BufferedDataSource&) = delete;

  BufferedDataSource(zmq_inproc::ZMQPuller* puller):
    DataSource<T>(), puller_(puller)
  {
    // nothing to see here
  }

  virtual ~BufferedDataSource()
  {
    delete puller_;
  }

  T* get_next(uint64_t& seq_index_out,
    google::protobuf::Arena** arena = nullptr) override
  {
    Payload<T> payload;
    puller_->pull_assert_size(&payload, sizeof(payload));
    if(payload.ptr == nullptr) {
      assert(payload.arena == nullptr);
      return nullptr;
    }
    seq_index_out = payload.seq_index;
    if(arena == nullptr)
    {
      if(payload.arena == nullptr)return payload.ptr;
      if(!warning_sent_) {
        io::log::LOG(io::log::WARNING) << "BufferedDataSource: "
          "arena not accepted, performing expensive copy";
        warning_sent_ = true;
      }
      T* data = new T(*payload.ptr);
      delete payload.arena;
      return data;
    }
    else if(*arena == nullptr)
    {
      *arena = payload.arena;
      return payload.ptr;
    }
    else
    {
      if(!warning_sent_) {
        io::log::LOG(io::log::WARNING) << "BufferedDataSource: "
          "pre-assigned arena, performing expensive copy";
        warning_sent_ = true;
      }
      T* data = google::protobuf::Arena::CreateMessage<T>(*arena);
      data->CopyFrom(*payload.ptr);
      if(payload.arena)delete payload.arena;
      else delete payload.ptr;
      return data;
    }
  }

private:
  bool warning_sent_ = false;
  zmq_inproc::ZMQPuller* puller_ = nullptr;
};

template<typename T> class MultiThreadDataSourceBuffer
{
public:
  CALIN_TYPEALIAS(data_type, typename DataSource<T>::data_type);

  MultiThreadDataSourceBuffer(const MultiThreadDataSourceBuffer&) = delete;
  MultiThreadDataSourceBuffer& operator=(const MultiThreadDataSourceBuffer&) = delete;

  MultiThreadDataSourceBuffer(DataSource<T>* source, unsigned buffer_size = 100,
      bool adopt_source = false):
    source_(source), adopt_source_(adopt_source),
    zmq_(new zmq_inproc::ZMQInprocPushPull(buffer_size)),
    reader_thread_(
      new std::thread(&MultiThreadDataSourceBuffer::reader_loop, this))
  {
    // Wait until BIND happens
    while(!reader_active_){ CALIN_SPINWAIT(); }
  }

  ~MultiThreadDataSourceBuffer()
  {
    stop_buffering_ = true;
    // This closes the ZMQ context prompting all readers and the writer to
    // unblock. It must happen before the join to the reader thread.
    delete zmq_;
    reader_thread_->join();
    delete reader_thread_;
    if(adopt_source_)delete source_;
  }

  BufferedDataSource<T>* new_data_source() {
    return new BufferedDataSource<T>(zmq_->new_puller()); }

  void stop_buffering() { stop_buffering_ = true; }

private:
  void reader_loop()
  {
    std::unique_ptr<zmq_inproc::ZMQPusher> pusher { zmq_->new_pusher() };
    reader_active_ = true;
    while(reader_active_)
    {
      // The reader stays active until the ZMQ context is closed in the
      // desctructor. However it stops getting data from the supplied source
      // when requested by the user. After this it just delivers "nullptr".
      // This simplifies the implementation of the BufferedDataSource which
      // otherwise has a problem either potentially getting blocked or missing
      // a data item.
      Payload<T> payload;
      if(!stop_buffering_)
        payload.ptr = source_->get_next(payload.seq_index, &payload.arena);
      if(!pusher->push(&payload, sizeof(payload)))
      {
        // Only occurs when context is closed
        if(payload.arena)delete payload.arena;
        else delete payload.ptr;
        reader_active_ = false;
      }
    }
  }

  DataSource<T>* source_ { nullptr };
  bool adopt_source_ { false };
  std::atomic<bool> stop_buffering_ { false };
  std::atomic<bool> reader_active_ { false };
  zmq_inproc::ZMQInprocPushPull* zmq_ = nullptr;
  std::thread* reader_thread_ = nullptr;
};

} } } // namespace calin::io::data_source
