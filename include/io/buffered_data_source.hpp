/*

   calin/io/buffered_data_source.hpp -- Stephen Fegan -- 2016-02-29

   A multi-threaded data source which buffers data from a parent source
   making it available

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <atomic>
#include <thread>
#include <cassert>

#include <util/spinlock.hpp>
#include <util/log.hpp>
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

template<typename T> inline T* unpack_payload(const Payload<T>& payload,
  uint64_t& seq_index_out, google::protobuf::Arena** arena,
  bool &warning_sent, const google::protobuf::internal::true_type&)
{
  if(payload.ptr == nullptr)return nullptr;
  seq_index_out = payload.seq_index;
  if(arena == nullptr)
  {
    if(payload.arena == nullptr)return payload.ptr;
    if(!warning_sent) {
      util::log::LOG(util::log::WARNING) << "unpack_payload: "
        "arena not accepted, performing expensive copy";
      warning_sent = true;
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
    if(!warning_sent) {
      util::log::LOG(util::log::WARNING) << "unpack_payload: "
        "pre-assigned arena, performing expensive copy";
      warning_sent = true;
    }
    T* data = google::protobuf::Arena::CreateMessage<T>(*arena);
    data->CopyFrom(*payload.ptr);
    if(payload.arena)delete payload.arena;
    else delete payload.ptr;
    return data;
  }
}

template<typename T> inline T* unpack_payload(const Payload<T>& payload,
  uint64_t& seq_index_out, google::protobuf::Arena** arena,
  bool &warning_sent, const google::protobuf::internal::false_type&)
{
  if(payload.arena != nullptr)
    throw std::runtime_error("unpack_payload: logic error, non-null payload arena for type that does not support arena");
  if(arena != nullptr and *arena != nullptr)
    throw std::runtime_error("unpack_payload: logic error, non-null user arena for type that does not support arena");
  seq_index_out = payload.seq_index;
  return payload.ptr;
}

template<typename T> bool pack_payload(Payload<T>& payload,
  T* ptr, uint64_t seq_index, google::protobuf::Arena* arena, bool adopt_data,
  const google::protobuf::internal::true_type&)
{
  payload.seq_index = seq_index;
  if(adopt_data) {
    payload.arena = arena;
    payload.ptr = ptr;
  } else if(ptr and arena) {
    payload.arena = new google::protobuf::Arena;
    T* data = google::protobuf::Arena::CreateMessage<T>(*payload.arena);
    data->CopyFrom(*payload.ptr);
    payload.ptr = data;
  } else if(ptr) {
    payload.arena = nullptr;
    payload.ptr = new T(*ptr);
  } else {
    payload.arena = nullptr;
    payload.ptr = nullptr;
  }
  return true;
}

template<typename T> bool pack_payload(Payload<T>& payload,
  T* ptr, uint64_t seq_index, google::protobuf::Arena* arena, bool adopt_data,
  const google::protobuf::internal::false_type&)
{
  if(arena != nullptr)
    throw std::runtime_error("pack_payload: logic error, non-null user arena for type that does not support arena");
  payload.arena = nullptr;
  payload.seq_index = seq_index;
  payload.ptr = (ptr==nullptr or adopt_data) ? ptr : new T(*ptr);
  return true;
}

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
    return unpack_payload<T>(payload, seq_index_out, arena, warning_sent_,
      google::protobuf::Arena::is_arena_constructable<T>());
  }

private:
  bool warning_sent_ = false;
  zmq_inproc::ZMQPuller* puller_ = nullptr;
};

template<typename T> class BufferedDataSink: public DataSink<T>
{
public:
  CALIN_TYPEALIAS(data_type, typename DataSink<T>::data_type);

  BufferedDataSink(const BufferedDataSink&) = delete;
  BufferedDataSink& operator=(const BufferedDataSink&) = delete;

  BufferedDataSink(zmq_inproc::ZMQPusher* pusher):
    DataSink<T>(), pusher_(pusher)
  {
    // nothing to see here
  }

  virtual ~BufferedDataSink()
  {
    delete pusher_;
  }

  bool put_next(T* data, uint64_t seq_index,
    google::protobuf::Arena* arena = nullptr, bool adopt_data = false) override
  {
    Payload<T> payload;
    if(not pack_payload<T>(payload, data, seq_index, arena,
      google::protobuf::Arena::is_arena_constructable<T>()))return false;
    return pusher_->push(&payload, sizeof(payload));
  }

private:
  zmq_inproc::ZMQPusher* pusher_ = nullptr;
};

template<typename T> class UnidirectionalDataSourcePump
{
public:
  CALIN_TYPEALIAS(data_type, typename DataSource<T>::data_type);

  UnidirectionalDataSourcePump(const UnidirectionalDataSourcePump&) = delete;
  UnidirectionalDataSourcePump& operator=(const UnidirectionalDataSourcePump&) = delete;

  UnidirectionalDataSourcePump(DataSource<T>* source, unsigned buffer_size = 100,
      bool adopt_source = false):
    source_(source), adopt_source_(adopt_source),
    zmq_(new zmq_inproc::ZMQInprocPushPull(buffer_size)),
    reader_thread_(
      new std::thread(&UnidirectionalDataSourcePump::reader_loop, this))
  {
    // Wait until BIND happens
    while(!reader_active_){ CALIN_SPINWAIT(); }
  }

  ~UnidirectionalDataSourcePump()
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
