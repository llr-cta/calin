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

#if GOOGLE_PROTOBUF_VERSION < 3006000
typedef google::protobuf::internal::true_type true_type;
typedef google::protobuf::internal::false_type false_type;
#else
typedef std::true_type true_type;
typedef std::false_type false_type;
#endif

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
  bool &warning_sent, const true_type&)
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
  bool &warning_sent, const false_type&)
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
  const true_type&)
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
  const false_type&)
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
    if(not pack_payload<T>(payload, data, seq_index, arena, adopt_data,
      google::protobuf::Arena::is_arena_constructable<T>()))return false;
    return pusher_->push(&payload, sizeof(payload));
  }

private:
  zmq_inproc::ZMQPusher* pusher_ = nullptr;
};

// Launch thread to read events from supplied DataSource, buffer them in 0MQ
// Push/Pull socket for use by other threads
template<typename T> class UnidirectionalBufferedDataSourcePump:
  public DataSourceFactory<T>
{
public:
  CALIN_TYPEALIAS(data_type, typename DataSource<T>::data_type);

  UnidirectionalBufferedDataSourcePump(const UnidirectionalBufferedDataSourcePump&) = delete;
  UnidirectionalBufferedDataSourcePump& operator=(const UnidirectionalBufferedDataSourcePump&) = delete;

  UnidirectionalBufferedDataSourcePump(DataSource<T>* source, unsigned buffer_size = 100,
      bool adopt_source = false):
    DataSourceFactory<T>(),
    source_(source), adopt_source_(adopt_source),
    zmq_(new zmq_inproc::ZMQInprocPushPull(buffer_size)),
    reader_thread_(
      new std::thread(&UnidirectionalBufferedDataSourcePump::reader_loop, this))
  {
    // Wait until BIND happens
    while(!reader_active_){ CALIN_SPINWAIT(); }
  }

  virtual ~UnidirectionalBufferedDataSourcePump()
  {
    stop_buffering_ = true;
    // This closes the ZMQ context prompting all readers and the writer to
    // unblock. It must happen before the join to the reader thread.
    delete zmq_;
    reader_thread_->join();
    delete reader_thread_;
    if(adopt_source_)delete source_;
  }

  BufferedDataSource<T>* new_data_source() override {
    return new BufferedDataSource<T>(zmq_->new_puller()); }

  void stop_buffering() { stop_buffering_ = true; }

private:
  void reader_loop()
  {
    try {
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
    } catch(const std::exception& x) {
      util::log::LOG(util::log::FATAL) << x.what();
      throw;
    }
  }

  DataSource<T>* source_ { nullptr };
  bool adopt_source_ { false };
  std::atomic<bool> stop_buffering_ { false };
  std::atomic<bool> reader_active_ { false };
  zmq_inproc::ZMQInprocPushPull* zmq_ = nullptr;
  std::thread* reader_thread_ = nullptr;
};

template<typename T> class BidirectionalBufferedDataSourcePump:
  public DataSourceFactory<T>, public DataSinkFactory<T>
{
public:
  CALIN_TYPEALIAS(data_type, typename DataSource<T>::data_type);

  BidirectionalBufferedDataSourcePump(const BidirectionalBufferedDataSourcePump&) = delete;
  BidirectionalBufferedDataSourcePump& operator=(const BidirectionalBufferedDataSourcePump&) = delete;

  BidirectionalBufferedDataSourcePump(DataSource<T>* source, DataSink<T>* sink,
      unsigned buffer_size = 100, bool sink_unsent_data = false,
      bool adopt_source = false, bool adopt_sink = false):
    DataSourceFactory<T>(), DataSinkFactory<T>(),
    source_(source), adopt_source_(adopt_source),
    sink_(sink), adopt_sink_(adopt_sink),
    sink_unsent_data_(sink_unsent_data),
    zmq_dn_(new zmq_inproc::ZMQInprocPushPull(buffer_size)),
    zmq_up_(new zmq_inproc::ZMQInprocPushPull(buffer_size, zmq_dn_)),
    ventilator_thread_(
      new std::thread(&BidirectionalBufferedDataSourcePump::ventilator_loop, this))
  {
    // Wait until BIND happens
    while(!ventilator_active_){ CALIN_SPINWAIT(); }
  }

  virtual ~BidirectionalBufferedDataSourcePump()
  {
    stop_ventilator_ = true;
    // This closes the ZMQ context prompting all readers and the writer to
    // unblock. It must happen before the join to the reader thread.
    delete zmq_up_;
    delete zmq_dn_;
    ventilator_thread_->join();
    delete ventilator_thread_;
    if(adopt_source_)delete source_;
    if(adopt_sink_)delete sink_;
  }

  BufferedDataSource<T>* new_data_source() override {
    return new BufferedDataSource<T>(zmq_dn_->new_puller(zmq_inproc::ZMQBindOrConnect::CONNECT)); }

  BufferedDataSink<T>* new_data_sink() override {
    return new BufferedDataSink<T>(zmq_up_->new_pusher(zmq_inproc::ZMQBindOrConnect::CONNECT)); }

  void stop_buffering() { stop_ventilator_ = true; }
  void stop_ventilator() { stop_ventilator_ = true; }

private:
  void ventilator_loop()
  {
    try {
      std::unique_ptr<zmq_inproc::ZMQPusher> pusher { zmq_dn_->new_pusher(zmq_inproc::ZMQBindOrConnect::BIND) };
      std::unique_ptr<zmq_inproc::ZMQPuller> puller { zmq_up_->new_puller(zmq_inproc::ZMQBindOrConnect::BIND) };
      ventilator_active_ = true;
      Payload<T> payload_push;
      payload_push.ptr = source_->get_next(payload_push.seq_index, &payload_push.arena);
      while(ventilator_active_)
      {
        // The reader stays active until the ZMQ context is closed in the
        // desctructor. However it stops getting data from the supplied source
        // when requested by the user. After this it just delivers "nullptr".
        // This simplifies the implementation of the BufferedDataSource which
        // otherwise has a problem either potentially getting blocked or missing
        // a data item.
        zmq_pollitem_t pollitems[] = { puller->pollitem(), pusher->pollitem() };
        if(zmq_poll(pollitems, 2, -1) < 0) {
          // Only occurs when context is closed
          ventilator_active_ = false;
        } else {
          if(pollitems[0].revents) {
            // Prioritize sata from puller - keep getting some while its available
            Payload<T> payload_pull;
            while(puller->pull_assert_size(&payload_pull, sizeof(payload_pull),
                /* dont_wait= */ true)) {
              sink_->put_next(payload_pull.ptr, payload_pull.seq_index,
                payload_pull.arena, /* adopt_data = */ true);
              payload_pull.ptr = nullptr;
              payload_pull.arena = nullptr;
            }
          }
          if(pollitems[1].revents) {
            if(pusher->push(&payload_push, sizeof(payload_push))) {
              payload_push.ptr = nullptr;
              payload_push.arena = nullptr;
              if(!stop_ventilator_)
                payload_push.ptr = source_->get_next(payload_push.seq_index, &payload_push.arena);
            } else {
              ventilator_active_ = false;
            }
          }
        }
      }
      if(payload_push.ptr) {
        if(sink_unsent_data_)
          sink_->put_next(payload_push.ptr, payload_push.seq_index,
            payload_push.arena, /* adopt_data = */ true);
        else if(payload_push.arena)delete payload_push  .arena;
        else delete payload_push.ptr;
      }
    } catch(const std::exception& x) {
      util::log::LOG(util::log::FATAL) << x.what();
      throw;
    }
  }

  DataSource<T>* source_ { nullptr };
  bool adopt_source_ { false };
  DataSink<T>* sink_ { nullptr };
  bool adopt_sink_ { false };
  bool sink_unsent_data_ = false;
  std::atomic<bool> stop_ventilator_ { false };
  std::atomic<bool> ventilator_active_ { false };
  zmq_inproc::ZMQInprocPushPull* zmq_dn_ = nullptr;
  zmq_inproc::ZMQInprocPushPull* zmq_up_ = nullptr;
  std::thread* ventilator_thread_ = nullptr;
};

} } } // namespace calin::io::data_source
