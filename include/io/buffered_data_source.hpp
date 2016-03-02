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

  BufferedDataSource(zmq_inproc::ZMQPuller* puller,
    std::atomic<bool>& reader_stopping): DataSource<T>(), puller_(puller),
    reader_stopping_(reader_stopping)
  {
    // nothing to see here
  }

  virtual ~BufferedDataSource()
  {
    delete puller_;
  }

  T* get_next() override
  {
    void* ptr = nullptr;
    puller_->pull_assert_size(&ptr, sizeof(ptr));
//    puller_->pull_assert_size(&ptr, sizeof(ptr), reader_stopping_);
//    log::LOG(log::INFO) << "Source: " << safe_downcast<T>(ptr);
    return safe_downcast<T>(ptr);
  }

private:
  zmq_inproc::ZMQPuller* puller_ = nullptr;
  std::atomic<bool>& reader_stopping_;
};

template<typename T> class MultiThreadDataSourceBuffer
{
public:
  CALIN_TYPEALIAS(data_type, typename DataSource<T>::data_type);

  MultiThreadDataSourceBuffer(const MultiThreadDataSourceBuffer&) = delete;
  MultiThreadDataSourceBuffer& operator=(const MultiThreadDataSourceBuffer&) = delete;

  MultiThreadDataSourceBuffer(DataSource<T>* source, bool adopt_source = false):
    source_(source), adopt_source_(adopt_source),
    zmq_(new zmq_inproc::ZMQInprocPushPull),
    reader_thread_(
      new std::thread(&MultiThreadDataSourceBuffer::reader_loop, this))
  {
    // Wait until BIND happens
    while(!reader_active_){ CALIN_SPINWAIT(); }
  }

  ~MultiThreadDataSourceBuffer()
  {
    stop_buffering_ = true;
    delete zmq_;
    reader_thread_->join();
    delete reader_thread_;
    if(adopt_source_)delete source_;
  }

  BufferedDataSource<T>* new_data_source(unsigned buffer_size) {
    return new BufferedDataSource<T>(zmq_->new_puller(), stop_buffering_); }

  void stop_buffering() { stop_buffering_ = true; }

private:
  void reader_loop()
  {
    std::unique_ptr<zmq_inproc::ZMQPusher> pusher { zmq_->new_pusher() };
    reader_active_ = true;
    while(reader_active_)
    {
      // The reader stays active until the ZMQ context is deleted in the
      // desctructor. However it stops getting data from the supplied source
      // when requested by the user.
      T* p = nullptr;
      if(!stop_buffering_)p = source_->get_next();
      if(!pusher->push(&p, sizeof(p)))
      {
        // only occurs if context closed
        delete p;
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
