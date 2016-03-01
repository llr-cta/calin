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
    std::atomic<bool>& reader_loop_finished): DataSource<T>(), puller_(puller),
    reader_loop_finished_(reader_loop_finished)
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
    puller_->pull_assert_size(&ptr, sizeof(ptr), reader_loop_finished_);
//    log::LOG(log::INFO) << "Source: " << safe_downcast<T>(ptr);
    return safe_downcast<T>(ptr);
  }

private:
  zmq_inproc::ZMQPuller* puller_ = nullptr;
  std::atomic<bool>& reader_loop_finished_;
};

template<typename T> class MultiThreadDataSourceBuffer
{
public:
  CALIN_TYPEALIAS(data_type, typename DataSource<T>::data_type);

  MultiThreadDataSourceBuffer(const MultiThreadDataSourceBuffer&) = delete;
  MultiThreadDataSourceBuffer& operator=(const MultiThreadDataSourceBuffer&) = delete;

  MultiThreadDataSourceBuffer(DataSource<T>* source, bool adopt_source = false):
    source_(source), adopt_source_(adopt_source),
    reader_thread_(&MultiThreadDataSourceBuffer::reader_loop, this)
  {
    // nothing to see here
  }

  ~MultiThreadDataSourceBuffer()
  {
    stop_buffering_ = true;
    std::unique_ptr<zmq_inproc::ZMQPuller> puller { zmq_.new_puller() };
    void* ptr = nullptr;
    while(puller->pull_assert_size(&ptr, sizeof(ptr), reader_loop_finished_))
      delete safe_downcast<T>(ptr);
    reader_thread_.join();
    if(adopt_source_)delete source_;
  }

  BufferedDataSource<T>* new_data_source(unsigned buffer_size) {
    return new BufferedDataSource<T>(zmq_.new_puller(), reader_loop_finished_); }

  void stop_buffering() { stop_buffering_ = true; }

private:
  void reader_loop()
  {
    std::unique_ptr<zmq_inproc::ZMQPusher> pusher { zmq_.new_pusher() };
    while(!stop_buffering_)
    {
      T* p = source_->get_next();
      pusher->push(&p, sizeof(p));
//      log::LOG(log::INFO) << "Reader: " << p;
    }
    reader_loop_finished_ = true;
  }

  DataSource<T>* source_ { nullptr };
  bool adopt_source_ { false };
  std::atomic<bool> stop_buffering_ { false };
  std::atomic<bool> reader_loop_finished_ { false };
  zmq_inproc::ZMQInprocPushPull zmq_;
  std::thread reader_thread_;
};

} } } // namespace calin::io::data_source
