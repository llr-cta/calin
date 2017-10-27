/*

   calin/unit_tests/io/test_data_source.cpp -- Stephen Fegan -- 2016-01-24

   Unit tests for data source classes

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

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <gtest/gtest.h>

#include <io/log.hpp>
#include <io/data_source.hpp>
#include <io/buffered_data_source.hpp>
#include <io/chained_data_source.hpp>
#include <unittest.pb.h>

using namespace calin::io::log;
using namespace calin::ix::unittest;
using namespace calin::io::data_source;

using UTSSM_RADS = RandomAccessDataSource<UnitTestSimpleSubMessage>;

class UnitTestIntegerDataSource: public UTSSM_RADS
{
public:
  UnitTestIntegerDataSource(unsigned count, int32_t start = 0,
    unsigned delay = 0):
    UTSSM_RADS(), count_(count), start_(start), delay_(delay) { }
  ~UnitTestIntegerDataSource() { }
  UnitTestSimpleSubMessage* get_next(uint64_t& seq_index_out,
      google::protobuf::Arena** arena = nullptr) override {
    if(index_>=count_)return nullptr;
    UnitTestSimpleSubMessage* m = new UnitTestSimpleSubMessage;
    m->set_ssm_i32(start_ + int32_t(index_));
    seq_index_out = index_;
    ++index_;
    if(delay_)std::this_thread::sleep_for(std::chrono::microseconds(delay_));
    return m;
  }
  uint64_t size() override { return count_; }
  void set_next_index(uint64_t index) override {
    index_ = std::min(count_, unsigned(index)); }
private:
  unsigned index_ = 0;
  unsigned count_ = 0;
  int32_t start_ = 0;
  unsigned delay_ = 0;
};

class UnitTestDataSourceOpener: public DataSourceOpener<UTSSM_RADS>
{
public:
  using DataSourceOpener<UTSSM_RADS>::data_source_type;
  UnitTestDataSourceOpener(unsigned nsource, unsigned count = 10000):
    DataSourceOpener<UTSSM_RADS>(), nsource_(nsource), count_(count) { }
  ~UnitTestDataSourceOpener() { }
  unsigned num_sources() override {
    return nsource_;
  }
  data_source_type* open(unsigned isource) override {
    if(isource>=nsource_)return nullptr;
    return new UnitTestIntegerDataSource(count_, count_*isource);
  }
private:
  unsigned nsource_;
  unsigned count_;
};

TEST(TestIntegerDataSource, Sequential) {
  unsigned N = 1000;
  UnitTestIntegerDataSource src(N,0);
  for(unsigned i=0;i<N;i++)
  {
    uint64_t seq_index = 0;
    auto* m = src.get_next(seq_index);
    ASSERT_NE(m, nullptr);
    EXPECT_EQ(m->ssm_i32(), int32_t(i));
    EXPECT_EQ(seq_index, uint64_t(i));
    delete m;
  }
  uint64_t seq_index = 0;

  auto* m = src.get_next(seq_index);
  ASSERT_EQ(m, nullptr);
}

TEST(TestIntegerDataSource, RandomAccess) {
  unsigned N = 1000;
  UnitTestIntegerDataSource src(N,0);
  for(unsigned i=0;i<10*N;i++)
  {
    unsigned index = random()%((11*N)/10);
    src.set_next_index(index);
    uint64_t seq_index = 0;
    auto* m = src.get_next(seq_index);
    if(index>=N) {
      ASSERT_EQ(m, nullptr);
    } else {
      ASSERT_NE(m, nullptr);
      EXPECT_EQ(m->ssm_i32(), int32_t(index));
      EXPECT_EQ(seq_index, uint64_t(index));
      delete m;
    }
  }
}

TEST(TestDataSourceOpener, Sequential) {
  unsigned N = 1000;
  unsigned M = 100;
  UnitTestDataSourceOpener opener(M,N);
  unsigned n=0;
  for(unsigned j=0;j<M;j++)
  {
    auto* src = opener.open(j);
    for(unsigned i=0;i<N;++i,++n)
    {
      uint64_t seq_index = 0;
      auto* m = src->get_next(seq_index);
      ASSERT_NE(m, nullptr);
      EXPECT_EQ(m->ssm_i32(), int32_t(n));
      delete m;
    }
    uint64_t seq_index = 0;
    auto* m = src->get_next(seq_index);
    ASSERT_EQ(m, nullptr);
    delete src;
  }
  ASSERT_EQ(opener.open(M), nullptr);
  ASSERT_EQ(opener.open(M+1), nullptr);
}

TEST(TestChainedRandomAccessDataSource, Sequental) {
  unsigned N = 1000;
  unsigned M = 100;
  UnitTestDataSourceOpener opener(M,N);
  BasicChainedRandomAccessDataSource<UTSSM_RADS> src(&opener, false);
  for(unsigned i=0;i<N*M;i++)
  {
    uint64_t seq_index = 0;
    auto* m = src.get_next(seq_index);
    ASSERT_NE(m, nullptr);
    EXPECT_EQ(m->ssm_i32(), int32_t(i));
    delete m;
  }
  uint64_t seq_index = 0;
  ASSERT_EQ(src.get_next(seq_index), nullptr);
  ASSERT_EQ(src.get_next(seq_index), nullptr);
}

TEST(TestChainedRandomAccessDataSource, RandomAccess) {
  unsigned N = 1000;
  unsigned M = 100;
  for(unsigned j=0; j<10; j++)
  {
    UnitTestDataSourceOpener opener(M,N);
    BasicChainedRandomAccessDataSource<UTSSM_RADS> src(&opener, false);

    for(unsigned i=0;i<N*M;i++)
    {
      unsigned index = random()%((11*N*M)/10);
      src.set_next_index(index);
      uint64_t seq_index = 0;
      auto* m = src.get_next(seq_index);
      if(index>=N*M) {
        ASSERT_EQ(m, nullptr);
      } else {
        ASSERT_NE(m, nullptr);
        EXPECT_EQ(m->ssm_i32(), int32_t(index));
        delete m;
      }
    }
  }
}

TEST(TestProtobufFile, WriteAndRead) {
  unsigned N = 1000;
  if(1)
  {
    UnitTestIntegerDataSource src(N,0);
    ProtobufFileDataSink<UnitTestSimpleSubMessage>
      file_out("unittest.proto_raw");
    for(unsigned i=0;i<N;i++) {
      uint64_t seq_index = 0;
      auto* m = src.get_next(seq_index);
      file_out.put_next(m, seq_index, nullptr, true);
    }
  }

  ProtobufFileDataSource<UnitTestSimpleSubMessage>
    file_in("unittest.proto_raw");
  for(unsigned i=0;i<N;i++)
  {
    uint64_t seq_index = 0;
    auto* m = file_in.get_next(seq_index);
    ASSERT_NE(m, nullptr);
    EXPECT_EQ(m->ssm_i32(), int32_t(i));
    delete m;
  }
  uint64_t seq_index = 0;
  ASSERT_EQ(file_in.get_next(seq_index), nullptr);
  ASSERT_EQ(file_in.get_next(seq_index), nullptr);
}

TEST(TestBufferedIntegerDataSource, Sequential) {
  // Test that data packets all arrive in correct order
  unsigned N = 1000;
  UnitTestIntegerDataSource src(N,0);
  MultiThreadDataSourceBuffer<UnitTestSimpleSubMessage> buffer(&src, N/10);
  auto* bsrc = buffer.new_data_source();

  for(unsigned i=0;i<N;i++)
  {
    uint64_t seq_index = 0;
    auto* m = bsrc->get_next(seq_index);
    ASSERT_NE(m, nullptr);
    EXPECT_EQ(m->ssm_i32(), int32_t(i));
    delete m;
  }
  uint64_t seq_index = 0;
  auto* m = bsrc->get_next(seq_index);
  ASSERT_EQ(m, nullptr);
  delete bsrc;
}

TEST(TestBufferedIntegerDataSource, SequentialWithStop) {
  // Test that no data packets are lost when stop is called
  unsigned N = 1000;
  UnitTestIntegerDataSource src(N,0);
  MultiThreadDataSourceBuffer<UnitTestSimpleSubMessage> buffer(&src,10);
  auto* bsrc = buffer.new_data_source();

  unsigned i=0;
  for(i=0;i<N/2;i++)
  {
    uint64_t seq_index = 0;
    auto* m = bsrc->get_next(seq_index);
    ASSERT_NE(m, nullptr);
    EXPECT_EQ(m->ssm_i32(), int32_t(i));
    delete m;
  }
  buffer.stop_buffering();
  uint64_t seq_index = 0;
  while(auto* m = bsrc->get_next(seq_index))
  {
    EXPECT_EQ(m->ssm_i32(), int32_t(i));
    delete m;
    i++;
  }
  delete bsrc;
  auto* m = src.get_next(seq_index);
  ASSERT_NE(m, nullptr);
  EXPECT_EQ(m->ssm_i32(), int32_t(i));
  delete m;
}

TEST(TestBufferedIntegerDataSource, MultiThreaded) {
  // Test that all data packets arrive, and that they come in order
  unsigned N = 10000;
  unsigned delay = 100;
  UnitTestIntegerDataSource src(N,0);
  MultiThreadDataSourceBuffer<UnitTestSimpleSubMessage> buffer(&src,10);

  std::vector<unsigned> ids(N);
  std::vector<std::thread> threads;

  for(unsigned ithread=0;ithread<10;ithread++)
    threads.emplace_back([&buffer,&ids,ithread,delay](){
      auto* bsrc = buffer.new_data_source();
      int last_index = -1;
      uint64_t seq_index = 0;
      while(auto* m = bsrc->get_next(seq_index))
      {
        EXPECT_GT(m->ssm_i32(), last_index);
        last_index = m->ssm_i32();
        ids[last_index] = ithread+1;
        delete m;
        if(delay)std::this_thread::sleep_for(std::chrono::microseconds(delay));
      }
      delete bsrc;
    });
  for(auto& i : threads)i.join();

  for(unsigned i=0;i<N;i++)EXPECT_GT(ids[i],0U);
}

TEST(TestBufferedIntegerDataSource, MultiThreadedWithStop) {
  // Test that no packets are lost in multi threaded scenario
  unsigned N = 10000;
  unsigned delay = 100;
  UnitTestIntegerDataSource src(N,0);
  MultiThreadDataSourceBuffer<UnitTestSimpleSubMessage> buffer(&src,10);

  std::vector<unsigned> ids(N);
  std::vector<std::thread> threads;

  for(unsigned ithread=0;ithread<10;ithread++)
    threads.emplace_back([&buffer,&ids,ithread,delay](){
      auto* bsrc = buffer.new_data_source();
      int last_index = -1;
      uint64_t seq_index = 0;
      while(auto* m = bsrc->get_next(seq_index))
      {
        EXPECT_GT(m->ssm_i32(), last_index);
        last_index = m->ssm_i32();
        ids[last_index] = ithread+1;
        delete m;
        if(delay)std::this_thread::sleep_for(std::chrono::microseconds(delay));
      }
      delete bsrc;
    });
  buffer.stop_buffering();
  for(auto& i : threads)i.join();

  uint64_t seq_index = 0;
  auto* m = src.get_next(seq_index);
  unsigned i = 0;
  if(m != nullptr)
    for(;i<unsigned(m->ssm_i32());i++)EXPECT_GT(ids[i],0U) << "With: i=" << i;
  for(;i<N;i++)EXPECT_EQ(ids[i],0U) << "With: i=" << i;
  delete m;
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
