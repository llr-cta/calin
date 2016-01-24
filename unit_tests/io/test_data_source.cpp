/*

   calin/unit_tests/io/test_data_source.cpp -- Stephen Fegan -- 2016-01-24

   Unit tests for data source classes

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

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <gtest/gtest.h>

#include <io/data_source.hpp>
#include <unittest.pb.h>

using namespace calin::ix::unittest;
using namespace calin::io::data_source;

using UTSSM_RADS = RandomAccessDataSource<UnitTestSimpleSubMessage>;

class UnitTestIntegerDataSource: public UTSSM_RADS
{
public:
  UnitTestIntegerDataSource(unsigned count, int32_t start = 0):
    UTSSM_RADS(), count_(count), start_(start) { }
  ~UnitTestIntegerDataSource() { }
  UnitTestSimpleSubMessage* get_next() override {
    if(index_>=count_)return nullptr;
    UnitTestSimpleSubMessage* m = new UnitTestSimpleSubMessage;
    m->set_ssm_i32(start_ + int32_t(index_));
    ++index_;
    return m;
  }
  uint64_t size() override { return count_; }
  void set_next_index(uint64_t index) override {
    index_ = std::min(count_, unsigned(index)); }
private:
  unsigned index_ = 0;
  unsigned count_ = 0;
  int32_t start_ = 0;
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
    auto* m = src.get_next();
    ASSERT_NE(m, nullptr);
    EXPECT_EQ(m->ssm_i32(), int32_t(i));
    delete m;
  }
  auto* m = src.get_next();
  ASSERT_EQ(m, nullptr);
}

TEST(TestIntegerDataSource, RandomAccess) {
  unsigned N = 1000;
  UnitTestIntegerDataSource src(N,0);
  for(unsigned i=0;i<10*N;i++)
  {
    unsigned index = random()%((11*N)/10);
    src.set_next_index(index);
    auto* m = src.get_next();
    if(index>=N) {
      ASSERT_EQ(m, nullptr);
    } else {
      ASSERT_NE(m, nullptr);
      EXPECT_EQ(m->ssm_i32(), int32_t(index));
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
      auto* m = src->get_next();
      ASSERT_NE(m, nullptr);
      EXPECT_EQ(m->ssm_i32(), int32_t(n));
      delete m;
    }
    auto* m = src->get_next();
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
  BasicChaninedRandomAccessDataSource<UTSSM_RADS> src(&opener, false);
  for(unsigned i=0;i<N*M;i++)
  {
    auto* m = src.get_next();
    ASSERT_NE(m, nullptr);
    EXPECT_EQ(m->ssm_i32(), int32_t(i));
    delete m;
  }
  ASSERT_EQ(src.get_next(), nullptr);
  ASSERT_EQ(src.get_next(), nullptr);
}

TEST(TestChainedRandomAccessDataSource, RandomAccess) {
  unsigned N = 1000;
  unsigned M = 100;
  for(unsigned j=0; j<10; j++)
  {
    UnitTestDataSourceOpener opener(M,N);
    BasicChaninedRandomAccessDataSource<UTSSM_RADS> src(&opener, false);

    for(unsigned i=0;i<N*M;i++)
    {
      unsigned index = random()%((11*N*M)/10);
      src.set_next_index(index);
      auto* m = src.get_next();
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
    for(unsigned i=0;i<N;i++)file_out.put_next(src.get_next());
  }

  ProtobufFileDataSource<UnitTestSimpleSubMessage>
    file_in("unittest.proto_raw");
  for(unsigned i=0;i<N;i++)
  {
    auto* m = file_in.get_next();
    ASSERT_NE(m, nullptr);
    EXPECT_EQ(m->ssm_i32(), int32_t(i));
    delete m;
  }
  ASSERT_EQ(file_in.get_next(), nullptr);
  ASSERT_EQ(file_in.get_next(), nullptr);
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
