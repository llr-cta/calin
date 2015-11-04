/* 

   calin/unit_tests/math/test_histogram.cpp -- Stephen Fegan -- 2015-03-15

   Unit tests for histogram classes

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <gtest/gtest.h>
#include "math/histogram.hpp"

using namespace calin::math;
using namespace calin::math::histogram;
using namespace calin::math::accumulator;
using namespace calin::ix::math;

TEST(TestSimpleHist, InsertIntegersFrom0ToN) {
  unsigned N { 1000 };
  SimpleHist myhist {1.0};
  for(unsigned i=0;i<N;++i)
    for(unsigned j=i;j<N;++j)
      myhist.insert(j);  
  ASSERT_EQ(myhist.size(), N);
  for(unsigned i=0;i<N;++i)
    EXPECT_EQ(myhist.xval_center(i), i);  
  for(unsigned i=0;i<N;++i)
    EXPECT_EQ(myhist.weight(i), i+1);
}

TEST(TestSimpleHist, InsertIntegersFromNTo0) {
  unsigned N { 1000 };
  SimpleHist myhist {1.0};
  for(unsigned i=0;i<N;++i)
    for(unsigned j=i;j<N;++j)
      myhist.insert(N-(j+1));
  ASSERT_EQ(myhist.size(), N);
  for(unsigned i=0;i<N;++i)
    EXPECT_EQ(myhist.xval_center(i), i);  
  for(unsigned i=0;i<N;++i)
    EXPECT_EQ(myhist.weight(i), N-i);
}

TEST(TestSimpleHist, InsertIntegersFrom0ToNWithNegativeDX) {
  unsigned N { 1000 };
  SimpleHist myhist {-1.0};
  for(unsigned i=0;i<N;++i)
    for(unsigned j=i;j<N;++j)
      myhist.insert(j);  
  ASSERT_EQ(myhist.size(), N);
  for(unsigned i=0;i<N;++i)
    EXPECT_EQ(myhist.xval_center(i), N-i-1);  
  for(unsigned i=0;i<N;++i)
    EXPECT_EQ(myhist.weight(i), N-i);
}

TEST(TestSimpleHist, LimitedInsertIntegersFrom0ToN) {
  unsigned N { 1000 };
  SimpleHist myhist {1.0,99.5,899.5};
  for(unsigned i=0;i<N;++i)
    for(unsigned j=i;j<N;++j)
      myhist.insert(j);  
  ASSERT_EQ(myhist.size(), N-200);
  for(unsigned i=0;i<N-200;++i)
    EXPECT_EQ(myhist.xval_center(i), i+100);  
  for(unsigned i=0;i<N-200;++i)
    EXPECT_EQ(myhist.weight(i), i+101);
  EXPECT_TRUE(myhist.is_limited());
  EXPECT_EQ(myhist.xval_limit_lo(),99.5);
  EXPECT_EQ(myhist.xval_limit_hi(),899.5);
  EXPECT_EQ(myhist.weight_overflow_lo(),100*101/2);
  EXPECT_EQ(myhist.weight_overflow_hi(),100*101/2+900*100);
}


TEST(TestSimpleHist, RangeFor) {
  unsigned N { 10000 };
  SimpleHist myhist {1.0};
  for(unsigned i=0;i<N;++i)
    myhist.insert(i,i);
  unsigned i=0;
  for(auto ibin : myhist)EXPECT_EQ(ibin.weight(), i++);
  ASSERT_EQ(i,N);
}

template<typename FB, typename FE, typename FC>
void iterator_test(FB fbegin, FE fend, FC fival)
{
  unsigned N { 10000 };
  SimpleHist myhist {1.0};
  for(unsigned i=0;i<N;++i)
    myhist.insert(i,i);

  // Test post increment
  unsigned i=0;
  for(auto ibin = fbegin(myhist); ibin!=fend(myhist); ibin++)
    EXPECT_EQ(ibin->weight(), fival(i++,N));
  ASSERT_EQ(i,N);

  // Test pre increment
  i=0;
  for(auto ibin = fbegin(myhist); ibin!=fend(myhist); )
  {
    EXPECT_EQ(ibin->weight(), fival(i,N));
    ++ibin;
    ++i;
  }
  ASSERT_EQ(i,N);

  // Test pre decrement
  i=N;
  for(auto ibin = fend(myhist); ibin!=fbegin(myhist); )
    EXPECT_EQ((--ibin)->weight(), fival(--i,N));
  
  // Test post decrement
  i=N;
  for(auto ibin = fend(myhist); ibin!=fbegin(myhist); )
  {
    ibin--;
    i--;
    EXPECT_EQ(ibin->weight(), fival(i,N));
  }

  // Test equality and assignment operators
  auto ibin1 = fbegin(myhist);
  auto ibin2 = ibin1;
  ASSERT_EQ(ibin1,ibin2);
  ASSERT_TRUE(ibin1 == ibin2);
  ASSERT_FALSE(ibin1 != ibin2);
  ASSERT_TRUE(ibin1 <= ibin2);
  ASSERT_TRUE(ibin2 <= ibin1);
  ASSERT_TRUE(ibin1 >= ibin2);
  ASSERT_TRUE(ibin2 >= ibin1);
  ASSERT_FALSE(ibin1 < ibin2);
  ASSERT_FALSE(ibin1 > ibin2);
  ASSERT_FALSE(ibin2 < ibin1);
  ASSERT_FALSE(ibin2 > ibin1);

  for(i=0; i<10; ++i)++ibin1;
  ASSERT_NE(ibin1,ibin2);
  ASSERT_FALSE(ibin1 == ibin2);
  ASSERT_TRUE(ibin1 != ibin2);
  ASSERT_FALSE(ibin1 <= ibin2);
  ASSERT_TRUE(ibin2 <= ibin1);
  ASSERT_TRUE(ibin1 >= ibin2);
  ASSERT_FALSE(ibin2 >= ibin1);
  ASSERT_FALSE(ibin1 < ibin2);
  ASSERT_TRUE(ibin1 > ibin2);
  ASSERT_TRUE(ibin2 < ibin1);
  ASSERT_FALSE(ibin2 > ibin1);
  ASSERT_EQ(ibin1-ibin2, 10);
  ASSERT_EQ(ibin2-ibin1, -10);

  ibin2 += 10;
  ASSERT_EQ(ibin1->ibin(),ibin2->ibin());
  ibin2 -= 10;
  ASSERT_EQ(fbegin(myhist)->ibin(),ibin2->ibin());
  ibin2 = ibin2 + 10;
  ASSERT_EQ(ibin1->ibin(),ibin2->ibin());
  ibin2 = ibin2 - 10;
  ASSERT_EQ(fbegin(myhist)->ibin(),ibin2->ibin());
}

TEST(TestSimpleHist, Iterator) {
  iterator_test([](SimpleHist& h){return h.begin();},
                [](SimpleHist& h){return h.end();},
                [](double i, double N) { return i; });
}

TEST(TestSimpleHist, ConstIterator) {
  iterator_test([](SimpleHist& h){return h.cbegin();},
                [](SimpleHist& h){return h.cend();},
                [](double i, double N) { return i; });
}

TEST(TestSimpleHist, ReverseIterator) {
  iterator_test([](SimpleHist& h){return h.rbegin();},
                [](SimpleHist& h){return h.rend();},
                [](double i, double N) { return N-i-1; });
}

TEST(TestSimpleHist, ConstReverseIterator) {
  iterator_test([](SimpleHist& h){return h.crbegin();},
                [](SimpleHist& h){return h.crend();},
                [](double i, double N) { return N-i-1; });
}

TEST(TestSimpleHist, Moments) {
  unsigned N { 10000 };
  SimpleHist myhist {1.0};
  for(unsigned i=0;i<N;++i)myhist.insert(i,i);
  double dN { static_cast<double>(N) };
  ASSERT_EQ(myhist.sum_w(), (dN-1)*dN/2);
  ASSERT_EQ(myhist.sum_wx(), (dN-1)*dN*(2*dN-1)/6);
  ASSERT_EQ(myhist.sum_wxx(), (dN-1)*(dN-1)*dN*dN/4);
}

TEST(TestSimpleHist, IterateOverBins) {
  unsigned N { 10000 };
  SimpleHist myhist {1.0};
  for(unsigned i=0;i<N;++i)myhist.insert(i,i);
  double dN { static_cast<double>(N) };
  double sw {};
  double swx {};
  double swxx {};
  myhist.visit_bins([&sw,&swx,&swxx](double x,double w){ sw+=w; swx+=w*x; swxx+=w*x*x; });
  ASSERT_EQ(sw, (dN-1)*dN/2);
  ASSERT_EQ(swx, (dN-1)*dN*(2*dN-1)/6);
  ASSERT_EQ(swxx, (dN-1)*(dN-1)*dN*dN/4);
}

TEST(TestSimpleHist, Integrate) {
  unsigned N { 10000 };
  SimpleHist myhist {1.0};
  for(unsigned i=0;i<N;++i)myhist.insert(i,i);
  double dN { static_cast<double>(N) };
  double sw { myhist.integrate([](double x,double w) { return w; }) };
  double swx { myhist.integrate([](double x,double w) { return w*x; }) };
  double swxx { myhist.integrate([](double x,double w) { return w*x*x; }) };
  ASSERT_EQ(sw, (dN-1)*dN/2);
  ASSERT_EQ(swx, (dN-1)*dN*(2*dN-1)/6);
  ASSERT_EQ(swxx, (dN-1)*(dN-1)*dN*dN/4);
}

TEST(TestSimpleHist, KahanInsertIntegersFrom0ToN) {
  unsigned N { 1000 };
  BasicHistogram1D<KahanAccumulator> myhist {1.0};
  for(unsigned i=0;i<N;++i)
    for(unsigned j=i;j<N;++j)
      myhist.insert(j);
  for(unsigned i=0;i<N;++i)
    EXPECT_EQ(myhist.weight(i), i+1);
}

TEST(TestSimpleHist, CopyAndEquality) {
  unsigned N { 10000 };
  SimpleHist myhist {1.0};
  for(unsigned i=0;i<N;++i)myhist.insert(i,i);
  SimpleHist myhist2 { myhist };
  ASSERT_EQ(myhist,myhist2);
  myhist2.insert(N,N);
  ASSERT_FALSE(myhist == myhist2);
  myhist2 = myhist;
  ASSERT_EQ(myhist,myhist2);
  myhist2.set_name("Hello there");
  ASSERT_FALSE(myhist == myhist2);
}
  
TEST(TestSimpleHist, ToAndFromProtobuf) {
  unsigned N { 300 };
  SimpleHist myhist {1.0,99.5,199.5};
  for(unsigned i=0;i<N;++i)myhist.insert(i,i);
  myhist.set_name("Test histogram");
  myhist.set_xval_units("xval units");
  myhist.set_weight_units("weight units");
  Histogram1DData* hist_data = myhist.getData();
  SimpleHist myhist2(*hist_data);
  ASSERT_EQ(myhist,myhist2);
  //std::cout << hist_data->DebugString();
  delete hist_data;
}

TEST(TestSimpleHist, KahanToAndFromProtobuf) {
  unsigned N { 10000 };
  BasicHistogram1D<KahanAccumulator> myhist {1.0};
  for(unsigned i=0;i<N;++i)myhist.insert(i,i);
  Histogram1DData* hist_data = myhist.getData();
  BasicHistogram1D<KahanAccumulator> myhist2(*hist_data);
  // THIS COULD FAIL IN GENERAL AS CORRECTIONS AREN'T STORED IN THE
  // PROTOBUF DATA - IT SHOULD BE OK WITH THE VALUES USED HERE
  ASSERT_EQ(myhist,myhist2);
  delete hist_data;
}

TEST(TestBinnedCDF, CreateFromHistogram) {
  unsigned N { 1000 };
  SimpleHist myhist {1.0};
  for(unsigned i=0;i<N;++i)
    myhist.insert(i,i);
  BinnedCDF mycdf(myhist);
  ASSERT_EQ(mycdf.size(), N);
  for(unsigned i=0;i<N;++i)
    EXPECT_EQ(mycdf.xval_center(i), i);  
  for(unsigned i=0;i<N;++i)
    EXPECT_EQ(mycdf.cumulative_right(i), double(i*(i+1)/2)/double(N*(N-1)/2));
}

double median_zero_to_Nminus1(unsigned N) {
  return (N%2==0)?0.5*double((N-1)/2+(N/2)):double(N/2);
}

TEST(TestBinnedCDF, Stats) {
  unsigned N { 1000 }; // Should be even
  //BasicHistogram1D<KahanAccumulator> myhist {1.0};
  SimpleHist myhist {1.0};
  for(unsigned i=0;i<N;++i)
    myhist.insert(i,1.0);
  BinnedCDF mycdf(myhist);
  EXPECT_EQ(mycdf.median(), median_zero_to_Nminus1(N));
  myhist.insert(N,1.0);
  BinnedCDF mycdf2(myhist);
  EXPECT_EQ(mycdf2.median(), median_zero_to_Nminus1(N+1));
  double dN { static_cast<double>(N) };
  double swx;
  double swxx;
  //mycdf.moments2<KahanAccumulator>(swx,swxx);
  mycdf.moments2(swx,swxx);
  ASSERT_LE(std::abs(swx - (dN-1.0)/2.0),1e-6);
  ASSERT_LE(std::abs(swxx - (dN-1)*(2.0*dN-1.0)/6.0),1e-6);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
