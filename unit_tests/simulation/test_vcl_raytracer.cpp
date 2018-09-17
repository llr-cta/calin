/*

   calin/unit_tests/io/test_vcl_raytracer.cpp -- Stephen Fegan -- 2018-09-17

   Unit tests for VCL raytracer

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
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
#include <iomanip>
#include <fstream>
#include <gtest/gtest.h>

#include <util/vcl.hpp>
#include <simulation/vcl_raytracer.hpp>

using namespace calin::util::vcl;
using namespace calin::simulation::vcl_raytracer;

template<typename VCLRealType> class TestVCLRaytracer :
  public VCLRealType, public testing::Test
{
public:
};

using RealTypes = ::testing::Types<VCL128FloatReal, VCL256FloatReal, VCL512FloatReal, VCL128DoubleReal, VCL256DoubleReal, VCL512DoubleReal>;

TYPED_TEST_CASE(TestVCLRaytracer, RealTypes);

TYPED_TEST(TestVCLRaytracer, Test) {
  ScopeRayTracer<TypeParam> raytracer;

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
