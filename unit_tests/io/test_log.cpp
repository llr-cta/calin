/* 

   calin/unit_tests/io/test_log.cpp -- Stephen Fegan -- 2015-05-13

   Unit tests for logg classe

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

#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>

#include "io/log.hpp"

using namespace calin::io::log;

TEST(TestLog, WriteStream) {
  default_logger()->add_cerr(true,true);
  LOG(FATAL) << "Fatal!\n";
  LOG(ERROR) << "Error!\n";
  LOG(WARNING) << "Warning!\n";
  LOG(INFO) << "Info!\n";
  LOG(SUCCESS) << "Success :-)\n";
  LOG(FAILURE) << "Failure :-(\n";
  LOG(VERBOSE) << "Multi-line\n   hello!\n";
  LOG(VERBOSE) << "Unfinished line";
  LOG(VERBOSE) << "Second unfinished line";
  LOG(VERBOSE) << "Unfinished multi-line\n   hello";
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

