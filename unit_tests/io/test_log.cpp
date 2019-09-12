/*

   calin/unit_tests/io/test_log.cpp -- Stephen Fegan -- 2015-05-13

   Unit tests for logging class

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

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
#include <fstream>
#include <iomanip>
#include <thread>
#include <gtest/gtest.h>

#include "util/log.hpp"

using namespace calin::util::log;

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

void hammer_log(unsigned id)
{
  for(unsigned i=0;i<10000;i++)
    LOG(INFO) << "The quick red fox jumped over the lazy brown dog. " << id;
}

TEST(TestLog, MultiThreaded) {
  default_logger()->clear_all_loggers_and_streams();
  default_logger()->add_file("testlog.txt", false, false);
  std::vector<std::thread> threads;
  for(unsigned i=0;i<10;i++)threads.emplace_back(hammer_log, i);
  for(auto& i : threads)i.join();
  default_logger()->clear_all_loggers_and_streams();
  std::ifstream stream("testlog.txt");
  std::string line;
  while(std::getline(stream, line))
    ASSERT_EQ(line.size(), 58U) << "Unexpected line size: " << line;
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
