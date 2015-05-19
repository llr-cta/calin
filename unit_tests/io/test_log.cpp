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

