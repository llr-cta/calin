#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>

#include "io/sql_transceiver.hpp"
#include "proto/unittest.pb.h"

using namespace calin::ix::unittest;
using namespace calin::io::sql_transceiver;

TEST(TestSQLTransceiver, MakeColumns) {
  auto tf = SQLTransceiver::list_all_table_columns("mytable",
                                                    UnitTestMessage::descriptor());
  for(auto itf : tf)
    std::cout << itf.first << ' ' << itf.second << '\n';
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
