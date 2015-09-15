#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>

#include "io/sql_transceiver.hpp"
#include "proto/unittest.pb.h"

using namespace calin::ix::unittest;
using namespace calin::io::sql_transceiver;

TEST(TestSQLTransceiver, ParentChildPointers) {
  SQLTransceiver xvr;
  SQLTransceiver::SQLTable* t =
      xvr.make_sqltable_tree("mytable", UnitTestMessage::descriptor());
  xvr.iterate_over_tables(t, [](const SQLTransceiver::SQLTable* t) {
      for(auto st : t->sub_tables) { EXPECT_EQ(st->parent_table, t); } });
}

TEST(TestSQLTransceiver, FieldTablePointers) {
  SQLTransceiver xvr;
  SQLTransceiver::SQLTable* t =
      xvr.make_sqltable_tree("mytable", UnitTestMessage::descriptor());
  xvr.iterate_over_tables(t, [](const SQLTransceiver::SQLTable* t) {
      for(auto f : t->fields) { EXPECT_EQ(f->table, t); } });
}  

TEST(TestSQLTransceiver, PruneEmptyTables_ParentChildPointers) {
  SQLTransceiver xvr;
  SQLTransceiver::SQLTable* t =
      xvr.make_sqltable_tree("mytable", UnitTestMessage::descriptor());
  xvr.prune_empty_tables(t);
  xvr.iterate_over_tables(t, [](const SQLTransceiver::SQLTable* t) {
      for(auto st : t->sub_tables) { EXPECT_EQ(st->parent_table, t); } });
}

TEST(TestSQLTransceiver, PruneEmptyTables_FieldTablePointers) {
  SQLTransceiver xvr;
  SQLTransceiver::SQLTable* t =
      xvr.make_sqltable_tree("mytable", UnitTestMessage::descriptor());
  xvr.prune_empty_tables(t);
  xvr.iterate_over_tables(t, [](const SQLTransceiver::SQLTable* t) {
      for(auto f : t->fields) { EXPECT_EQ(f->table, t); } });
}  

TEST(TestSQLTransceiver, PruneEmptyTables_NoFieldLeftBehind) {
  SQLTransceiver xvr;
  SQLTransceiver::SQLTable* t =
      xvr.make_sqltable_tree("mytable", UnitTestMessage::descriptor());
  auto tf_unpruned = xvr.list_all_table_columns(t);
  xvr.prune_empty_tables(t);
  auto tf_pruned = xvr.list_all_table_columns(t);
  ASSERT_EQ(tf_unpruned, tf_pruned);
}

TEST(TestSQLTransceiver, PruneEmptyTables_SomeTablesDeleted) {
  SQLTransceiver xvr;
  SQLTransceiver::SQLTable* t =
      xvr.make_sqltable_tree("mytable", UnitTestMessage::descriptor());
  unsigned n_unpruned = 0;
  xvr.iterate_over_tables(t, [&n_unpruned](const SQLTransceiver::SQLTable* t) {
      n_unpruned++; });
  xvr.prune_empty_tables(t);
  unsigned n_pruned = 0;
  xvr.iterate_over_tables(t, [&n_pruned](const SQLTransceiver::SQLTable* t) {
      n_pruned++; });
  ASSERT_GT(n_unpruned, n_pruned);
}

TEST(TestSQLTransceiver, MakeColumns) {
  SQLTransceiver xvr;
  SQLTransceiver::SQLTable* t =
      xvr.make_sqltable_tree("mytable", UnitTestMessage::descriptor());
  xvr.propagate_keys(t);
  auto tf = xvr.list_all_table_columns(t);
  for(auto itf : tf)
    std::cout << itf.first << ' ' << itf.second << '\n';
}

TEST(TestSQLTransceiver, CreateTable) {
  SQLTransceiver xvr;
  xvr.create_tables("mytable", UnitTestMessage::descriptor());
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
