#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>

#include "io/sql_transceiver.hpp"
#include "io/sqlite3_transceiver.hpp"
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

TEST(TestSQLTransceiver, PruneEmptyTables_NoEmptyTables) {
  SQLTransceiver xvr;
  SQLTransceiver::SQLTable* t =
      xvr.make_sqltable_tree("mytable", UnitTestMessage::descriptor());
  xvr.prune_empty_tables(t);
  xvr.iterate_over_tables(t, [](const SQLTransceiver::SQLTable* t) {
      EXPECT_GT(t->fields.size(), 0); } );
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
  ASSERT_LT(n_pruned, n_unpruned);
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

TEST(TestSQLTransceiver, CreateTableWithKey) {
  SQLTransceiver xvr;
  xvr.create_tables("mytable", UnitTestMessage::descriptor(),
                    UnitTestKey::descriptor(), "", true);
}

TEST(TestSQLTransceiver, Insert) {
  SQLTransceiver xvr;
  UnitTestMessage m_data;
  UnitTestKey m_key;
  m_key.set_user_key_i32(1234);
  m_key.set_user_key_string("key string value");
  m_key.mutable_user_key_ssm()->set_ssm_i32(2345);

  m_data.set_i32(1);
  m_data.set_i64(2);
  m_data.set_f(3.14);
  m_data.set_s("string4");
  m_data.mutable_ssm()->set_ssm_i32(5111);
  m_data.mutable_ssm_inline()->set_ssm_i32(6111);

  m_data.mutable_csm()->set_csm_i32(7111);
  m_data.mutable_csm()->mutable_csm_ssm()->set_ssm_i32(712111);
  m_data.mutable_csm_inline()->set_csm_i32(8111);
  m_data.mutable_csm_inline()->mutable_csm_ssm()->set_ssm_i32(812111);
  m_data.mutable_ism()->set_ism_i32(9111);
  m_data.mutable_ism()->mutable_ism_ssm()->set_ssm_i32(912111);
  m_data.mutable_ism_inline()->set_ism_i32(10111);
  m_data.mutable_ism_inline()->mutable_ism_ssm()->set_ssm_i32(1012111);

#if 0
  UnitTestInlinedSubMessage ism = 9;
  UnitTestInlinedSubMessage ism_inline = 10 [(CFO).sql.inline_message = true];
  UnitTestVectorSubMessage vsm = 11;
  UnitTestVectorSubMessage vsm_inline = 12 [(CFO).sql.inline_message = true];
#endif

  
  xvr.insert("mytable", &m_data, &m_key, true);
}

TEST(TestSQLite3Transceiver, CreateTableWithKey) {
  SQLite3Transceiver xvr("test_db.sqlite", SQLite3Transceiver::TRUNCATE_RW);
  //SQLite3Transceiver xvr(":memory:", SQLite3Transceiver::TRUNCATE_RW);
  xvr.create_tables("mytable", UnitTestMessage::descriptor(),
                    UnitTestKey::descriptor(), "", true);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
