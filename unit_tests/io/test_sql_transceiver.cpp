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
  m_data.mutable_ssm()->set_ssm_i32(51);
  m_data.mutable_ssm_inline()->set_ssm_i32(611);

  m_data.mutable_csm()->set_csm_i32(71);
  m_data.mutable_csm()->mutable_csm_ssm()->set_ssm_i32(721);
  m_data.mutable_csm_inline()->set_csm_i32(81);
  m_data.mutable_csm_inline()->mutable_csm_ssm()->set_ssm_i32(821);
  m_data.mutable_ism()->set_ism_i32(91);
  m_data.mutable_ism()->mutable_ism_ssm()->set_ssm_i32(921);
  m_data.mutable_ism_inline()->set_ism_i32(101);
  m_data.mutable_ism_inline()->mutable_ism_ssm()->set_ssm_i32(1021);
#if 0
  for(unsigned i=0;i<10;i++)
  {
    m_data.add_vsm()->set_ssm_i32(11001+i*10);
    m_data.add_vsm_inline()->set_ssm_i32(12001+i*10);
  }
#endif
  for(unsigned i=0;i<10;i++)
  {
    m_data.add_vec_i32(i*1000000 + 101);
  }

#if 0
    repeated int64   vec_i64  = 102;
  repeated float   vec_f    = 103;
  repeated string  vec_s    = 104;
  repeated UnitTestSimpleSubMessage vec_ssm = 105;
  repeated UnitTestSimpleSubMessage vec_ssm_inline = 106 [(CFO).sql.inline_message = true];
#endif

  m_data.set_oo_s("OO string test");


  m_data.mutable_oosm_inline()->set_oosm_i32(1234567);
  m_data.mutable_oosm_inline()->set_oosm_s("OO test string 2");
  m_data.mutable_oosm_inline()->mutable_oosm_ssm_inline()->set_ssm_i32(987654);
      
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
