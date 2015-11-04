/* 

   calin/unit_tests/io/test_sql_transceiver.cpp -- Stephen Fegan -- 2015-09-13

   Unit tests for sql_transceiver and sqlite3_transceiver classes

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
#include <fstream>
#include <gtest/gtest.h>
#include <google/protobuf/text_format.h>

#include "io/sql_transceiver.hpp"
#include "io/sqlite3_transceiver.hpp"
#include "proto/unittest.pb.h"
#include <google/protobuf/text_format.h>

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
      EXPECT_GT(t->fields.size(), 0U); } );
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

  xvr.iterate_over_tables(t, [](const SQLTransceiver::SQLTable* t) {
      for(auto f : t->fields) { 
        std::cout << t->table_name << ' ' << f->field_name;
        if(f->field_d)std::cout << ' ' << f->field_d->number();
        std::cout << '\n';
      } });
}

TEST(TestSQLTransceiver, CreateTable) {
  SQLTransceiver xvr;
  xvr.create_tables("mytable", UnitTestMessage::descriptor());
}

TEST(TestSQLTransceiver, CreateTableWithKey) {
  SQLTransceiver xvr;
  xvr.create_tables("mytable", UnitTestMessage::descriptor(),
                    UnitTestKey::descriptor(), "Instance of mytable");
}

TEST(TestSQLite3Transceiver, CreateTableWithKey) {
  SQLite3Transceiver xvr("test_db.sqlite", SQLite3Transceiver::TRUNCATE_RW);
  //SQLite3Transceiver xvr(":memory:", SQLite3Transceiver::TRUNCATE_RW);
  xvr.create_tables("mytable", UnitTestMessage::descriptor(),
                    UnitTestKey::descriptor(), "Instance of mytable");
}

TEST(TestSQLite3Transceiver, Insert) {
  SQLite3Transceiver xvr("test_db.sqlite", SQLite3Transceiver::TRUNCATE_RW);
  //SQLite3Transceiver xvr(":memory:", SQLite3Transceiver::TRUNCATE_RW);
  xvr.create_tables("mytable", UnitTestMessage::descriptor(),
                    UnitTestKey::descriptor(), "Instance of mytable");
  
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
    m_data.add_vec_i64(i*1000000 + 102);
    m_data.add_vec_f((i*1000000 + 103)+0.1);
    m_data.add_vec_s(std::string("string ")+std::to_string(i*1000000 + 104));
    m_data.add_vec_ssm()->set_ssm_i32(i*1000000 + 105);
    m_data.add_vec_ssm_inline()->set_ssm_i32(i*1000000 + 106);
    
    m_data.add_vec_d(i*1000000 + 113.5);
    m_data.add_vec_ui32(i*1000000 + 114);
    m_data.add_vec_ui64(i*1000000 + 115);
    m_data.add_vec_si32(i*1000000 + 116);
    m_data.add_vec_si64(i*1000000 + 117);
    m_data.add_vec_fi32(i*1000000 + 118);
    m_data.add_vec_fi64(i*1000000 + 119);
    m_data.add_vec_sfi32(i*1000000 + 120);
    m_data.add_vec_sfi64(i*1000000 + 121);
    m_data.add_vec_b(i%2 == 0);
    m_data.add_vec_bb(std::string("string ")+std::to_string(i*1000000 + 123));

    if(i%2==0)m_data.add_vec_e(UnitTestMessage::RUNNING);
    else m_data.add_vec_e(UnitTestMessage::STARTED);
        
    auto ivec_csm = m_data.add_vec_csm();
    ivec_csm->set_csm_i32(i*1000000 + 107);
    ivec_csm->mutable_csm_ssm()->set_ssm_i32(i*1000000 + 1107);

    auto ivec_csm_inline = m_data.add_vec_csm_inline();
    ivec_csm_inline->set_csm_i32(i*1000000 + 108);
    ivec_csm_inline->mutable_csm_ssm()->set_ssm_i32(i*1000000 + 1108);

    auto ivec_ism = m_data.add_vec_ism();
    ivec_ism->set_ism_i32(i*1000000 + 109);
    ivec_ism->mutable_ism_ssm()->set_ssm_i32(i*1000000 + 1109);

    auto ivec_ism_inline = m_data.add_vec_ism_inline();
    ivec_ism_inline->set_ism_i32(i*1000000 + 110);
    ivec_ism_inline->mutable_ism_ssm()->set_ssm_i32(i*1000000 + 1110);

    auto ivec_vsm = m_data.add_vec_vsm();
    for(unsigned j=0;j<i;j++)
    {
      ivec_vsm->add_vsm_vec_i32(i*1000000 + j*10000 + 1111);
      ivec_vsm->add_vsm_vec_ssm()->set_ssm_i32(i*1000000 + j*10000 + 2111);
    }

    m_data.add_vec_oosm()->set_oosm_i32(i*1000000 + 4031);
    m_data.add_vec_oosm()->set_oosm_s(std::string("oo string ") +
                                      std::to_string(i*1000000 + 4032));
    m_data.add_vec_oosm()->mutable_oosm_ssm()->set_ssm_i32(i*1000000 + 40331);
    m_data.add_vec_oosm()->mutable_oosm_ssm_inline()->
        set_ssm_i32(i*1000000 + 40341);

    std::string map_key = std::string("map key ")+std::to_string(i);
    (*m_data.mutable_map_i32())[std::string("i32 ")+map_key] = i*1000000 + 201;
    (*m_data.mutable_map_ssm())[std::string("ssm ")+map_key].
        set_ssm_i32(i*1000000 + 202); 
    auto &map_csm = (*m_data.mutable_map_csm())[std::string("csm ")+map_key];
    map_csm.set_csm_i32(i*1000000 + 203); 
    map_csm.mutable_csm_ssm()->set_ssm_i32(i*100000 + 1203);

    auto &map_ism = (*m_data.mutable_map_ism())[std::string("ism ")+map_key];
    map_ism.set_ism_i32(i*1000000 + 204); 
    map_ism.mutable_ism_ssm()->set_ssm_i32(i*100000 + 1204); 

    auto &map_vsm = (*m_data.mutable_map_vsm())[std::string("vsm ")+map_key];
    for(unsigned j=0;j<=i;j++)
    {
      map_vsm.add_vsm_vec_i32(i*1000000 + j*10000 + 211);
      map_vsm.add_vsm_vec_ssm()->set_ssm_i32(i*1000000 + j*10000 + 2211);
    }
  }

  m_data.set_d(1300000.1);
  m_data.set_ui32(1400000);
  m_data.set_ui64(1500000);
  m_data.set_si32(1600000);
  m_data.set_si64(1700000);
  m_data.set_fi32(1800000);
  m_data.set_fi64(1900000);
  m_data.set_sfi32(200000);
  m_data.set_sfi64(210000);
  m_data.set_b(true);
  m_data.set_bb("Bytes bytes bytes");
  m_data.set_e(UnitTestMessage::RUNNING);

  
  m_data.set_oo_s("OO string test");

  m_data.mutable_oosm_inline()->set_oosm_i32(1234567);
  m_data.mutable_oosm_inline()->set_oosm_s("OO test string 2");
  m_data.mutable_oosm_inline()->mutable_oosm_ssm_inline()->set_ssm_i32(987654);

  uint64_t oid;

  std::ofstream mystream("test_db_insert.txt");
  std::string str;
  google::protobuf::TextFormat::PrintToString(m_data, &str);
  mystream << str;

  xvr.insert("mytable", oid, &m_data, &m_key);
  //std::cerr << oid << '\n';
}

TEST(TestSQLite3Transceiver, RetreiveByOID) {
  SQLite3Transceiver xvr("test_db.sqlite", SQLite3Transceiver::READ_ONLY);
  //SQLite3Transceiver xvr(":memory:", SQLite3Transceiver::TRUNCATE_RW);

#if 0
  xvr.create_tables("mytable", UnitTestMessage::descriptor(),
                    UnitTestKey::descriptor(), "Instance of mytable");
#endif
  
  UnitTestMessage m_data;
  UnitTestKey m_key;

  xvr.retrieve_by_oid("mytable", 1, &m_data, &m_key);

  std::ofstream mystream("test_db_select.txt");
  std::string str;
  google::protobuf::TextFormat::PrintToString(m_data, &str);
  mystream << str;
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
