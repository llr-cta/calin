#include "io/log.hpp"
#include "io/sqlite3_transceiver.hpp"

using namespace calin::io::log;
using namespace calin::io::sql_transceiver;

SQLite3Transceiver::SQLite3Transceiver(const std::string& filename):
    SQLTransceiver(), db_(), inherit_db_()
{
  sqlite3_open_v2(filename.c_str(), &db_,
                   SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE,
                   NULL);
}

SQLite3Transceiver::~SQLite3Transceiver()
{
  if(inherit_db_)sqlite3_close(db_);
}

bool SQLite3Transceiver::
create_tables(const std::string& table_name,
              const google::protobuf::Descriptor* d_data,
              const google::protobuf::Descriptor* d_key,
              const std::string& instance_desc,
              bool write_sql_to_log)
{
  SQLTable* t = make_keyed_sqltable_tree(table_name, d_data, d_key, false);
  iterate_over_tables(t, [this,write_sql_to_log](const SQLTable* it) {
      sqlite3_stmt* stmt { nullptr };
      int retcode { 0 };
      retcode = sqlite3_prepare_v2(db_, sql_create_table(it).c_str(), -1,
                                   &stmt, nullptr);
      if(retcode != SQLITE_OK)
      {
        LOG(ERROR) << "SQLITE error: " << sqlite3_errmsg(db_);
	    goto finalize;
      }
      retcode = sqlite3_step(stmt);
      if(retcode == SQLITE_ERROR)
        LOG(ERROR) << "SQLITE error: " << sqlite3_errmsg(db_) << '\n';
      else if(retcode == SQLITE_MISUSE)
        LOG(ERROR) << "SQLITE misuse" << '\n';
   finalize:
      sqlite3_finalize(stmt);      
      if(write_sql_to_log)
        LOG(VERBOSE) << sql_create_table(it) << '\n';
    });
  delete t;
  return true;
}
