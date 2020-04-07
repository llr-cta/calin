/*

   calin/util/log.hpp -- Stephen Fegan -- 2015-05-12

   Class providing logging capabilities

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

#include <Python.h>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include <time.h>
#include <memory>

#include <util/log.hpp>

using namespace calin::util::log;

Logger::~Logger()
{
  //nothing to see here
}

namespace {

const char* ansi_color_string(Level level)
{
  switch(level)
  {
    case Level::FATAL:    return "\x1b[37;41;97;101;1m";
    case Level::ERROR:    return "\x1b[30;45;105;1m";
    case Level::WARNING:  return "\x1b[30;43;103;1m";
      //case Level::SUCCESS:  return "\x1b[37;42;97;1m";
    case Level::SUCCESS:  return "\x1b[37;42;97;38;5;15;1m";
    case Level::FAILURE:  return "\x1b[31;91;47;1m";
      //case Level::FAILURE:  return "\x1b[30;41;101;1m";
    case Level::INFO:     return "\x1b[37;46;97;1m";
    case Level::VERBOSE:  return "";
    case Level::DISCARD:  return "";
    default: assert(0); return "";
  }
}

const char* level_string(Level level)
{
  switch(level)
  {
    case Level::FATAL:    return "FATAL";
    case Level::ERROR:    return "ERROR";
    case Level::WARNING:  return "WARNING";
    case Level::FAILURE:  return "FAILURE";
    case Level::SUCCESS:  return "SUCCESS";
    case Level::INFO:     return "INFO";
    case Level::VERBOSE:  return "";
    case Level::DISCARD:  return "";
    default: assert(0); return "";
  }
}

const char* ansi_reset_string()
{
  return "\x1b[0m";
}

template <typename Writer> void
write_message_lines(Writer &&writer,
                    const char* timestamp_string,
                    const char* this_level_string,
                    const char* apply_color_string,
                    const char* reset_color_string,
                    const std::string& message)
{
  size_t spos = 0;
  while(spos < message.size())
  {
    size_t epos = message.find_first_of('\n',spos);
    unsigned n = (epos==std::string::npos ? message.size() : epos) - spos;

    if(timestamp_string and *timestamp_string) {
      writer(timestamp_string, std::strlen(timestamp_string));
      writer(" ",1);
    }
    if(apply_color_string and *apply_color_string and
       this_level_string and *this_level_string) {
      writer(apply_color_string, std::strlen(apply_color_string));
      writer(" ",1);
      writer(this_level_string, std::strlen(this_level_string));
      writer(" ",1);
      writer(reset_color_string, std::strlen(reset_color_string));
      writer(" ",1);
    }
    else if(this_level_string and *this_level_string) {
      writer("[",1);
      writer(this_level_string, std::strlen(this_level_string));
      writer("] ",2);
    }

    writer(message.c_str() + spos, n);
    writer("\n",1);

    spos += n+1;
  }
}

std::unique_ptr<ProtobufLogger> default_protobuf_logger_ { new ProtobufLogger };

} // anonymous namespace

ProtobufLogger* calin::util::log::default_protobuf_logger()
{
  return default_protobuf_logger_.get();
}

void calin::util::log::prune_default_protobuf_log()
{
  default_protobuf_logger()->prune_log();
}

MultiLogger::~MultiLogger()
{
  //nothing to see here
}

void MultiLogger::clear_all_loggers_and_streams()
{
  sub_loggers_.clear();
  sub_streams_.clear();
}

void MultiLogger::log_message(Level level, const std::string& message,
  calin::util::timestamp::Timestamp timestamp)
{
  if(message.empty() || level==DISCARD)return;

  std::lock_guard<decltype(lockable_)> lockable_guard(lockable_);

  if(sub_loggers_.empty() and sub_streams_.empty())
  {
    nolock_add_logger(default_protobuf_logger(), false);
    if(PythonLogger::is_python_initialised())
      nolock_add_logger(new PythonLogger, true);
    else
      nolock_add_stream(&std::cout, false, true, true);
  }

  for(auto& sl : sub_loggers_)
  {
    try
    {
      sl.logger_->log_message(level, message, timestamp);
    }
    catch(...)
    {
      // for better or worse - ignore all errors
    }
  }

  std::string timestamp_string = timestamp.as_string();
  const char* this_level_string = level_string(level);
  const char* apply_color_string = ansi_color_string(level);
  const char* reset_color_string = ansi_reset_string();

  for(auto& ss : sub_streams_)
  {
    try
    {
      std::ostream* stream { ss.stream_ };
      write_message_lines([stream](const char* c, unsigned n) {
          stream->write(c,n); },
        ss.apply_timestamp_?timestamp_string.c_str():nullptr, this_level_string,
        ss.use_colors_?apply_color_string:nullptr, reset_color_string,
        message);
    }
    catch(...)
    {
      // for better or worse - ignore all errors
    }
  }
}

void MultiLogger::add_logger(Logger* logger, bool adopt_logger)
{
  std::lock_guard<decltype(lockable_)> lockable_guard(lockable_);
  nolock_add_logger(logger, adopt_logger);
}

void MultiLogger::add_stream(std::ostream* stream, bool adopt_stream,
                             bool apply_timestamp, bool use_colors)
{
  std::lock_guard<decltype(lockable_)> lockable_guard(lockable_);
  nolock_add_stream(stream, adopt_stream, apply_timestamp, use_colors);
}

void MultiLogger::nolock_add_logger(Logger* logger, bool adopt_logger)
{
  sub_loggers_.emplace_back(logger, adopt_logger);
}

void MultiLogger::nolock_add_stream(std::ostream* stream, bool adopt_stream,
                                      bool apply_timestamp, bool use_colors)
{
  sub_streams_.emplace_back(stream, adopt_stream,
    apply_timestamp, use_colors);
}

void MultiLogger::add_cout(bool apply_timestamp, bool use_colors)
{
  add_stream(&std::cout, false, apply_timestamp, use_colors);
}

void MultiLogger::add_cerr(bool apply_timestamp, bool use_colors)
{
  add_stream(&std::cerr, false, apply_timestamp, use_colors);
}

void MultiLogger::add_clog(bool apply_timestamp, bool use_colors)
{
  add_stream(&std::clog, false, apply_timestamp, use_colors);
}

bool MultiLogger::add_file(const std::string filename,
                           bool apply_timestamp, bool use_colors)
{
  std::ofstream* stream = new std::ofstream(filename.c_str());
  if(*stream)
  {
    add_stream(stream, true, apply_timestamp, use_colors);
    return true;
  }
  else
  {
    delete stream;
    return false;
  }
}

ProtobufLogger::~ProtobufLogger()
{
  // nothing to see here
}

void ProtobufLogger::
log_message(Level level, const std::string& message,
  calin::util::timestamp::Timestamp timestamp)
{
  auto* m = log_.add_message();
  switch(level)
  {
  case FATAL:
    m->set_level(calin::ix::util::log::LogMessage::FATAL); break;
  case ERROR:
    m->set_level(calin::ix::util::log::LogMessage::ERROR); break;
  case WARNING:
    m->set_level(calin::ix::util::log::LogMessage::WARNING); break;
  case INFO:
    m->set_level(calin::ix::util::log::LogMessage::INFO); break;
  case SUCCESS:
    m->set_level(calin::ix::util::log::LogMessage::SUCCESS); break;
  case FAILURE:
    m->set_level(calin::ix::util::log::LogMessage::FAILURE); break;
  case VERBOSE:
    m->set_level(calin::ix::util::log::LogMessage::VERBOSE); break;
  case DISCARD:
    m->set_level(calin::ix::util::log::LogMessage::DISCARD); break;
  }
  m->set_message(message);
  timestamp.as_proto(m->mutable_timestamp());
}

PythonLogger::~PythonLogger()
{
  // nothing to see here
}

void PythonLogger::
log_message(Level level, const std::string& message,
  calin::util::timestamp::Timestamp timestamp)
{
  if(message.empty() || level==DISCARD)return;

  const char* this_level_string = level_string(level);
  const char* apply_color_string = ansi_color_string(level);
  const char* reset_color_string = ansi_reset_string();

#if 1
  /* Would seem like we should get the GIL before calling the write but
     this seems to block if we're calling from a different (non Python)
     thread. Best to just leave it for now */
  PyGILState_STATE gstate;
  gstate = PyGILState_Ensure();
#endif

  if(use_stderr_)
    write_message_lines([](const char* c, unsigned n) { while(n) {
          unsigned nn = std::min(n,1000U); PySys_WriteStderr("%.*s", nn, c);
          c += nn; n -= nn; } },
      nullptr, this_level_string, apply_color_string, reset_color_string,
      message);
  else
    write_message_lines([](const char* c, unsigned n) { while(n) {
          unsigned nn = std::min(n,1000U); PySys_WriteStdout("%.*s", nn, c);
          c += nn; n -= nn; } },
      nullptr, this_level_string, apply_color_string, reset_color_string,
      message);

#if 1
  /* Release the thread. No Python API allowed beyond this point. */
  PyGILState_Release(gstate);
#endif
}

bool PythonLogger::is_python_initialised()
{
  return Py_IsInitialized() != 0;
}
