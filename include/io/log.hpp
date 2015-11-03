/* 

   calin/io/log.hpp -- Stephen Fegan -- 2015-05-12

   Class providing logging capabilities

   Copyright 2015, Stephen Fegan <sfegan@gmail.com>

   This file is part of "calin"
   
   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.
    
   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#pragma once

#include <string>
#include <streambuf>
#include <iostream>
#include <sstream>
#include <vector>

namespace calin { namespace io { namespace log {

enum Level { FATAL, ERROR, WARNING, INFO, SUCCESS, FAILURE, VERBOSE, DISCARD };

struct TimeStamp
{
  std::string string() const;
  static TimeStamp now();
  double seconds_since(const TimeStamp& then);
  uint64_t sec;
  uint32_t usec;
};

class Logger
{
 public:
  virtual ~Logger();
  virtual void log_message(Level level, const std::string& message,
                           TimeStamp timestamp = TimeStamp::now()) = 0;
};

class MultiLogger: public Logger
{
 public:
  virtual ~MultiLogger();
  void log_message(Level level, const std::string& message,
                   TimeStamp timestamp = TimeStamp::now()) override;

  void add_logger(Logger* logger, bool adopt_logger = false);
  void add_stream(std::ostream* stream, bool adopt_stream = false,
                  bool apply_timestamp = false, bool use_colors = false);

  void add_cout(bool apply_timestamp = false, bool use_colors = false);
  void add_cerr(bool apply_timestamp = false, bool use_colors = false);
  void add_clog(bool apply_timestamp = false, bool use_colors = false);
  bool add_file(const std::string filename,
                bool apply_timestamp = false, bool use_colors = false);
  
 protected:
  void lock();
  void unlock();

  void nolock_add_logger(Logger* logger, bool adopt_logger);
  void nolock_add_stream(std::ostream* stream, bool adopt_stream,
                           bool apply_timestamp, bool use_colors);

  
  class sub_logger
  {
   public:
    sub_logger(Logger* logger_, bool adopt_logger_):
        logger(logger_), adopt_logger(adopt_logger_) { }
    Logger* logger = nullptr;
    bool adopt_logger = false;
  };
    
  struct sub_stream
  {
    sub_stream(std::ostream* stream_, bool adopt_stream_, bool apply_timestamp_,
               bool use_colors_):
        stream(stream_), adopt_stream(adopt_stream_),
        apply_timestamp(apply_timestamp_), use_colors(use_colors_) { }
    std::ostream* stream = nullptr;
    bool adopt_stream = false;
    bool apply_timestamp = false;
    bool use_colors = false;
  };

  std::vector<sub_logger> sub_loggers_;
  std::vector<sub_stream> sub_streams_;
};

inline MultiLogger* default_logger()
{
  static MultiLogger s_logger;
  return &s_logger;
}

class PythonLogger: public Logger
{
 public:
  PythonLogger(bool use_stderr = false): Logger(), use_stderr_(use_stderr) { }
  virtual ~PythonLogger();
  void log_message(Level level, const std::string& message,
                   TimeStamp timestamp = TimeStamp::now()) override;
  static bool is_python_initialised();
 protected:
  bool use_stderr_ { false };
};

enum class LoggerStreamEnabledState { LEAVE_AS_IS, ENABLE, DISABLE, TOGGLE };

class LoggerStreamManipulorSetEnabledState
{
 public:
  LoggerStreamManipulorSetEnabledState(LoggerStreamEnabledState state):
      state_(state) { /* nothing to see here */ }
  LoggerStreamEnabledState state() const { return state_; }
 private:
  LoggerStreamEnabledState state_;
};

inline LoggerStreamManipulorSetEnabledState enable_logging() {
  return LoggerStreamManipulorSetEnabledState(LoggerStreamEnabledState::ENABLE); }

inline LoggerStreamManipulorSetEnabledState disable_logging() {
  return LoggerStreamManipulorSetEnabledState(LoggerStreamEnabledState::DISABLE); }

inline LoggerStreamManipulorSetEnabledState toggle_logging() {
  return LoggerStreamManipulorSetEnabledState(LoggerStreamEnabledState::TOGGLE); }

inline LoggerStreamManipulorSetEnabledState enable_logging_if(bool clause) {
  return LoggerStreamManipulorSetEnabledState(clause?
      LoggerStreamEnabledState::ENABLE:
      LoggerStreamEnabledState::LEAVE_AS_IS); }

inline LoggerStreamManipulorSetEnabledState disable_logging_if(bool clause) {
  return LoggerStreamManipulorSetEnabledState(clause?
      LoggerStreamEnabledState::DISABLE:
      LoggerStreamEnabledState::LEAVE_AS_IS); }

class LoggerStream
{
 public:
  LoggerStream(Logger* logger, Level level):
      logger_(logger), level_(level), message_(new std::ostringstream) { }
  LoggerStream(LoggerStream&& o):
      logger_(o.logger_), level_(o.level_), message_(o.message_) {
    o.message_ = nullptr; }
  ~LoggerStream()
  {
    if(message_ != nullptr && !message_->str().empty())
      logger_->log_message(level_, message_->str());
    delete message_;
  }
  template<typename T> LoggerStream& operator<< (const T& t)
  {
    if(enabled_)*message_ << t;
    return *this;
  }
  std::ostream& stream() { return *message_; }
  void set_logging_enabled_state(LoggerStreamEnabledState state)
  {
    switch(state) {
      case LoggerStreamEnabledState::LEAVE_AS_IS: break;
      case LoggerStreamEnabledState::ENABLE:  enabled_ = true; break;
      case LoggerStreamEnabledState::DISABLE: enabled_ = false; break;
      case LoggerStreamEnabledState::TOGGLE:  enabled_ = !enabled_; break;
    }
  }
  LoggerStream& operator<< (const LoggerStreamManipulorSetEnabledState& m)
  {
    set_logging_enabled_state(m.state());
    return *this;
  }

  void enable_logging() { enabled_ = true; }
  void disable_logging() { enabled_ = false; }
  bool is_logging_enabled() const { return enabled_; }
 protected:
  Logger* logger_;
  Level level_;
  std::ostringstream* message_;
  bool enabled_ { true };
};

class VaporStream
{
 public:
  template<typename T> VaporStream& operator<< (const T&)
  {
    return *this;
  }
};

inline LoggerStream LOG(Level level, Logger* logger = default_logger())
{  
  return LoggerStream(logger, level);
}

#if 0
template<typename Streamer>
const char* printv(Streamer& streamer, const char* format)
{
  
}

template<typename Streamer, typename T>
const char* printv(Streamer& streamer, const char* format, const T& x)
{
  std::ostringstream str;
  
  while(*format)
  {
    if(*format == '%')
    {
      format++;
      if(*format == '-')
      {
        str << std::left;
        format++;
      }
    }
    else
    {
      str << *format;
      format++;
    }
}

template<typename Streamer, typename T, typename... Params>
const char* printv(Streamer& streamer, const char* format, const T& x,
                   const Params & ... params)
{
  format = printv(streamer, format, x);
  if(*format) format = printv(streamer, format, params...);
  return format;
}
#endif

inline unsigned num_digits(unsigned x)
{
  unsigned n = 1;
  while(x/=10)n++;
  return n;
}

} } } // namespace calin::io::log

