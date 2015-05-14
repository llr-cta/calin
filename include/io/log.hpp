/* 

   calin/io/log.hpp -- Stephen Fegan -- 2015-05-12

   Class providing logging capabilities

*/

#pragma once

#include <string>
#include <streambuf>
#include <iostream>
#include <sstream>
#include <vector>

namespace calin { namespace io { namespace log {

enum Level { FATAL, ERROR, WARNING, INFO, VERBOSE };

struct TimeStamp
{
  std::string string() const;
  static TimeStamp now();
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

  static const char* level_string(Level level);
  static const char* color_string(Level level);
  
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
    (*message_) << t;
    return *this;
  }
 protected:
  Logger* logger_;
  Level level_;
  std::ostringstream* message_;
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

} } } // namespace calin::io::log

