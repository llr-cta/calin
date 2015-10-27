/* 

   calin/io/log.hpp -- Stephen Fegan -- 2015-05-12

   Class providing logging capabilities

*/

#include <Python.h>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include <time.h>

#include "io/log.hpp"

using namespace calin::io::log;

std::string TimeStamp::string() const
{
  time_t ts = time_t(sec);
  struct tm the_tm;
  localtime_r(&ts, &the_tm);
  char buffer[] = "1999-12-31 23:59:59";
  strftime(buffer, sizeof(buffer)-1, "%Y-%m-%dT%H:%M:%S", &the_tm);
  std::string str(buffer);
  uint32_t ms = usec/1000;
  if(ms<10) { str += ".00"; }
  else if(ms<100) { str += ".0"; }
  else { str += "."; }
  str += std::to_string(ms);
  return str;
}

TimeStamp TimeStamp::now()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return { uint64_t(tv.tv_sec), uint32_t(tv.tv_usec) };
}

double TimeStamp::seconds_since(const TimeStamp& then)
{
  double dt { sec>then.sec?double(sec-then.sec):-double(then.sec-sec) };
  dt += (usec>then.usec?double(usec-then.usec):-double(then.usec-usec))*1e-6;
  return dt;
}
  
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
      writer("]",1);
    }      
        
    writer(message.c_str() + spos, n);
    writer("\n",1);
    
    spos += n+1;
  }
}

} // anonymous namespace

MultiLogger::~MultiLogger()
{
  //nothing to see here
}

void MultiLogger::log_message(Level level, const std::string& message,
                              TimeStamp timestamp)
{
  if(message.empty() || level==DISCARD)return;
  
  lock();

  if(sub_loggers_.empty() and sub_streams_.empty())
  {
    if(PythonLogger::is_python_initialised())
      nolock_add_logger(new PythonLogger, true);
    else
      nolock_add_stream(&std::cout, false, false, true);
  }
  
  for(auto& sl : sub_loggers_)
  {
    try
    {
      sl.logger->log_message(level, message, timestamp);
    }
    catch(...)
    {
      // for better or worse - ignore all errors      
    }
  }

  std::string timestamp_string = timestamp.string();
  const char* this_level_string = level_string(level);
  const char* apply_color_string = ansi_color_string(level);
  const char* reset_color_string = ansi_reset_string();
  
  for(auto& ss : sub_streams_)
  {
    try
    {
      std::ostream* stream { ss.stream };
      write_message_lines([stream](const char* c, unsigned n) {
          stream->write(c,n); },
        ss.apply_timestamp?timestamp_string.c_str():nullptr, this_level_string,
        ss.use_colors?apply_color_string:nullptr, reset_color_string,
        message);
    }
    catch(...)
    {
      // for better or worse - ignore all errors
    }
  }
  
  unlock();
}

void MultiLogger::add_logger(Logger* logger, bool adopt_logger)
{
  lock();
  nolock_add_logger(logger, adopt_logger);
  unlock();
}

void MultiLogger::add_stream(std::ostream* stream, bool adopt_stream,
                             bool apply_timestamp, bool use_colors)
{
  lock();
  nolock_add_stream(stream, adopt_stream, apply_timestamp, use_colors);
  unlock();
}

void MultiLogger::nolock_add_logger(Logger* logger, bool adopt_logger)
{
  sub_logger sl { logger, adopt_logger };
  sub_loggers_.push_back(sl);
}

void MultiLogger::nolock_add_stream(std::ostream* stream, bool adopt_stream,
                                      bool apply_timestamp, bool use_colors)
{
  sub_stream ss = { stream, adopt_stream, apply_timestamp, use_colors };
  sub_streams_.push_back(ss);
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

void MultiLogger::lock()
{

}

void MultiLogger::unlock()
{

}

PythonLogger::~PythonLogger()
{
  // nothing to see here
}

void PythonLogger::
log_message(Level level, const std::string& message, TimeStamp timestamp)
{
  if(message.empty() || level==DISCARD)return;
  
  const char* this_level_string = level_string(level);
  const char* apply_color_string = ansi_color_string(level);
  const char* reset_color_string = ansi_reset_string();
  
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
}

bool PythonLogger::is_python_initialised()
{
  return Py_IsInitialized() != 0;
}

