/* 

   calin/io/log.hpp -- Stephen Fegan -- 2015-05-12

   Class providing logging capabilities

*/

#include <fstream>
#include <Python.h>

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

Logger::~Logger()
{
  //nothing to see here
}

MultiLogger::~MultiLogger()
{
  //nothing to see here
}

void MultiLogger::log_message(Level level, const std::string& message,
                              TimeStamp timestamp)
{
  if(message.empty())return;
  
  lock();

  if(sub_loggers_.empty() and sub_streams_.empty())
  {
    if(PythonLogger::is_python_initialised())
      nolock_add_logger(new PythonLogger, true);
    else
      nolock_add_stream(&std::cout, false, false, false);
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
  const char* apply_color_string = color_string(level);
  const char* reset_color_string = "\x1b[0m";

  size_t spos = 0;
  while(spos < message.size())
  {
    size_t epos = message.find_first_of('\n',spos);
    unsigned n = (epos==std::string::npos ? message.size() : epos) - spos;

    for(auto& ss : sub_streams_)
    {
      try
      {
        if(ss.apply_timestamp) (*ss.stream) << timestamp_string << ' ';
        if(ss.use_colors and *apply_color_string and *this_level_string)
        {
          (*ss.stream) << apply_color_string << ' ' << this_level_string
                       << ' ' << reset_color_string << ' ';
        }
        else if(*this_level_string)
        {
          (*ss.stream) << '[' << this_level_string << ']' << ' ';
        }
        
        ss.stream->write(message.c_str() + spos, n);
        ss.stream->write("\n",1);
      }
      catch(...)
      {
        // for better or worse - ignore all errors
      }
    }

    spos += n+1;    
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

const char* MultiLogger::color_string(Level level)
{
  switch(level)
  {
    case Level::FATAL:    return "\x1b[37;41;97;101;1m";
    case Level::ERROR:    return "\x1b[30;45;105;1m";
    case Level::WARNING:  return "\x1b[30;43;103;1m";
    case Level::INFO:     return "\x1b[37;46;97;1m";
    case Level::VERBOSE:  return "";
  }
}

const char* MultiLogger::level_string(Level level)
{
  switch(level)
  {
    case Level::FATAL:    return "FATAL";
    case Level::ERROR:    return "ERROR";
    case Level::WARNING:  return "WARNING";
    case Level::INFO:     return "INFO";
    case Level::VERBOSE:  return "";
  }
}

PythonLogger::~PythonLogger()
{
  // nothing to see here
}

void PythonLogger::
log_message(Level level, const std::string& message, TimeStamp timestamp)
{
  if(use_stderr_)
    PySys_WriteStderr("%.1000s", message.c_str());
  else
    PySys_WriteStdout("%.1000s", message.c_str());
}

bool PythonLogger::is_python_initialised()
{
  return Py_IsInitialized() != 0;
}

