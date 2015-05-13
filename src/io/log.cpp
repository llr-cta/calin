/* 

   calin/io/log.hpp -- Stephen Fegan -- 2015-05-12

   Class providing logging capabilities

*/

#include <fstream>

#include "io/log.hpp"

using namespace calin::io::log;

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

  std::string timestamp_string = "1999-01-01 12:34:56"; //timestamp.string();
  //const char* level_string = level_string(level);
  const char* apply_color_string = color_string(level);
  const char* reset_color_string = "\x1b[39;0m";

  size_t spos = 0;
  while(spos < message.size())
  {
    size_t epos = message.find_first_of('\n',spos);
    unsigned n = (epos==std::string::npos ? message.size() : epos) - spos;

    for(auto& ss : sub_streams_)
    {
      try
      {
        if(ss.apply_timestamp) (*ss.stream) << timestamp_string;
        if(ss.use_colors and *apply_color_string)
          (*ss.stream) << ": " << apply_color_string;
        else (*ss.stream) << /* ' ' << level_string << */ ": ";
      
        ss.stream->write(message.c_str() + spos, n);
        
        if(ss.use_colors and *apply_color_string)
          (*ss.stream) << reset_color_string;

        ss.stream->write("\n",1);

        spos += n+1;
      }
      catch(...)
      {
        // for better or worse - ignore all errors
      }
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

const char* MultiLogger::color_string(Level level)
{
  switch(level)
  {
    case Level::FATAL:    return "\x1b[31;4m";
    case Level::ERROR:    return "\x1b[35m";
    case Level::WARNING:  return "\x1b[33m";
    case Level::INFO:     return "\x1b[32m";
    case Level::VERBOSE:  return "";
  }
}

const char* MultiLogger::level_string(Level level)
{
  switch(level)
  {
    case Level::FATAL:    return "FATAL";
    case Level::ERROR:    return "ERROR";
    case Level::WARNING:  return "WARN";
    case Level::INFO:     return "INFO";
    case Level::VERBOSE:  return "VERBOSE";
  }
}
