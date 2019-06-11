/*

   calin/util/file.cpp -- Stephen Fegan -- 2016-01-20

   File utility functions, e.g. test whether files exist, expand filenames etc

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

   Based on original, copyright 2006, Stephen Fegan, see notice below

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

/*! \file VSFileUtility.cpp

  Miscellaneous utility file functions

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       11/11/2006
*/

#include <unistd.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <pwd.h>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cstring>

#include <google/protobuf/descriptor.h>
#include <google/protobuf/message.h>
#include <google/protobuf/util/json_util.h>

#include <util/string.hpp>
#include <util/file.hpp>

#include <libgen.h>

// Some libgen.h implementations seem to typedef "basename" to point
// elsewhere.. this undoes this on those systems so calin::util::file::basename
// is not accidentally renamed

namespace {

char* system_basename(char* x) {
  return ::basename(x);
}

char* system_dirname(char* x) {
  return ::dirname(x);
}

} // anonymous namespace

#ifdef basename
#undef basename
#endif

#ifdef dirname
#undef dirname
#endif

using namespace calin::util;

bool calin::util::file::exists(const std::string& filename)
{
  if(filename.empty())return false;
  struct stat statbuf;
  if(stat(filename.c_str(),&statbuf)<0)return false;
  return true;
}

bool calin::util::file::is_file(const std::string& filename)
{
  if(filename.empty())return false;
  struct stat statbuf;
  if(stat(filename.c_str(),&statbuf)<0)return false;
  return S_ISREG(statbuf.st_mode);
}

bool calin::util::file::is_directory(const std::string& filename)
{
  if(filename.empty())return false;
  struct stat statbuf;
  if(stat(filename.c_str(),&statbuf)<0)return false;
  return S_ISDIR(statbuf.st_mode);
}

bool calin::util::file::is_readable(const std::string& filename)
{
  if(filename.empty())return false;
  return access(filename.c_str(),R_OK)==0;
}

bool calin::util::file::is_writable(const std::string& filename)
{
  if(filename.empty())return false;
  return access(filename.c_str(),W_OK)==0;
}

bool calin::util::file::can_write_file(const std::string& filename)
{
  if(filename.empty())return false;

  // Test whether the file exists and is writable OR does not exist
  // but the directory in which it would be created is writable
  if(exists(filename))
    return is_file(filename) and is_writable(filename);
  else
  {
    std::string directory = filename;
    unsigned dpos = 0;
    for(unsigned ipos=0;ipos<filename.length();ipos++)
      if(filename[ipos]=='/')dpos=ipos;
    if((dpos==0)&&(filename[dpos]!='/'))directory=".";
    else directory=directory.substr(0,dpos);
    return is_directory(directory) and is_writable(directory);
  }
  assert(0);
  return false;
}

std::string calin::util::file::dirname(const std::string& filename)
{
  char* cfilename = ::strdup(filename.c_str());
  for(unsigned iend=::strlen(cfilename); iend>0 and cfilename[iend-1]=='/'; iend--)
    cfilename[iend-1]='\0';
  std::string dirname = system_dirname(cfilename);
  ::free(cfilename);
  return dirname;
}

std::string calin::util::file::basename(const std::string& filename, const std::string& suffix)
{
  char* cfilename = ::strdup(filename.c_str());
  for(unsigned iend=::strlen(cfilename); iend>0 and cfilename[iend-1]=='/'; iend--)
    cfilename[iend-1]='\0';
  char* cbasename = system_basename(cfilename);
  if(not suffix.empty()) {
    const char* suffix_cstr = suffix.c_str();
    unsigned cbasenamelen = ::strlen(cbasename);
    unsigned suffixlen = ::strlen(suffix_cstr);
    if(cbasenamelen>suffixlen and
        strcmp(cbasename + cbasenamelen - suffixlen, suffix_cstr)==0)
      cbasename[cbasenamelen - suffixlen] = '\0';
  }
  std::string basename = cbasename;
  ::free(cfilename);
  return basename;
}

std::string calin::util::file::strip_extension(const std::string& filename)
{
  size_t ifind = filename.rfind('.');
  if(ifind == 0 or ifind == std::string::npos)return {};
  return filename.substr(0,ifind);
}

void calin::util::file::expand_filename_in_place(std::string& filename)
{
  /* Do leading tilde expansion, replace with home directory */
  if((!filename.empty())&&(filename[0]=='~'))
  {
    unsigned ipos = filename.find('/',1);
    std::string user = filename.substr(1,ipos-1);
    if(user.empty())
	  {
  	  char* login = getlogin();
  	  if(login)user = login;
  	  else
	    {
	      struct passwd* pwd = getpwuid(getuid());
	      if(pwd)user = pwd->pw_name;
	    }
    }

    if(!user.empty())
	  {
  	  struct passwd* pwd = getpwnam(user.c_str());
  	  if(pwd)filename.replace(0,ipos,pwd->pw_dir);
  	}
  }
  else if(filename.size()>1 and filename[0]=='\\' and filename[1]=='~')
  {
    filename.replace(0,2,"~");
  }
}

unsigned calin::util::file::extract_longest_number_from_filename(const std::string& filename)
{
  std::string::const_iterator longest_num_start = filename.begin();

  for(std::string::const_iterator ichar = filename.begin();
      ichar!=filename.end();ichar++)if(*ichar=='/')longest_num_start=ichar;
  if(*longest_num_start=='/')longest_num_start++;

  std::string::const_iterator longest_num_end = longest_num_start;
  std::string::const_iterator number_begin = longest_num_start;

  for(std::string::const_iterator ichar = longest_num_start;
      ichar!=filename.end(); ichar++)
  {
    if(isdigit(*ichar))
    {
  	  if(number_begin==filename.end())number_begin=ichar;
  	}
    else if(number_begin!=filename.end())
  	{
  	  if((ichar-number_begin)>(longest_num_end-longest_num_start))
  	    longest_num_start = number_begin, longest_num_end = ichar;
  	  number_begin = filename.end();
  	}
  }

  if(number_begin!=filename.end()&&
     ((filename.end()-number_begin)>(longest_num_end-longest_num_start)))
    longest_num_start = number_begin, longest_num_end = filename.end();

  unsigned n = 0;

  if(longest_num_end != longest_num_start)
    n = std::atoi(std::string(longest_num_start,longest_num_end).c_str());

  return n;
}

unsigned calin::util::file::extract_first_number_from_filename(const std::string& filename)
{
  std::string::const_iterator num_start = filename.begin();

  for(std::string::const_iterator ichar = filename.begin();
      ichar!=filename.end();ichar++)if(*ichar=='/')num_start=ichar;
  if(*num_start=='/')num_start++;

  while(num_start!=filename.end() and not isdigit(*num_start))++num_start;
  if(num_start==filename.end())return 0;

  std::string::const_iterator num_end = num_start;
  ++num_end;
  while(num_end!=filename.end() and isdigit(*num_end))++num_end;

  return calin::util::string::unsigned_from_string(std::string(num_start,num_end));
}

namespace {
  std::string::size_type advance_if_not_npos(std::string::size_type a, std::string::size_type b) {
    if(b == std::string::npos)return a;
    else return std::max(a,b);
  }
}

unsigned calin::util::file::extract_run_number_from_filename(const std::string& filename)
{
  std::string::size_type num_start = advance_if_not_npos(0, filename.find_last_of('/'));
  num_start = advance_if_not_npos(num_start, filename.find("Run", num_start));
  num_start = filename.find_first_of("0123456789", num_start);
  if(num_start == std::string::npos)return 0;
  return calin::util::string::unsigned_from_string(filename.substr(num_start));
}

std::vector<std::string> calin::util::file::file_fragments(std::string filename,
  const std::string& extension, unsigned fragment_stride)
{
  expand_filename_in_place(filename);
  if(not is_file(filename))return { };

  std::vector<std::string> filenames;
  fragment_stride = std::max(1U, fragment_stride);
  filenames.emplace_back(filename);

  auto ifind = filename.rfind(extension);
  if(ifind == filename.size()-extension.size())
  {
    unsigned padded_size = 0;
    filename = filename.substr(0, ifind);

    unsigned istart = 0;
    if(not is_file(filename+".1"+extension))
    {
      ifind = filename.rfind('.');
      if(ifind != std::string::npos and
        std::all_of(filename.begin() + ifind + 1, filename.end(), ::isdigit))
      {
        istart = std::stoi(filename.substr(ifind + 1));
        padded_size = filename.size()-ifind-1;
        filename = filename.substr(0, ifind);
      }
    }

    bool fragment_found = true;
    for(unsigned i=istart+fragment_stride; fragment_found; i+=fragment_stride)
    {
      fragment_found = false;
      std::string fragment_i { std::to_string(i) };
      do {
        std::string filename_i { filename+"."+fragment_i+extension };
        if(is_file(filename_i)) {
          filenames.emplace_back(filename_i);
          fragment_found = true;
        } else {
          fragment_i = std::string("0") + fragment_i;
        }
      }while(not fragment_found and fragment_i.size() <= padded_size);
    }
  }
  return filenames;
}

void calin::util::file::replace_question_with_number(std::string& filename, unsigned n)
{
  std::string::size_type iquestion = filename.find('?');
  if(iquestion != filename.npos)
    filename.replace(iquestion,1, std::to_string(n));
}

void calin::util::file::save_protobuf_to_json_file(const std::string& filename,
  const google::protobuf::Message* message)
{
  google::protobuf::util::JsonPrintOptions opt;
  opt.add_whitespace = true;
  opt.always_print_primitive_fields = true;
  std::string s;
  if(!google::protobuf::util::MessageToJsonString(*message, &s, opt).ok())
    throw std::runtime_error("Could not encode message as JSON");
  std::ofstream stream(filename);
  if(!stream)throw std::runtime_error("Could not open file: "+filename);
  stream << s;
  if(!stream)throw std::runtime_error("Error writing file: "+filename);
}

void calin::util::file::load_protobuf_from_json_file(const std::string& filename,
  google::protobuf::Message* message)
{
  std::ifstream stream(filename);
  if(!stream)throw std::runtime_error("Could not open file: "+filename);
  std::stringstream buffer;
  buffer << stream.rdbuf();
  google::protobuf::util::JsonParseOptions opt;
  opt.ignore_unknown_fields = false;
  if(!google::protobuf::util::
    JsonStringToMessage(buffer.str(), message, opt).ok())
    throw std::runtime_error("Could not decode JSON file: "+filename);
}
