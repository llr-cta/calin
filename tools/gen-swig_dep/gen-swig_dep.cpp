/*

   calin/tools/gen-swig_dep.cpp -- Stephen Fegan -- 2016-11-18

   Parser for generating SWIG dependencies interface files

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <pwd.h>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cerrno>

using namespace std;

namespace {

string chomp(const string& s_in)
{
  string::size_type ifind = s_in.find_first_not_of(" \t");
  if(ifind == string::npos) {
    return string();
  } else {
    return s_in.substr(ifind);
  }
}

string gen_filename(string swig_filename,
  const string& srcdir, const string& bindir, string extension = ".swig_dep")
{
  if(swig_filename.size()>=srcdir.size() and
      swig_filename.substr(0,srcdir.size()) == srcdir)
    swig_filename = bindir + swig_filename.substr(srcdir.size());
  else if(swig_filename.empty() or swig_filename[0] != '/')
    swig_filename = bindir + "/" + swig_filename;
  return swig_filename.substr(0,swig_filename.find_last_of('.')) + extension;
}

void write_dep(std::ostream& depstream, const string& dep, string& sep,
  const string & src_dir, const string & bin_dir)
{
  static const string PBI = ".pb.i";
  static const string PROTO = ".proto";
  string cmake_sep = ";";

  string test_file = src_dir + "/swig/" + dep;
  if(access(test_file.c_str(),R_OK)==0) {
    depstream << sep << test_file;
    sep = cmake_sep;
    return;
  }

  test_file = src_dir + "/include/" + dep;
  if(access(test_file.c_str(),R_OK)==0) {
    depstream << sep << test_file;
    sep = cmake_sep;
    return;
  }

  if(dep.size()>PBI.size() and
    dep.substr(dep.size()-PBI.size(), PBI.size()) == PBI)
  {
    string basename = dep.substr(0, dep.size()-PBI.size());
    test_file = src_dir + "/proto/" + basename + PROTO;
    if(access(test_file.c_str(),R_OK)==0) {
      //      depstream << sep << bin_dir + "/proto/" + dep;
      if(basename.find_last_of('/') != std::string::npos) {
        depstream << sep << "protoc-gen-swig-" << basename.substr(basename.find_last_of('/')+1);
      } else {
        depstream << sep << "protoc-gen-swig-" << basename;
      }
      sep = cmake_sep;
      return;
    }
  }
}


} // anonymous namespace

void usage(const string& progname)
{
  std::cerr << "Usage: " << progname << " swig_file src_dir bin_dir\n";
  exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
{
  string progname(*argv);
  argv++, argc--;

  if(argc != 3)usage(progname);

  string swigname(*argv);
  argv++, argc--;

  string srcdir(*argv);
  argv++, argc--;

  string bindir(*argv);
  argv++, argc--;

  std::ifstream stream(swigname.c_str());
  if(!stream.good())
  {
    std::cerr <<
      "Could not open file " << swigname << ": " << strerror(errno) << '\n';
    exit(EXIT_FAILURE);
  }

  std::ofstream depstream(gen_filename(swigname, srcdir, bindir).c_str());

  string line;
  string sep;
  while(std::getline(stream, line))
  {
    line = chomp(line);
    if(line.size()>=7 and line.substr(0,7) == "%import")
      line = chomp(line.substr(7));
    else if(line.size()>=8 and line.substr(0,8) == "%include")
      line = chomp(line.substr(8));
    else continue;
    if(line.empty())continue;
    size_t closing_delim = string::npos;
    if(line[0]=='"')closing_delim = line.find_first_of('"', 1);
    else if(line[0]=='<')closing_delim = line.find_first_of('>', 1);
    if(closing_delim != string::npos) {
      write_dep(depstream, line.substr(1,closing_delim-1), sep, srcdir, bindir);
    }
  }
}
