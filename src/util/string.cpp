/*

   calin/io/until.cpp -- Stephen Fegan -- 2015-11-13

   Various IO related utility functions

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

#include <sstream>

#include <util/string.hpp>

// Adapted from code copyright 2008 Evan Teran and others under the CC license
// From: http://stackoverflow.com/questions/236129/split-a-string-in-c

std::vector<std::string>& calin::util::string::
split(const std::string &s, char delim, std::vector<std::string> &elems)
{
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) { elems.push_back(item); }
  return elems;
}

std::vector<std::string> calin::util::string::
split(const std::string &s, char delim)
{
  std::vector<std::string> elems;
  split(s, delim, elems);
  return elems;
}

std::string calin::util::string::
join(std::vector<std::string>::const_iterator begin,
     std::vector<std::string>::const_iterator end,
     const std::string& sep)
{
  std::string result;
  auto it = begin;
  if (it != end) {
    result.append(*it);
    for(++it ; it!=end; ++it)result.append(sep).append(*it); }
  return result;
}

std::string calin::util::string::
join(const std::vector<std::string>& vec, const std::string& sep)
{
  return join(vec.cbegin(), vec.cend(), sep);
}

std::string calin::util::string::string_escape(const std::string& s_in)
{
  std::string s_out;
  for(auto c : s_in) {
    if(c=='\t') { s_out += "\\t"; }
    else if(c=='\n') { s_out += "\\n"; }
    else if(c=='\r') { s_out += "\\r"; }
    else if(c=='\r') { s_out += "\\r"; }
    else if(c=='\\') { s_out += "\\\\"; }
    else if(c=='"') { s_out += "\\\""; }
    else s_out += c;
  }
  return s_out;
}
