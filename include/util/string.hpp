/*

   calin/util/string.hpp -- Stephen Fegan -- 2015-11-13

   Various string functions

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

#include <vector>
#include <string>

namespace calin { namespace util { namespace string {

// Copyright 2008 Evan Teran and others under the CC license
// From: http://stackoverflow.com/questions/236129/split-a-string-in-c
// Compatible with GPL

std::vector<std::string>&
split(const std::string &s, char delim, std::vector<std::string> &elems);

std::vector<std::string> split(const std::string &s, char delim);

std::string join(std::vector<std::string>::const_iterator begin,
                 std::vector<std::string>::const_iterator end,
                 const std::string& sep);
                 
std::string join(const std::vector<std::string>& vec, const std::string& sep);

} } } // namespace calin::util::string
