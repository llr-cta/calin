/*

   calin/util/string.hpp -- Stephen Fegan -- 2015-11-13

   Various string functions

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

#include <vector>
#include <string>
#include <cstdint>
#include <sstream>
#include <iomanip>

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

std::string chomp(const std::string& s_in);

std::string string_escape(const std::string& s_in);

std::string to_lower(const std::string& s_in);
std::string to_upper(const std::string& s_in);

std::string reflow(const std::string& s_in,
  unsigned width, const std::string& indent,
  unsigned line1_width, const std::string& line1_indent);

inline std::string reflow(const std::string& s_in,
  unsigned width, unsigned indent,
  unsigned line1_width, unsigned line1_indent)
{
  return reflow(s_in, width, std::string(indent, ' '),
            line1_width, std::string(line1_indent, ' '));
}

inline std::string reflow(const std::string& s_in,
  unsigned width=80, unsigned indent=0)
{
  return reflow(s_in, width, indent, width, indent);
}

// http://stackoverflow.com/questions/7276826/c-format-number-with-commas

template<class T>
std::string to_string_with_commas(T value, unsigned precision = 2)
{
  struct Numpunct: public std::numpunct<char>{
  protected:
    virtual char do_thousands_sep() const{return ',';}
    virtual std::string do_grouping() const{return "\03";}
  };
  std::stringstream ss;
  ss.imbue({std::locale(), new Numpunct});
  ss << std::setprecision(precision) << std::fixed << value;
  return ss.str();
}

bool from_string(const std::string& s, bool& x);

template<typename T> bool from_string(const std::string& s, T& x)
{
  std::istringstream stream(s);
  stream >> std::boolalpha >> x;
  return stream.eof();
}

template<typename T> bool from_string(const std::string& s, std::vector<T>& x)
{
  bool all_good = true;
  std::vector<std::string> bits = split(s, ',');
  for(auto ibit : bits)
  {
    T ix;
    all_good &= from_string(ibit, ix);
    x.emplace_back(ix);
  }
  return all_good;
}

template<typename T> std::string to_string(const T& x)
{
  std::ostringstream stream;
  stream << std::boolalpha << x;
  return stream.str();
}

template<typename T> std::string to_string(const std::vector<T>& x)
{
  std::string s;
  for(const auto& ix : x)
  {
    if(!s.empty())s += ",";
    s += ix;
  }
  return s;
}

#define CALIN_DEFINE_T_FROM_STRING(T, fn) \
inline T fn(const std::string& s, bool* good = nullptr) \
{ \
  T t; \
  bool is_good = from_string(s, t); \
  if(good) *good = is_good; \
  return t; \
}

CALIN_DEFINE_T_FROM_STRING(bool, bool_from_string)
CALIN_DEFINE_T_FROM_STRING(int8_t, int8_from_string)
CALIN_DEFINE_T_FROM_STRING(uint8_t, uint8_from_string)
CALIN_DEFINE_T_FROM_STRING(int16_t, int16_from_string)
CALIN_DEFINE_T_FROM_STRING(uint16_t, uint16_from_string)
CALIN_DEFINE_T_FROM_STRING(int32_t, int32_from_string)
CALIN_DEFINE_T_FROM_STRING(uint32_t, uint32_from_string)
CALIN_DEFINE_T_FROM_STRING(int64_t, int64_from_string)
CALIN_DEFINE_T_FROM_STRING(uint64_t, uint64_from_string)
CALIN_DEFINE_T_FROM_STRING(float, float_from_string)
CALIN_DEFINE_T_FROM_STRING(double, double_from_string)

CALIN_DEFINE_T_FROM_STRING(int, int_from_string)
CALIN_DEFINE_T_FROM_STRING(unsigned, unsigned_from_string)

#undef CALIN_DEFINE_T_FROM_STRING

#define CALIN_DEFINE_T_TO_STRING(T, fn) \
inline std::string fn(const T& val) \
{ \
  return to_string(val); \
}

CALIN_DEFINE_T_TO_STRING(bool, bool_to_string)
CALIN_DEFINE_T_TO_STRING(int8_t, int8_to_string)
CALIN_DEFINE_T_TO_STRING(uint8_t, uint8_to_string)
CALIN_DEFINE_T_TO_STRING(int16_t, int16_to_string)
CALIN_DEFINE_T_TO_STRING(uint16_t, uint16_to_string)
CALIN_DEFINE_T_TO_STRING(int32_t, int32_to_string)
CALIN_DEFINE_T_TO_STRING(uint32_t, uint32_to_string)
CALIN_DEFINE_T_TO_STRING(int64_t, int64_to_string)
CALIN_DEFINE_T_TO_STRING(uint64_t, uint64_to_string)
CALIN_DEFINE_T_TO_STRING(float, float_to_string)
CALIN_DEFINE_T_TO_STRING(double, double_to_string)

CALIN_DEFINE_T_TO_STRING(int, int_to_string)
CALIN_DEFINE_T_TO_STRING(unsigned, unsigned_to_string)

#undef CALIN_DEFINE_T_TO_STRING

#define CALIN_DEFINE_T_TO_STRING_WITH_COMMAS(T, fn) \
inline std::string fn(const T& val, unsigned precision = 2) \
{ \
  return to_string_with_commas(val, precision); \
}

CALIN_DEFINE_T_TO_STRING_WITH_COMMAS(bool, bool_to_string_with_commas)
CALIN_DEFINE_T_TO_STRING_WITH_COMMAS(int8_t, int8_to_string_with_commas)
CALIN_DEFINE_T_TO_STRING_WITH_COMMAS(uint8_t, uint8_to_string_with_commas)
CALIN_DEFINE_T_TO_STRING_WITH_COMMAS(int16_t, int16_to_string_with_commas)
CALIN_DEFINE_T_TO_STRING_WITH_COMMAS(uint16_t, uint16_to_string_with_commas)
CALIN_DEFINE_T_TO_STRING_WITH_COMMAS(int32_t, int32_to_string_with_commas)
CALIN_DEFINE_T_TO_STRING_WITH_COMMAS(uint32_t, uint32_to_string_with_commas)
CALIN_DEFINE_T_TO_STRING_WITH_COMMAS(int64_t, int64_to_string_with_commas)
CALIN_DEFINE_T_TO_STRING_WITH_COMMAS(uint64_t, uint64_to_string_with_commas)
CALIN_DEFINE_T_TO_STRING_WITH_COMMAS(float, float_to_string_with_commas)
CALIN_DEFINE_T_TO_STRING_WITH_COMMAS(double, double_to_string_with_commas)

CALIN_DEFINE_T_TO_STRING_WITH_COMMAS(int, int_to_string_with_commas)
CALIN_DEFINE_T_TO_STRING_WITH_COMMAS(unsigned, unsigned_to_string_with_commas)

#undef CALIN_DEFINE_T_TO_STRING_WITH_COMMAS

// https://stackoverflow.com/questions/5100718/integer-to-hex-string-in-c

template <typename I> std::string number_to_hexstring(I w, size_t hex_len = sizeof(I)<<1)
{
  static const char* digits = "0123456789ABCDEF";
  std::string rc(hex_len,'0');
  for (size_t i=0, j=(hex_len-1)*4 ; i<hex_len; ++i,j-=4)
      rc[i] = digits[(w>>j) & 0x0f];
  return rc;
}

inline std::string instance_identifier(void* ptr)
{
  return number_to_hexstring(reinterpret_cast<intptr_t>(ptr));
}

} } } // namespace calin::util::string
