/*

   calin/util/timestamp.hpp -- Stephen Fegan -- 2017-10-29

   Simple (non-precision) timestamp class

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <sys/time.h>
#include <time.h>

#include <util/timestamp.hpp>

using namespace calin::util::timestamp;

std::string Timestamp::as_string(bool utc) const
{
  time_t ts = time_t(unix_sec());
  struct tm the_tm;
  char buffer[80];
  if(utc) {
    gmtime_r(&ts, &the_tm);
    strftime(buffer, sizeof(buffer), "%FT%T.000Z", &the_tm);
  } else {
    localtime_r(&ts, &the_tm);
    strftime(buffer, sizeof(buffer), "%FT%T.000%z", &the_tm);
  }
  std::string str(buffer);
  uint32_t ms = unix_nsec()/1000000;
  if(ms<10) { str.replace(22,1,std::to_string(ms)); }
  else if(ms<100) { str.replace(21,2,std::to_string(ms)); }
  else { str.replace(20,3,std::to_string(ms)); }
  return str;
}

calin::ix::util::timestamp::Timestamp*
Timestamp::as_proto(calin::ix::util::timestamp::Timestamp* x, bool utc) const
{
  if(x == nullptr)x = new calin::ix::util::timestamp::Timestamp;
  x->set_unix_sec(unix_sec());
  x->set_unix_nsec(unix_nsec());
  x->set_nsec(nsec_);
  x->set_printable(as_string(utc));
  return x;
}

calin::ix::util::timestamp::Timestamp*
Timestamp::as_proto(bool utc, calin::ix::util::timestamp::Timestamp* x) const
{
  return this->as_proto(x, utc);
}

Timestamp Timestamp::now()
{
  struct timespec tp;
  clock_gettime(CLOCK_REALTIME, &tp);
  return { tp };
}

double Timestamp::seconds_since(const Timestamp& then) const
{
  double dt = double(nanoseconds_since(then))*1e-9;
  return dt;
}

int64_t Timestamp::nanoseconds_since(const Timestamp& then) const
{
  return nsec_ - then.nsec_;
}
