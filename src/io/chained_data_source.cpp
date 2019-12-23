/*

   calin/io/chained_data_source.cpp -- Stephen Fegan -- 2019-12-17

   A supplier of generic data, probably Protobufs, chaining many low level
   DataSources together

   Copyright 2019, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <io/chained_data_source.hpp>

using namespace calin::io::data_source;

FragmentList::~FragmentList()
{
  // nothing to see here
}
