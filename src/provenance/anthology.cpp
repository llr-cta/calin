/*

   calin/provenance/anthology.cpp -- Stephen Fegan -- 2017-10-30

   Anthology of all provenance system info

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <util/timestamp.hpp>
#include <provenance/anthology.hpp>
#include <provenance/chronicle.hpp>
#include <provenance/system_info.hpp>
#include <util/log.hpp>

using namespace calin::provenance::anthology;

calin::ix::provenance::anthology::Anthology*
calin::provenance::anthology::get_current_anthology(calin::ix::provenance::anthology::Anthology* x)
{
  if(x == nullptr) x = new calin::ix::provenance::anthology::Anthology;
  calin::util::timestamp::Timestamp::now().as_proto(x->mutable_timestamp());
  x->mutable_default_log()->CopyFrom(calin::util::log::default_protobuf_logger()->log_messages());
  calin::provenance::system_info::copy_the_build_info(x->mutable_build_info());
  calin::provenance::system_info::copy_the_host_info(x->mutable_host_info());
  calin::provenance::chronicle::copy_the_chronicle(x->mutable_chronicle());
  return x;
}
