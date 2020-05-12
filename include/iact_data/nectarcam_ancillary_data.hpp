/*

   calin/iact_data/nectarcam_ancillary_data.hpp -- Stephen Fegan -- 2020-05-11

   Classes to extract NectarCAM ancillary data from SQL database

   Copyright 2020, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#pragma once

#include <string>

#include <io/sqlite3_serializer.hpp>
#include <iact_data/nectarcam_ancillary_data.pb.h>

namespace calin { namespace iact_data { namespace nectarcam_ancillary_data {

#ifndef SWIG
calin::ix::iact_data::nectarcam_ancillary_data::NectarCAMAncillaryData*
retrieve_nectarcam_ancillary_data(const std::string& db_file,
  int camera_id, int64_t start_time_sec, int64_t end_time_sec,
  calin::ix::iact_data::nectarcam_ancillary_data::NectarCAMAncillaryData* data = nullptr,
  bool log_sql = false);

calin::ix::iact_data::nectarcam_ancillary_data::NectarCAMAncillaryData*
retrieve_nectarcam_ancillary_data(const std::string& db_file,
  int camera_id, const std::string& start_time_sec, const std::string& end_time_sec,
  calin::ix::iact_data::nectarcam_ancillary_data::NectarCAMAncillaryData* data = nullptr,
  bool log_sql = false);
#else
void retrieve_nectarcam_ancillary_data(const std::string& db_file,
  int camera_id, int64_t start_time_sec, int64_t end_time_sec,
  calin::ix::iact_data::nectarcam_ancillary_data::NectarCAMAncillaryData* data,
  bool log_sql = false);
calin::ix::iact_data::nectarcam_ancillary_data::NectarCAMAncillaryData*
retrieve_nectarcam_ancillary_data(const std::string& db_file,
  int camera_id, int64_t start_time_sec, int64_t end_time_sec);

void retrieve_nectarcam_ancillary_data(const std::string& db_file,
  int camera_id, const std::string& start_time_sec, const std::string& end_time_sec,
  calin::ix::iact_data::nectarcam_ancillary_data::NectarCAMAncillaryData* data,
  bool log_sql = false);
calin::ix::iact_data::nectarcam_ancillary_data::NectarCAMAncillaryData*
retrieve_nectarcam_ancillary_data(const std::string& db_file,
  int camera_id, const std::string& start_time_sec, const std::string& end_time_sec);
#endif

} } } // namespace calin::iact_data::nectarcam_data_source
