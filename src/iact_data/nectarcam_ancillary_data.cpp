/*

   calin/iact_data/nectarcam_ancillary_data.cpp -- Stephen Fegan -- 2020-05-11

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

#include <util/log.hpp>
#include <io/sqlite3_serializer.hpp>
#include <iact_data/nectarcam_ancillary_data.hpp>

using namespace calin::iact_data::nectarcam_ancillary_data;
using namespace calin::util::log;

namespace {
  calin::ix::iact_data::nectarcam_ancillary_data::NectarCAMAncillaryData*
  do_retrieve_nectarcam_ancillary_data(const std::string db_file,
    const google::protobuf::Message* start, const google::protobuf::Message* end,
    calin::ix::iact_data::nectarcam_ancillary_data::NectarCAMAncillaryData* data,
    bool log_sql)
  {
    if(data == nullptr) {
      data = new calin::ix::iact_data::nectarcam_ancillary_data::NectarCAMAncillaryData();
    } else {
      data->Clear();
    }

    auto* sql = new calin::io::sql_serializer::SQLite3Serializer(db_file,
      calin::io::sql_serializer::SQLite3Serializer::READ_ONLY_NON_CALIN_DB);
    sql->set_write_sql_to_log(log_sql);
    sql->set_write_errors_to_log(log_sql);

    data->set_database(db_file);

    sql->register_externally_created_table("monitoring_drawer_temperatures",
      calin::ix::iact_data::nectarcam_ancillary_data::FEBTemperatureMeasurement::descriptor());
    sql->register_externally_created_table("monitoring_channel_currents",
      calin::ix::iact_data::nectarcam_ancillary_data::HVPACurrentMeasurement::descriptor());
    sql->register_externally_created_table("monitoring_channel_voltages",
      calin::ix::iact_data::nectarcam_ancillary_data::HVPAVoltageMeasurement::descriptor());
    sql->register_externally_created_table("monitoring_ecc",
      calin::ix::iact_data::nectarcam_ancillary_data::ECCMeasurement::descriptor());

    std::vector<uint64_t> oids;
    bool good = true;

    try {
      oids = sql->select_oids_in_range("monitoring_drawer_temperatures", start, end);
      for(auto oid : oids) {
        calin::ix::iact_data::nectarcam_ancillary_data::FEBTemperatureMeasurement meas;
        meas.Clear();
        if(sql->retrieve_by_oid("monitoring_drawer_temperatures", oid, &meas)) {
          data->increment_num_feb_temperature_measurements();
          (*data->mutable_feb_temperature())[meas.drawer()].add_measurement()->CopyFrom(meas);
        } else {
          good = false;
        }
      }
    } catch(...) {
      good = false;
    }
    if(not good) {
      LOG(WARNING) << "Error retrieving FEB temperature measurements from ancillary DB";
    }

    good = true;
    try {
      oids = sql->select_oids_in_range("monitoring_channel_currents", start, end);
      for(auto oid : oids) {
        calin::ix::iact_data::nectarcam_ancillary_data::HVPACurrentMeasurement meas;
        meas.Clear();
        if(sql->retrieve_by_oid("monitoring_channel_currents", oid, &meas)) {
          data->increment_num_hvpa_current_measurements();
          (*data->mutable_hvpa_current())[meas.drawer()*7 + meas.channel()].add_measurement()->CopyFrom(meas);
        } else {
          good = false;
        }
      }
    } catch(...) {
      good = false;
    }
    if(not good) {
      LOG(WARNING) << "Error retrieving FPM current measurements from ancillary DB";
    }

    good = true;
    try {
      oids = sql->select_oids_in_range("monitoring_channel_voltages", start, end);
      for(auto oid : oids) {
        calin::ix::iact_data::nectarcam_ancillary_data::HVPAVoltageMeasurement meas;
        meas.Clear();
        if(sql->retrieve_by_oid("monitoring_channel_voltages", oid, &meas)) {
          data->increment_num_hvpa_voltage_measurements();
          (*data->mutable_hvpa_voltage())[meas.drawer()*7 + meas.channel()].add_measurement()->CopyFrom(meas);
        } else {
          good = false;
        }
      }
    } catch(...) {
      good = false;
    }
    if(not good) {
      LOG(WARNING) << "Error retrieving FPM voltage measurements from ancillary DB";
    }

    good = true;
    try {
      oids = sql->select_oids_in_range("monitoring_ecc", start, end);
      for(auto oid : oids) {
        calin::ix::iact_data::nectarcam_ancillary_data::ECCMeasurement meas;
        meas.Clear();
        if(sql->retrieve_by_oid("monitoring_ecc", oid, &meas)) {
          data->increment_num_ecc_measurements();
          data->mutable_ecc_measurements()->add_measurement()->CopyFrom(meas);
        } else {
          good = false;
        }
      }
      data->set_ecc_data_version(1);
    } catch(...) {
      good = false;
    }
    if(not good) {
      good = true;
      try {
        oids = sql->select_oids_in_range("monitoring_ecc", start, end);
        for(auto oid : oids) {
          calin::ix::iact_data::nectarcam_ancillary_data::ECCMeasurement meas;
          calin::ix::iact_data::nectarcam_ancillary_data::ECCMeasurement_OLD old_meas;
          meas.Clear();
          old_meas.Clear();
          if(sql->retrieve_by_oid("monitoring_ecc", oid, &old_meas)) {
            data->increment_num_ecc_measurements();
#define COPY_OLD_2_NEW(x) meas.set_##x(old_meas.x())
            COPY_OLD_2_NEW(camera);
            COPY_OLD_2_NEW(time);
            COPY_OLD_2_NEW(state);
            COPY_OLD_2_NEW(temp_avg);
          	COPY_OLD_2_NEW(temp_01);
          	COPY_OLD_2_NEW(temp_02);
          	COPY_OLD_2_NEW(temp_03);
          	COPY_OLD_2_NEW(temp_04);
          	COPY_OLD_2_NEW(temp_05);
          	COPY_OLD_2_NEW(temp_06);
          	COPY_OLD_2_NEW(temp_09);
          	COPY_OLD_2_NEW(temp_10);
          	COPY_OLD_2_NEW(temp_11);
          	COPY_OLD_2_NEW(temp_12);
          	COPY_OLD_2_NEW(current_psb1);
          	COPY_OLD_2_NEW(current_psb2);
            COPY_OLD_2_NEW(fan_speed_03);
          	COPY_OLD_2_NEW(fan_speed_04);
          	COPY_OLD_2_NEW(air_pressure_diff);
          	COPY_OLD_2_NEW(water_pressure_01);
          	COPY_OLD_2_NEW(water_pressure_02);
          	COPY_OLD_2_NEW(humidity);
          	COPY_OLD_2_NEW(back_doors_open);
          	COPY_OLD_2_NEW(shutter_status);
          	COPY_OLD_2_NEW(light_01);
          	COPY_OLD_2_NEW(light_02);
          	COPY_OLD_2_NEW(light_03);
          	COPY_OLD_2_NEW(light_04);
          } else {
            good = false;
          }
        }
      } catch(...) {
        good = false;
      }
    }
    if(not good) {
      LOG(WARNING) << "Error retrieving ECC measurements from ancillary DB";
    }

    delete sql;
    return data;
  }
}

calin::ix::iact_data::nectarcam_ancillary_data::NectarCAMAncillaryData*
calin::iact_data::nectarcam_ancillary_data::
retrieve_nectarcam_ancillary_data(const std::string& db_file,
  int camera_id, int64_t start_time_sec, int64_t end_time_sec,
  calin::ix::iact_data::nectarcam_ancillary_data::NectarCAMAncillaryData* data,
  bool log_sql)
{
  calin::ix::iact_data::nectarcam_ancillary_data::SelectByTimeAndCamera start;
  start.set_camera(camera_id);
  start.set_time(start_time_sec);
  calin::ix::iact_data::nectarcam_ancillary_data::SelectByTimeAndCamera end;
  end.set_camera(camera_id);
  end.set_time(end_time_sec);
  return do_retrieve_nectarcam_ancillary_data(db_file, &start, &end, data, log_sql);
}

calin::ix::iact_data::nectarcam_ancillary_data::NectarCAMAncillaryData*
calin::iact_data::nectarcam_ancillary_data::
retrieve_nectarcam_ancillary_data(const std::string& db_file,
  int camera_id, const std::string& start_time_sec, const std::string& end_time_sec,
  calin::ix::iact_data::nectarcam_ancillary_data::NectarCAMAncillaryData* data,
  bool log_sql)
{
  calin::ix::iact_data::nectarcam_ancillary_data::SelectByTimeStringAndCamera start;
  start.set_camera(camera_id);
  start.set_time(start_time_sec);
  calin::ix::iact_data::nectarcam_ancillary_data::SelectByTimeStringAndCamera end;
  end.set_camera(camera_id);
  end.set_time(end_time_sec);
  return do_retrieve_nectarcam_ancillary_data(db_file, &start, &end, data, log_sql);
}
