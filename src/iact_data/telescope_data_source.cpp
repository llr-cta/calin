/*

   calin/iact_data/telescope_data_source.cpp -- Stephen Fegan -- 2016-01-08

   A supplier of single telescope data, for example:
   ix::iact_data::telescope_event::TelescopeEvent

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

//#define CALIN_TELESCOPE_DATA_SOURCE_NO_EXTERN
#include <iact_data/telescope_data_source.hpp>

using namespace calin::util::log;
using namespace calin::iact_data::telescope_data_source;

template class calin::io::data_source::DataSource<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;
template class calin::io::data_source::RandomAccessDataSource<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;
template class calin::io::data_source::ProtobufFileDataSource<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;
template class calin::io::data_source::BufferedDataSource<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;
template class calin::io::data_source::UnidirectionalBufferedDataSourcePump<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;

template class calin::io::data_source::DataSink<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;
template class calin::io::data_source::ProtobufFileDataSink<
  calin::ix::iact_data::telescope_event::TelescopeEvent>;

TelescopeDataSourceWithRunConfig::
~TelescopeDataSourceWithRunConfig()
{
  // nothing to see here
}

TelescopeRandomAccessDataSourceWithRunConfig::
~TelescopeRandomAccessDataSourceWithRunConfig()
{
  // nothing to see here
}

void calin::iact_data::telescope_data_source::
report_run_configuration_problems(
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
  calin::util::log::Logger* logger)
{
  unsigned cam_nmod = run_config->camera_layout().module_size();
  if(run_config->configured_module_id_size())
  {
    bool module_id_strictly_increasing = true;
    std::map<unsigned, unsigned> mod_map;
    unsigned mod_id = run_config->configured_module_id(0);
    if(mod_id >= cam_nmod)
      LOG(ERROR,logger) << "Expected module #" << mod_id << " not valid in camera";
    else
      mod_map[mod_id] = 0;
    for(int irank=1; irank<run_config->configured_module_id_size();irank++)
    {
      unsigned next_mod_id = run_config->configured_module_id(irank);
      if(next_mod_id <= mod_id)module_id_strictly_increasing=false;
      mod_id = next_mod_id;
      if(mod_id >= cam_nmod)
        LOG(ERROR,logger) << "Expected module #" << mod_id << " not valid in camera";
      else if(mod_map.find(mod_id) == mod_map.end()) {
        mod_map[mod_id] = irank;
      } else {
        LOG(WARNING,logger) << "Module #" << mod_id << " has multiple expected ranks: "
          << irank  << " and " << mod_map[mod_id];
      }
    }
    // if(!module_id_strictly_increasing)
    //   LOG(INFO,logger) << "Module ids not strictly increasing.";
    if(mod_map.size() != cam_nmod) {
      auto log = LOG(INFO,logger);
      log << mod_map.size() << " modules configured, " << cam_nmod-mod_map.size()
        << " missing:";
      unsigned imod_print = 0;
      for(unsigned imod=0;imod<cam_nmod;imod++) {
        if(mod_map.find(imod) == mod_map.end()) {
          if(imod_print == 10) {
            log << " ...";
            break;
          } else {
            log << ' ' << imod;
            imod_print++;
          }
        }
      }
    }
  }
}
