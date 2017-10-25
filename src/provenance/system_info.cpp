/*

   calin/provenance/system_info.cpp -- Stephen Fegan -- 2016-11-06

   Provenance information about build-time and run-time system environment

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

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

#include <calin_global_config.hpp>
#include <provenance/system_info.hpp>

namespace {

calin::ix::provenance::system_info::BuildInfo* new_build_info()
{
  calin::ix::provenance::system_info::BuildInfo* info =
    new calin::ix::provenance::system_info::BuildInfo;
  info->set_install_prefix(CALIN_BUILD_INSTALL_PREFIX);
  info->set_proto_install_dir(CALIN_PROTO_INSTALL_DIR);
  info->set_proto_header_install_dir(CALIN_PROTO_HEADER_INSTALL_DIR);
  info->set_data_install_dir(CALIN_DATA_INSTALL_DIR);
  info->set_lib_install_dir(CALIN_LIB_INSTALL_DIR);
  info->set_bin_install_dir(CALIN_BIN_INSTALL_DIR);
  info->set_header_install_dir(CALIN_HEADER_INSTALL_DIR);
  info->set_python_install_dir(CALIN_PYTHON_INSTALL_DIR);
  info->set_build_system(CALIN_BUILD_SYSTEM);
  info->set_build_type(CALIN_BUILD_TYPE);
  info->set_build_arch(CALIN_BUILD_ARCH);
  info->set_build_c_compiler_id(CALIN_BUILD_C_COMPILER_ID);
  info->set_build_cxx_compiler_id(CALIN_BUILD_CXX_COMPILER_ID);
  info->set_build_system_fqdn(CALIN_BUILD_SYSTEM_FQDN);
  info->set_build_system_user(CALIN_BUILD_USER);
  info->set_build_date(CALIN_BUILD_DATE);
  info->set_git_origin_url(CALIN_GIT_ORIGIN_URL);
  info->set_git_branch(CALIN_GIT_BRANCH);
  info->set_git_commit_sha1(CALIN_GIT_COMMIT_SHA1);
  info->set_git_commit_date(CALIN_GIT_COMMIT_DATE);
  info->set_git_repo_status(CALIN_GIT_REPO_STATUS);
  return info;
};

calin::ix::provenance::system_info::BuildInfo* singleton_build_info_ =
  new_build_info();

} // anonymous namespace

const calin::ix::provenance::system_info::BuildInfo*
calin::provenance::system_info::build_info()
{
  return singleton_build_info_;
}
