/*

   calin/provenance/system_info.cpp -- Stephen Fegan -- 2016-11-06

   Provenance information about build-time and run-time system environment

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

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

#include <unistd.h>
#include <thread>
#include <sys/utsname.h>
#include <cpuid.h>

#include <calin_global_config.hpp>
#include <provenance/system_info.hpp>
#include <io/log.hpp>

extern char **environ;

#include "cpuid_bits.h"

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

bool append_4chars(std::string& t, unsigned u)
{
  for(unsigned i=0; i<4; ++i) {
    char c = u & 0xFF;
    if(c == '\0')return false;
    t += c;
    u >>= 8;
  }
  return true;
}

calin::ix::provenance::system_info::HostAndProcessInfo* new_host_info()
{
  calin::ix::provenance::system_info::HostAndProcessInfo* info =
    new calin::ix::provenance::system_info::HostAndProcessInfo;
  info->set_process_id(::getpid());
  info->set_user_id(::getuid());
  if(::getenv("USER"))info->set_user_name(::getenv("USER"));
  calin::io::log::TimeStamp ts = calin::io::log::TimeStamp::now();
  info->set_process_start_time_unix_sec(ts.sec);
  info->set_process_start_time_unix_usec(ts.usec);
  info->set_process_start_time(ts.string());
  // char hostname_buffer[256];
  // gethostname(hostname_buffer, 255);
  // info->set_host_name(hostname_buffer);
  info->set_hardware_concurrency(std::thread::hardware_concurrency());

  for(char** ienv = environ; *ienv; ++ienv) {
    std::string env(*ienv);
    auto iequals = env.find_first_of('=');
    if(iequals == std::string::npos)
      (*info->mutable_environment())[env] = "";
    else if(iequals == env.size()-1)
      (*info->mutable_environment())[env.substr(0,iequals)] = "";
    else
      (*info->mutable_environment())[env.substr(0,iequals)] = env.substr(iequals+1);
  }

  struct utsname uname_info;
  ::uname(&uname_info);
  info->set_uname_sysname(uname_info.sysname);
  // info->set_uname_nodename(uname_info.nodename);
  info->set_host_name(uname_info.nodename);
  info->set_uname_release(uname_info.release);
  info->set_uname_version(uname_info.version);
  info->set_uname_machine(uname_info.machine);

  unsigned a,b,c,d;
  __cpuid (0 /* vendor string */, a, b, c, d);

  unsigned max_frame = a;
  std::string vendor_string_12char;
  append_4chars(vendor_string_12char, b) &&
  append_4chars(vendor_string_12char, d) &&
  append_4chars(vendor_string_12char, c);
  info->set_cpu_vendor_string(vendor_string_12char);

  __cpuid (0x80000000U /* get highest extended function support */, a, b, c, d);
  unsigned max_eframe = a;
  if(max_eframe >= 0x80000004U)
  {
    std::string processor_brand_48char;
    bool more = true;
    for(unsigned iframe = 0x80000002; more && iframe<=0x80000004; iframe++)
    {
      __cpuid (iframe /* processor brand string */ , a, b, c, d);
      more = append_4chars(processor_brand_48char, a) &&
        append_4chars(processor_brand_48char, b) &&
        append_4chars(processor_brand_48char, c) &&
        append_4chars(processor_brand_48char, d);
    }
    info->set_cpu_processor_brand(processor_brand_48char);
  }

  if(max_frame >= 1)
  {
    __cpuid (1 /* Processor Info and Feature Bits */, a, b, c, d);

    info->set_cpu_model(((a>>4)&0x0F) + ((a>>12)&0xF0));
    info->set_cpu_family(((a>>8)&0x0F) + ((a>>20)&0x0F));
    info->set_cpu_stepping(a & 0x0F);
    info->set_cpu_has_fpu(d & 1);
    info->set_cpu_has_mmx(d & bit_MMX);
    info->set_cpu_has_sse(d & bit_SSE);
    info->set_cpu_has_sse2(d & bit_SSE2);
    info->set_cpu_has_sse3(c & bit_SSE3);
    info->set_cpu_has_ssse3(c & bit_SSSE3);
    info->set_cpu_has_sse4_1(c & bit_SSE4_1);
    info->set_cpu_has_sse4_2(c & bit_SSE4_2);
    info->set_cpu_has_pclmulqdq(c & bit_PCLMUL);
    info->set_cpu_has_avx(c & bit_AVX);
    info->set_cpu_has_fma3(c & bit_FMA);
  }

  if(max_frame >= 7)
  {
    __cpuid_count (7, 0 /* Extended Features */, a, b, c, d);

    info->set_cpu_has_avx2(b & bit_AVX2);
    info->set_cpu_has_bmi1(b & bit_BMI);
    info->set_cpu_has_bmi2(b & bit_BMI2);
    info->set_cpu_has_adx(b & bit_ADX);
    info->set_cpu_has_avx512f(b & bit_AVX512F);
    info->set_cpu_has_avx512dq(b & bit_AVX512DQ);
    info->set_cpu_has_avx512ifma(b & bit_AVX512IFMA);
    info->set_cpu_has_avx512pf(b & bit_AVX512PF);
    info->set_cpu_has_avx512er(b & bit_AVX512ER);
    info->set_cpu_has_avx512cd(b & bit_AVX512CD);
    info->set_cpu_has_avx512bw(b & bit_AVX512BW);
    info->set_cpu_has_avx512vl(b & bit_AVX512VL);

    info->set_cpu_has_avx512vbmi(c & bit_AVX512VBMI);
    info->set_cpu_has_avx512vbmi2(c & bit_AVX512VBMI2);
    info->set_cpu_has_avx512vnni(c & bit_AVX512VNNI);
    info->set_cpu_has_avx512bitalg(c & bit_AVX512BITALG);
    info->set_cpu_has_avx512vpopcntdq(c & bit_AVX512VPOPCNTDQ);

    info->set_cpu_has_avx512_4vnniw(d & bit_AVX512_4VNNIW);
    info->set_cpu_has_avx512_4fmaps(d & bit_AVX512_4FMAPS);
  }

  if(max_eframe >= 0x80000001)
  {
    __cpuid (0x80000001 /* Extended Processor Info and Feature Bits */, a, b, c, d);

    info->set_cpu_has_fma4(c & bit_FMA4);
  }

#if 0

bool cpu_has_fma4                                        = 222 [
  (CFO).desc = "CPU indicates that it has FMA-4 (cpuid)." ];

#endif

  return info;
}

calin::ix::provenance::system_info::HostAndProcessInfo* singleton_host_info_ =
  new_host_info();

} // anonymous namespace

const calin::ix::provenance::system_info::BuildInfo*
calin::provenance::system_info::build_info()
{
  return singleton_build_info_;
}

const calin::ix::provenance::system_info::HostAndProcessInfo*
calin::provenance::system_info::host_info()
{
  return singleton_host_info_;
}
