/*

   calin/provenance/system_info.cpp -- Stephen Fegan -- 2016-11-06

   Provenance information about build-time and run-time system environment

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

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
#include <sys/utsname.h>
#include <cpuid.h>

#include <thread>
#include <memory>
#include <time.h>

#include <calin_global_config.hpp>
#include <provenance/system_info.hpp>
#include <util/log.hpp>
#include <util/timestamp.hpp>
#include <util/vcl.hpp>

extern char **environ;

#include <provenance/cpuid_bits.h>

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
  info->set_unit_test_install_dir(CALIN_UNIT_TEST_INSTALL_DIR);
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

#ifdef __x86_64__
  info->set_compiled_with_x86_64(true);
#endif
#ifdef __SSE__
  info->set_compiled_with_sse(true);
#endif
#ifdef __SSE2__
  info->set_compiled_with_sse2(true);
#endif
#ifdef __SSE3__
  info->set_compiled_with_sse3(true);
#endif
#ifdef __SSSE3__
  info->set_compiled_with_ssse3(true);
#endif
#ifdef __SSE4_1__
  info->set_compiled_with_sse4_1(true);
#endif
#ifdef __SSE4_2__
  info->set_compiled_with_sse4_2(true);
#endif
#ifdef __AVX__
  info->set_compiled_with_avx(true);
#endif
#ifdef __AVX2__
  info->set_compiled_with_avx2(true);
#endif
#ifdef __FMA__
  info->set_compiled_with_fma3(true);
#endif
#ifdef __FMA4__
  info->set_compiled_with_fma4(true);
#endif
#ifdef __AVX512__
  info->set_compiled_with_avx512(true);
#endif
#ifdef __AVX512F__
  info->set_compiled_with_avx512f(true);
#endif
#ifdef __AVX512DQ__
  info->set_compiled_with_avx512dq(true);
#endif
#ifdef __AVX512IFMA__
  info->set_compiled_with_avx512ifma(true);
#endif
#ifdef __AVX512PF__
  info->set_compiled_with_avx512pf(true);
#endif
#ifdef __AVX512ER__
  info->set_compiled_with_avx512er(true);
#endif
#ifdef __AVX512CD__
  info->set_compiled_with_avx512cd(true);
#endif
#ifdef __AVX512BW__
  info->set_compiled_with_avx512bw(true);
#endif
#ifdef __AVX512VL__
  info->set_compiled_with_avx512vl(true);
#endif
#ifdef __AVX512VBMI__
  info->set_compiled_with_avx512vbmi(true);
#endif
#ifdef __AVX512VBMI2__
  info->set_compiled_with_avx512vbmi2(true);
#endif
#ifdef __AVX512VNNI__
  info->set_compiled_with_avx512vnni(true);
#endif
#ifdef __AVX512BITALG__
  info->set_compiled_with_avx512bitalg(true);
#endif
#ifdef __AVX512VPOPCNTDQ__
  info->set_compiled_with_avx512vpopcntdq(true);
#endif
#ifdef __AVX5124VNNIW__
  info->set_compiled_with_avx512_4vnniw(true);
#endif
#ifdef __AVX5124MAPS__
  info->set_compiled_with_avx512_4maps(true);
#endif

  return info;
};

std::unique_ptr<calin::ix::provenance::system_info::BuildInfo> singleton_build_info_ {
  new_build_info() };

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
  std::set<std::string> env_whitelist;
  env_whitelist.insert("USER");
  env_whitelist.insert("PATH");
  env_whitelist.insert("LD_LIBRARY_PATH");
  env_whitelist.insert("PYTHONPATH");
  env_whitelist.insert("DYLD_LIBRARY_PATH");
  env_whitelist.insert("G4DATADIR");
  env_whitelist.insert("G4DATADIR");
  env_whitelist.insert("G4NEUTRONHPDATA");
  env_whitelist.insert("G4LEDATA");
  env_whitelist.insert("G4LEVELGAMMADATA");
  env_whitelist.insert("G4RADIOACTIVEDATA");
  env_whitelist.insert("G4SAIDXSDATA");
  env_whitelist.insert("G4PARTICLEXSDATA");
  env_whitelist.insert("G4ABLADATA");
  env_whitelist.insert("G4INCLDATA");
  env_whitelist.insert("G4PIIDATA");
  env_whitelist.insert("G4ENSDFSTATEDATA");
  env_whitelist.insert("G4REALSURFACEDATA");
  env_whitelist.insert("G4TENDL");

  calin::ix::provenance::system_info::HostAndProcessInfo* info =
    new calin::ix::provenance::system_info::HostAndProcessInfo;
  info->set_process_id(::getpid());
  info->set_user_id(::getuid());
  if(::getenv("USER"))info->set_user_name(::getenv("USER"));
  calin::util::timestamp::Timestamp ts = calin::util::timestamp::Timestamp::now();
  ts.as_proto(info->mutable_process_start_time());
  // char hostname_buffer[256];
  // gethostname(hostname_buffer, 255);
  // info->set_host_name(hostname_buffer);
  info->set_hardware_concurrency(std::thread::hardware_concurrency());

  for(char** ienv = environ; *ienv; ++ienv) {
    std::string env(*ienv);
    auto iequals = env.find_first_of('=');
    std::string evar;
    std::string eval;
    if(iequals == std::string::npos) {
      evar = env; eval = "";
    } else if(iequals == env.size()-1) {
      evar = env.substr(0,iequals); eval = "";
    } else {
      evar = env.substr(0,iequals); eval = env.substr(iequals+1);
    }
    if(env_whitelist.count(evar))
      (*info->mutable_environment())[evar] = eval;
  }

  char* cwd = ::getcwd(NULL, 0);
  info->set_current_working_directory(cwd);
  ::free(cwd);

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

  unsigned simd_vec_size = 0;
  if(info->cpu_has_avx512f())simd_vec_size=512/8;
  else if(info->cpu_has_avx())simd_vec_size=256/8;
  else if(info->cpu_has_sse2())simd_vec_size=128/8;
  else simd_vec_size = 64/8;
  info->set_simd_vec_size(simd_vec_size);
  unsigned log2_simd_vec_size = 0;
  while(simd_vec_size>>=1)++log2_simd_vec_size;
  info->set_log2_simd_vec_size(log2_simd_vec_size);

  return info;
}

std::unique_ptr<calin::ix::provenance::system_info::HostAndProcessInfo> singleton_host_info_ {
  new_host_info() };

} // anonymous namespace

const calin::ix::provenance::system_info::BuildInfo*
calin::provenance::system_info::the_build_info()
{
  return singleton_build_info_.get();
}

const calin::ix::provenance::system_info::HostAndProcessInfo*
calin::provenance::system_info::the_host_info()
{
  return singleton_host_info_.get();
}

calin::ix::provenance::system_info::BuildInfo*
calin::provenance::system_info::copy_the_build_info(calin::ix::provenance::system_info::BuildInfo* x)
{
  if(x == nullptr) x = new calin::ix::provenance::system_info::BuildInfo;
  x->CopyFrom(*calin::provenance::system_info::the_build_info());
  return x;
}

calin::ix::provenance::system_info::HostAndProcessInfo*
calin::provenance::system_info::copy_the_host_info(calin::ix::provenance::system_info::HostAndProcessInfo* x)
{
  if(x == nullptr) x = new calin::ix::provenance::system_info::HostAndProcessInfo;
  x->CopyFrom(*calin::provenance::system_info::the_host_info());
  return x;
}

std::string calin::provenance::system_info::build_info_string(
  const calin::ix::provenance::system_info::BuildInfo* build_info)
{
  if(build_info == nullptr) {
    build_info = calin::provenance::system_info::the_build_info();
  }
  std::ostringstream L;
  L << "Git origin : " << build_info->git_origin_url();
  if(not build_info->git_branch().empty()) {
    L << "  Branch : " << build_info->git_branch();
  }
  L << '\n'
    << "Git commit : " << build_info->git_commit_sha1() << " from " << build_info->git_commit_date();
  if(build_info->git_repo_status() == "dirty") {
    L << " (dirty)";
  }
  L << '\n';
  std::string origin = build_info->git_origin_url();
  if(origin.length()>=4 and origin.substr(origin.length()-4)==".git") {
    origin = origin.substr(0, origin.length()-4);
  }
  L << "URL : " << origin << "/commit/" << build_info->git_commit_sha1() << '\n';
  if(build_info->compiled_with_sse()) {
    L << "Built with :";
    if(build_info->compiled_with_avx512f())L << " AVX-512F";
    if(build_info->compiled_with_avx512dq())L << "/DQ";
    if(build_info->compiled_with_avx512ifma())L << "/IFMA";
    if(build_info->compiled_with_avx512pf())L << "/PF";
    if(build_info->compiled_with_avx512er())L << "/ER";
    if(build_info->compiled_with_avx512cd())L << "/CD";
    if(build_info->compiled_with_avx512bw())L << "/BW";
    if(build_info->compiled_with_avx512vl())L << "/VL";
    if(build_info->compiled_with_avx512vbmi())L << "/VBMI";
    if(build_info->compiled_with_avx512vbmi2())L << "/VBMI2";
    if(build_info->compiled_with_avx512vnni())L << "/VNNI";
    if(build_info->compiled_with_avx512bitalg())L << "/BITALG";
    if(build_info->compiled_with_avx512vpopcntdq())L << "/VPOPCNTDQ";
    if(build_info->compiled_with_avx512_4vnniw())L << "/4VNNIW";
    if(build_info->compiled_with_avx512_4fmaps())L << "/4FMAPS";
    if(build_info->compiled_with_avx2())L << " AVX2";
    if(build_info->compiled_with_avx())L << " AVX";
    if(build_info->compiled_with_fma3())L << " FMA3";
    if(build_info->compiled_with_fma4())L << " FMA4";
    if(build_info->compiled_with_sse4_2())L << " SSE4.2";
    if(build_info->compiled_with_sse4_1())L << " SSE4.1";
    if(build_info->compiled_with_ssse3())L << " SSSE3";
    if(build_info->compiled_with_sse3())L << " SSE3";
    if(build_info->compiled_with_sse2())L << " SSE2";
    if(build_info->compiled_with_sse())L << " SSE";
  }
  return L.str();
}

std::string calin::provenance::system_info::system_info_string(
  const calin::ix::provenance::system_info::HostAndProcessInfo* host_info)
{
  if(host_info == nullptr) {
    host_info = calin::provenance::system_info::the_host_info();
  }
  std::ostringstream L;
  L << host_info->uname_sysname() << " " << host_info->uname_release() << " on "
    << host_info->cpu_processor_brand() << '\n';
  if(host_info->cpu_has_sse()) {
    L << "CPU has :";
    if(host_info->cpu_has_avx512f())L << " AVX-512F";
    if(host_info->cpu_has_avx512dq())L << "/DQ";
    if(host_info->cpu_has_avx512ifma())L << "/IFMA";
    if(host_info->cpu_has_avx512pf())L << "/PF";
    if(host_info->cpu_has_avx512er())L << "/ER";
    if(host_info->cpu_has_avx512cd())L << "/CD";
    if(host_info->cpu_has_avx512bw())L << "/BW";
    if(host_info->cpu_has_avx512vl())L << "/VL";
    if(host_info->cpu_has_avx512vbmi())L << "/VBMI";
    if(host_info->cpu_has_avx512vbmi2())L << "/VBMI2";
    if(host_info->cpu_has_avx512vnni())L << "/VNNI";
    if(host_info->cpu_has_avx512bitalg())L << "/BITALG";
    if(host_info->cpu_has_avx512vpopcntdq())L << "/VPOPCNTDQ";
    if(host_info->cpu_has_avx512_4vnniw())L << "/4VNNIW";
    if(host_info->cpu_has_avx512_4fmaps())L << "/4FMAPS";
    if(host_info->cpu_has_avx2())L << " AVX2";
    if(host_info->cpu_has_avx())L << " AVX";
    if(host_info->cpu_has_fma3())L << " FMA3";
    if(host_info->cpu_has_fma4())L << " FMA4";
    if(host_info->cpu_has_sse4_2())L << " SSE4.2";
    if(host_info->cpu_has_sse4_1())L << " SSE4.1";
    if(host_info->cpu_has_ssse3())L << " SSSE3";
    if(host_info->cpu_has_sse3())L << " SSE3";
    if(host_info->cpu_has_sse2())L << " SSE2";
    if(host_info->cpu_has_sse())L << " SSE";
  }
  return L.str();
}

void calin::provenance::system_info::write_system_info_to_log(calin::util::log::Level level,
  calin::util::log::Logger* logger)
{
  auto L = LOG(level, logger);
  L << system_info_string();
}

bool calin::provenance::system_info::has_avx()
{
  unsigned a,b,c,d;
  __cpuid (0 /* vendor string */, a, b, c, d);
  unsigned max_frame = a;
  if(max_frame >= 1)
  {
    __cpuid (1 /* Processor Info and Feature Bits */, a, b, c, d);
    return c & bit_AVX;
  }
  return false;
}

bool calin::provenance::system_info::has_fma3()
{
  unsigned a,b,c,d;
  __cpuid (0 /* vendor string */, a, b, c, d);
  unsigned max_frame = a;
  if(max_frame >= 1)
  {
    __cpuid (1 /* Processor Info and Feature Bits */, a, b, c, d);
    return c & bit_FMA;
  }
  return false;
}

bool calin::provenance::system_info::has_avx2()
{
  unsigned a,b,c,d;
  __cpuid (0 /* vendor string */, a, b, c, d);
  unsigned max_frame = a;
  if(max_frame >= 7)
  {
    __cpuid_count (7, 0 /* Extended Features */, a, b, c, d);
    return b & bit_AVX2;
  }
  return false;
}

bool calin::provenance::system_info::has_avx512f()
{
  unsigned a,b,c,d;
  __cpuid (0 /* vendor string */, a, b, c, d);
  unsigned max_frame = a;
  if(max_frame >= 7)
  {
    __cpuid_count (7, 0 /* Extended Features */, a, b, c, d);
    return b & bit_AVX512F;
  }
  return false;
}

bool calin::provenance::system_info::vcl_uses_avx()
{
#if MAX_VECTOR_SIZE >= 256
#if INSTRSET >= 7
  return true;
#endif
#endif
  return false;
}

bool calin::provenance::system_info::vcl_uses_avx2()
{
#if MAX_VECTOR_SIZE >= 256
#if INSTRSET >= 8
  return true;
#endif
#endif
  return false;
}

bool calin::provenance::system_info::vcl_uses_avx512f()
{
#if MAX_VECTOR_SIZE >= 512
#if INSTRSET >= 9
  return true;
#endif
#endif
  return false;
}
