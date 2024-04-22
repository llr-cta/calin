/*

   calin/math/fftw_util.cpp -- Stephen Fegan -- 2016-03-29

   Utility functions for fftw HC types

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

#include <cstdlib>

#include <fftw3.h>

#include <util/vcl.hpp>
#include <util/file.hpp>
#include <util/memory.hpp>
#include <math/fftw_util.hpp>

void calin::math::fftw_util::
hcvec_fftfreq(double* ovec, unsigned nsample, double d, bool imaginary_negative)
{
  double norm = 1.0/(d*nsample);
  double* r = ovec;
  double* c = ovec + nsample;
  unsigned k=0;
  *r++ = k++;
  c--;
  if(imaginary_negative) {
    while(r<c) {
      *c-- = -(*r++ = k++ * norm);
    }
  } else {
    while(r<c) {
      *c-- = *r++ = k++ * norm;
    }
  }
  if(r==c) {
    *r++ = k * norm;
  }
}

Eigen::VectorXd calin::math::fftw_util::
hcvec_fftfreq(unsigned nsample, double d, bool imaginary_negative)
{
  Eigen::VectorXd ovec(nsample);
  hcvec_fftfreq(ovec.data(), nsample, d, imaginary_negative);
  return ovec;
}

void calin::math::fftw_util::
hcvec_fftindex(int* ovec, unsigned nsample, bool imaginary_negative)
{
  int* r = ovec;
  int* c = ovec + nsample;
  int k=0;
  *r++ = k++;
  c--;
  if(imaginary_negative) {
    while(r<c) {
      *c-- = -(*r++ = k++);
    }
  } else {
    while(r<c) {
      *c-- = *r++ = k++;
    }
  }
  if(r==c) {
    *r++ = k;
  }
}

Eigen::VectorXi calin::math::fftw_util::
hcvec_fftindex(unsigned nsample, bool imaginary_negative)
{
  Eigen::VectorXi ovec(nsample);
  hcvec_fftindex(ovec.data(), nsample, imaginary_negative);
  return ovec;
}

#if INSTRSET >= 7
void calin::math::fftw_util::hcvec_scale_and_multiply(double* ovec, const double* ivec1,
  const double* ivec2, unsigned nsample, double scale)
{
  hcvec_scale_and_multiply_vcl<calin::util::vcl::VCLDoubleReal<calin::util::vcl::VCL256Architecture> >(ovec,
    ivec1, ivec2, nsample, scale);
  // hcvec_scale_and_multiply<double>(ovec, ivec1, ivec2, nsample, scale);
}

void calin::math::fftw_util::hcvec_multiply_and_add_real(double* ovec, const double* ivec1,
  const double* ivec2, double real_addand, unsigned nsample)
{
  hcvec_multiply_and_add_real_vcl<calin::util::vcl::VCLDoubleReal<calin::util::vcl::VCL256Architecture> >(ovec,
    ivec1, ivec2, real_addand, nsample);
  // hcvec_multiply_and_add_real<double>(ovec, ivec1, ivec2, real_addand, nsample);
}

void calin::math::fftw_util::hcvec_polynomial(double* ovec, const double* ivec,
  const std::vector<double>& p, unsigned nsample)
{
  hcvec_polynomial_vcl<calin::util::vcl::VCLDoubleReal<calin::util::vcl::VCL256Architecture> >(ovec, ivec, p, nsample);
  //hcvec_polynomial<double>(ovec, ivec, p, nsample);
}

void calin::math::fftw_util::hcvec_multi_stage_polynomial(double* ovec, double* ivec,
  const std::vector<const std::vector<double>*>& stage_p, unsigned nsample)
{
  hcvec_multi_stage_polynomial_vcl<calin::util::vcl::VCLDoubleReal<calin::util::vcl::VCL256Architecture> >(ovec, ivec, stage_p, nsample);
}

void calin::math::fftw_util::hcvec_gaussian_dft(double* ovec, double mean, double sigma, unsigned nsample)
{
  hcvec_gaussian_dft_vcl<calin::util::vcl::VCLDoubleReal<calin::util::vcl::VCL256Architecture> >(ovec, mean, sigma, nsample);
}

void calin::math::fftw_util::hcvec_2gaussian_dft(double* ovec, double mean, double sigma, double split, unsigned nsample)
{
  hcvec_2gaussian_dft_vcl<calin::util::vcl::VCLDoubleReal<calin::util::vcl::VCL256Architecture> >(ovec, mean, sigma, split, nsample);
  // hcvec_gaussian_dft_vcl<calin::util::vcl::VCLDoubleReal<calin::util::vcl::VCL256Architecture> >(ovec, mean, sigma, nsample);
}

void calin::math::fftw_util::hcvec_delta_dft(double* ovec, double x0, unsigned nsample)
{
  hcvec_delta_dft_vcl<calin::util::vcl::VCLDoubleReal<calin::util::vcl::VCL256Architecture> >(ovec, x0, nsample);
}

void calin::math::fftw_util::hcvec_delta_idft(double* ovec, double k0, double phase0, unsigned nsample)
{
  hcvec_delta_idft_vcl<calin::util::vcl::VCLDoubleReal<calin::util::vcl::VCL256Architecture> >(ovec, k0, phase0, nsample);
}

void calin::math::fftw_util::hcvec_delta_iq_idft(
  double* oivec, double* oqvec, double k0, double phase0, unsigned nsample)
{
  hcvec_delta_iq_idft_vcl<calin::util::vcl::VCLDoubleReal<calin::util::vcl::VCL256Architecture> >(oivec, oqvec, k0, phase0, nsample);
}
#endif

Eigen::VectorXd calin::math::fftw_util::fftw_r2hc(const Eigen::VectorXd& x,
  calin::ix::math::fftw_util::FFTWPlanningRigor fftw_rigor)
{
  int plan_flags = proto_planning_enum_to_fftw_flag(fftw_rigor);

  uptr_fftw_data x_copy { fftw_alloc_real(x.size()), fftw_free };
  assert(x_copy);

  uptr_fftw_data f { fftw_alloc_real(x.size()), fftw_free };
  assert(f);

  // Prepare the DFT
  uptr_fftw_plan r2hc_plan = {
    fftw_plan_r2r_1d(x.size(), x_copy.get(), f.get(),
                    FFTW_R2HC , plan_flags), fftw_destroy_plan };
  assert(r2hc_plan);

  std::copy(x.data(), x.data()+x.size(), x_copy.get());

  fftw_execute(r2hc_plan.get());

  return Eigen::Map<Eigen::VectorXd>(f.get(), x.size());
}

Eigen::VectorXd calin::math::fftw_util::fftw_hc2r(const Eigen::VectorXd& f,
  calin::ix::math::fftw_util::FFTWPlanningRigor fftw_rigor)
{
  int plan_flags = proto_planning_enum_to_fftw_flag(fftw_rigor);

  uptr_fftw_data f_copy { fftw_alloc_real(f.size()), fftw_free };
  assert(f_copy);

  uptr_fftw_data x { fftw_alloc_real(f.size()), fftw_free };
  assert(x);

  // Prepare the DFT
  uptr_fftw_plan hc2r_plan = {
    fftw_plan_r2r_1d(f.size(), f_copy.get(), x.get(),
                    FFTW_HC2R , plan_flags), fftw_destroy_plan };
  assert(hc2r_plan);

  std::copy(f.data(), f.data()+f.size(), f_copy.get());

  fftw_execute(hc2r_plan.get());

  return Eigen::Map<Eigen::VectorXd>(x.get(), f.size());
}

Eigen::VectorXd calin::math::fftw_util::hcvec_multiply_and_add_real(const Eigen::VectorXd& ivec1, const Eigen::VectorXd& ivec2, double real_addand)
{
  if(ivec1.size() != ivec2.size()) {
    throw std::runtime_error("Input vectors must have same number of elements");
  }
  Eigen::VectorXd ovec = ivec1;
  hcvec_multiply_and_add_real(ovec.data(), ivec1.data(), ivec2.data(), real_addand, ovec.size());
  return ovec;
}

Eigen::VectorXd calin::math::fftw_util::hcvec_scale_and_add_real(const Eigen::VectorXd& ivec, double scale, double real_addand)
{
  Eigen::VectorXd ovec = ivec;
  hcvec_scale_and_add_real(ovec.data(), scale, real_addand, ovec.size());
  return ovec;
}

double calin::math::fftw_util::hcvec_sum_real(const Eigen::VectorXd& ivec)
{
  return hcvec_sum_real(ivec.data(), ivec.size());
}

double calin::math::fftw_util::hcvec_avg_real(const Eigen::VectorXd& ivec)
{
  return hcvec_avg_real(ivec.data(), ivec.size());
}

Eigen::VectorXd calin::math::fftw_util::hcvec_gaussian_dft(double mean, double sigma, unsigned nsample, bool vcl)
{
  Eigen::VectorXd ovec(nsample);
  if(vcl) {
    hcvec_gaussian_dft(ovec.data(), mean, sigma, nsample);
  } else {
    hcvec_gaussian_dft<double>(ovec.data(), mean, sigma, nsample);
  }
  return ovec;
}

Eigen::VectorXd calin::math::fftw_util::hcvec_2gaussian_dft(double mean, double sigma, double split, unsigned nsample, bool vcl)
{
  Eigen::VectorXd ovec(nsample);
  if(vcl) {
    hcvec_2gaussian_dft(ovec.data(), mean, sigma, split, nsample);
  } else {
    hcvec_2gaussian_dft<double>(ovec.data(), mean, sigma, split, nsample);
  }
  return ovec;
}

Eigen::VectorXd calin::math::fftw_util::hcvec_delta_dft(double x0, unsigned nsample, bool vcl)
{
  Eigen::VectorXd ovec(nsample);
  if(vcl) {
    hcvec_delta_dft(ovec.data(), x0, nsample);
  } else {
    hcvec_delta_dft<double>(ovec.data(), x0, nsample);
  }
  return ovec;
}

Eigen::VectorXd calin::math::fftw_util::hcvec_delta_idft(double k0, double phase0, unsigned nsample, bool vcl)
{
  Eigen::VectorXd ovec(nsample);
  if(vcl) {
    hcvec_delta_idft(ovec.data(), k0, phase0, nsample);
  } else {
    hcvec_delta_idft<double>(ovec.data(), k0, phase0, nsample);
  }
  return ovec;
}

Eigen::VectorXd calin::math::fftw_util::hcvec_delta_idft_by_index(unsigned index, unsigned nsample, bool vcl)
{
  if(index >= nsample) {
    throw std::out_of_range("hcvec_delta_idft_by_index: index must be smaller than nsample");
  }
  double k0;
  double phase0;
  if(index < hcvec_num_real(nsample)) {
    k0 = index;
    phase0 = 0;
  } else {
    k0 = nsample - index;
    phase0 = M_PI_2;
  }
  Eigen::VectorXd ovec(nsample);
  if(vcl) {
    hcvec_delta_idft(ovec.data(), k0, phase0, nsample);
  } else {
    hcvec_delta_idft<double>(ovec.data(), k0, phase0, nsample);
  }
  return ovec;
}

void calin::math::fftw_util::hcvec_delta_iq_idft(Eigen::VectorXd& oivec, Eigen::VectorXd& oqvec,
  double k0, double phase0, unsigned nsample, bool vcl)
{
  oivec.resize(nsample);
  oqvec.resize(nsample);
  if(vcl) {
    hcvec_delta_iq_idft(oivec.data(), oqvec.data(), k0, phase0, nsample);
  } else {
    hcvec_delta_iq_idft<double>(oivec.data(), oqvec.data(), k0, phase0, nsample);
  }
}

void calin::math::fftw_util::hcvec_delta_iq_idft_by_index(Eigen::VectorXd& oivec, Eigen::VectorXd& oqvec,
  unsigned index, unsigned nsample, bool vcl)
{
  if(index >= hcvec_num_real(nsample)) {
    throw std::out_of_range("hcvec_delta_iq_idft_by_index: index must be smaller than number of frequencies");
  }
  oivec.resize(nsample);
  oqvec.resize(nsample);
  if(vcl) {
    hcvec_delta_iq_idft(oivec.data(), oqvec.data(), double(index), 0.0, nsample);
  } else {
    hcvec_delta_iq_idft<double>(oivec.data(), oqvec.data(), double(index), 0.0, nsample);
  }
}

Eigen::VectorXd calin::math::fftw_util::hcvec_psd_weight(unsigned nsample)
{
  Eigen::VectorXd ovec(nsample);
  hcvec_psd_weight(ovec.data(), nsample);
  return ovec;
}

Eigen::VectorXd calin::math::fftw_util::hcvec_to_psd(const Eigen::VectorXd& ivec, double dc_cpt)
{
  unsigned nsample = ivec.size();
  Eigen::VectorXd ovec(hcvec_num_real(nsample));
  hcvec_to_psd(ovec.data(), ivec.data(), nsample, dc_cpt);
  return ovec;
}

Eigen::VectorXd calin::math::fftw_util::hcvec_to_psd_no_square(const Eigen::VectorXd& ivec)
{
  unsigned nsample = ivec.size();
  Eigen::VectorXd ovec(hcvec_num_real(nsample));
  hcvec_to_psd_no_square(ovec.data(), ivec.data(), nsample);
  return ovec;
}

bool calin::math::fftw_util::load_wisdom_from_file(std::string filename)
{
  calin::util::file::expand_filename_in_place(filename);
  return fftw_import_wisdom_from_filename(filename.c_str()) != 0;
}

bool calin::math::fftw_util::
load_wisdom_from_proto(const calin::ix::math::fftw_util::FFTWWisdom& proto)
{
  return fftw_import_wisdom_from_string(proto.wisdom().c_str()) != 0;
}

bool calin::math::fftw_util::save_wisdom_to_file(std::string filename)
{
  calin::util::file::expand_filename_in_place(filename);
  return fftw_export_wisdom_to_filename(filename.c_str()) != 0;
}

bool calin::math::fftw_util::
save_wisdom_to_proto(calin::ix::math::fftw_util::FFTWWisdom& proto)
{
  char* data = fftw_export_wisdom_to_string();
  if(!data)return false;
  proto.set_wisdom(data);
  std::free(data);
  return true;
}

int calin::math::fftw_util::
proto_planning_enum_to_fftw_flag(calin::ix::math::fftw_util::FFTWPlanningRigor x)
{
  switch(x) {
    case calin::ix::math::fftw_util::ESTIMATE:    return FFTW_ESTIMATE;
    case calin::ix::math::fftw_util::MEASURE:     return FFTW_MEASURE;
    case calin::ix::math::fftw_util::PATIENT:     return FFTW_PATIENT;
    case calin::ix::math::fftw_util::EXHAUSTIVE:  return FFTW_EXHAUSTIVE;
    default: throw std::invalid_argument("Unknown planning rigor type.");
  }
}

Eigen::VectorXd calin::math::fftw_util::fftw_codelet_r2hc(const Eigen::VectorXd& x)
{
  using Real = calin::util::vcl::VCLMaxDoubleReal;
  auto* r = calin::util::memory::aligned_calloc<typename Real::real_vt>(x.size());
  auto* c = calin::util::memory::aligned_calloc<typename Real::real_vt>(x.size());
  for(unsigned i=0; i<x.size(); i++) {
    r[i] = x[i];
  }
  calin::math::fftw_util::FFTWCodelet<Real> codelet;
  codelet.r2hc(x.size(), r, c);
  Eigen::VectorXd f(x.size());
  for(unsigned i=0; i<x.size(); i++) {
    f[i] = c[i][0];
  }
  ::free(r);
  ::free(c);
  return f;
}

Eigen::VectorXd calin::math::fftw_util::fftw_codelet_hc2r(const Eigen::VectorXd& f)
{
  using Real = calin::util::vcl::VCLMaxDoubleReal;
  auto* r = calin::util::memory::aligned_calloc<typename Real::real_vt>(f.size());
  auto* c = calin::util::memory::aligned_calloc<typename Real::real_vt>(f.size());
  for(unsigned i=0; i<f.size(); i++) {
    c[i] = f[i];
  }
  calin::math::fftw_util::FFTWCodelet<Real> codelet;
  codelet.hc2r(f.size(), r, c);
  Eigen::VectorXd x(f.size());
  for(unsigned i=0; i<f.size(); i++) {
    x[i] = r[i][0];
  }
  ::free(r);
  ::free(c);
  return x;
}

Eigen::VectorXd calin::math::fftw_util::fftw_codelet_r2hc_float(const Eigen::VectorXd& x)
{
  using Real = calin::util::vcl::VCLMaxFloatReal;
  auto* r = calin::util::memory::aligned_calloc<typename Real::real_vt>(x.size());
  auto* c = calin::util::memory::aligned_calloc<typename Real::real_vt>(x.size());
  for(unsigned i=0; i<x.size(); i++) {
    r[i] = x[i];
  }
  calin::math::fftw_util::FFTWCodelet<Real> codelet;
  codelet.r2hc(x.size(), r, c);
  Eigen::VectorXd f(x.size());
  for(unsigned i=0; i<x.size(); i++) {
    f[i] = c[i][0];
  }
  ::free(r);
  ::free(c);
  return f;
}

Eigen::VectorXd calin::math::fftw_util::fftw_codelet_hc2r_float(const Eigen::VectorXd& f)
{
  using Real = calin::util::vcl::VCLMaxFloatReal;
  auto* r = calin::util::memory::aligned_calloc<typename Real::real_vt>(f.size());
  auto* c = calin::util::memory::aligned_calloc<typename Real::real_vt>(f.size());
  for(unsigned i=0; i<f.size(); i++) {
    c[i] = f[i];
  }
  calin::math::fftw_util::FFTWCodelet<Real> codelet;
  codelet.hc2r(f.size(), r, c);
  Eigen::VectorXd x(f.size());
  for(unsigned i=0; i<f.size(); i++) {
    x[i] = r[i][0];
  }
  ::free(r);
  ::free(c);
  return x;
}
