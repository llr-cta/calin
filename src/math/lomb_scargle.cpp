/*

   calin/math/lomb_scargle.cpp -- Stephen Fegan -- 2017-12-05

   Calculate Lomb-Scargle periodogram

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

#include <stdexcept>
#include <cmath>

#include <provenance/system_info.hpp>
#include <math/lomb_scargle.hpp>
#include <util/memory.hpp>

namespace {

inline void validate(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti)
{
  if(xi.size() != ti.size())
    throw std::runtime_error("Lomb-Scargle: xi and ti must be same size, "
      + std::to_string(xi.size()) + " != " + std::to_string(ti.size()));
}

} // anonymous namespace

std::pair<double, double> calin::math::lomb_scargle::
amplitudes(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti, double freq)
{
  validate(xi,ti);

  unsigned ni = xi.size();
  double xmean = 0;
  for(unsigned i=0; i<ni; i++){
    xmean += xi[i];
  }
  xmean /= double(ni);

  const double omega = 2*M_PI*freq;
  double CC = 0;
  double SS = 0;
  double CS = 0;
  double XC = 0;
  double XS = 0;
  for(unsigned i=0; i<ni; i++) {
    const double theta = omega*ti[i];
    const double x = xi[i] - xmean;
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    CC += c*c;
    SS += s*s;
    CS += c*s;
    XC += x*c;
    XS += x*s;
  }
  const double A = (XC*SS-XS*CS)/(SS*CC-CS*CS);
  const double B = (XS*CC-XC*CS)/(SS*CC-CS*CS);
  return {A,B};
}

double calin::math::lomb_scargle::
power(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti, double freq)
{
  validate(xi,ti);

  unsigned ni = xi.size();
  double xmean = 0;
  for(unsigned i=0; i<ni; i++){
    xmean += xi[i];
  }
  xmean /= double(ni);

  const double omega = 2*M_PI*freq;
  double CC = 0;
  double SS = 0;
  double CS = 0;
  double XC = 0;
  double XS = 0;
  for(unsigned i=0; i<ni; i++) {
    const double theta = omega*ti[i];
    const double x = xi[i] - xmean;
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    CC += c*c;
    SS += s*s;
    CS += c*s;
    XC += x*c;
    XS += x*s;
  }
  const double A = (XC*SS-XS*CS)/(SS*CC-CS*CS);
  const double B = (XS*CC-XC*CS)/(SS*CC-CS*CS);
  return XC*A+XS*B;
}

Eigen::VectorXd calin::math::lomb_scargle::
periodogram(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti,
  double freq_lo, double freq_hi, double delta_freq)
{
  return periodogram_fast(xi, ti, freq_lo, freq_hi, delta_freq);
}

Eigen::VectorXd calin::math::lomb_scargle::
periodogram_slow(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti,
  double freq_lo, double freq_hi, double delta_freq)
{
  validate(xi,ti);
  unsigned nfreq = unsigned((freq_hi-freq_lo)/delta_freq);
  Eigen::VectorXd periodogram(nfreq);
  for(unsigned ifreq=0; ifreq<nfreq; ifreq++) {
    double freq = freq_lo + delta_freq*ifreq;
    periodogram[ifreq] = power(xi, ti, freq);
  }
  return periodogram;
}

Eigen::VectorXd calin::math::lomb_scargle::
periodogram_fast(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti,
  double freq_lo, double freq_hi, double delta_freq, unsigned renormalize_nfreq)
{
  validate(xi,ti);
  const unsigned nfreq = unsigned((freq_hi-freq_lo)/delta_freq);
  const unsigned ni = xi.size();

  Eigen::VectorXd cdf(ni);
  Eigen::VectorXd sdf(ni);
  const double delta_omega = 2*M_PI*delta_freq;
  double xmean = 0;
  for(unsigned i=0; i<ni; i++) {
    xmean += xi[i];
    const double delta_theta = delta_omega*ti[i];
    cdf[i] = std::cos(delta_theta);
    sdf[i] = std::sin(delta_theta);
  }

  Eigen::VectorXd xi_centered(ni);
  xmean /= double(ni);
  for(unsigned i=0; i<ni; i++) {
    xi_centered[i] = xi[i]-xmean;
  }

  Eigen::VectorXd cf(ni);
  Eigen::VectorXd sf(ni);
  Eigen::VectorXd periodogram(nfreq);
  for(unsigned ifreq=0; ifreq<nfreq; ifreq++) {

    // Renormalize the cosine and sine vectors as requested
    if(ifreq==0 or (renormalize_nfreq!=0 and ifreq%renormalize_nfreq==0)) {
      const double omega = 2*M_PI*(freq_lo + delta_freq*ifreq);
      for(unsigned i=0; i<ni; i++) {
        const double theta = omega*ti[i];
        cf[i] = std::cos(theta);
        sf[i] = std::sin(theta);
      }
    }

    double CC = 0;
    double SS = 0;
    double CS = 0;
    double XC = 0;
    double XS = 0;
    for(unsigned i=0; i<ni; i++) {
      const double c = cf[i];
      const double s = sf[i];

      CC += c*c;
      SS += s*s;
      CS += c*s;
      const double x = xi_centered[i];
      XC += x*c;
      XS += x*s;

      cf[i] = c*cdf[i] - s*sdf[i];
      sf[i] = c*sdf[i] + s*cdf[i];
    }
    const double A = (XC*SS-XS*CS)/(SS*CC-CS*CS);
    const double B = (XS*CC-XC*CS)/(SS*CC-CS*CS);
    periodogram[ifreq] = XC*A+XS*B;
  }
  return periodogram;
}

Eigen::VectorXd calin::math::lomb_scargle::
frequencies(const Eigen::VectorXd& periodogram, double freq_lo, double delta_freq)
{
  unsigned nfreq = periodogram.size();
  Eigen::VectorXd freq(nfreq);
  for(unsigned ifreq=0; ifreq<nfreq; ifreq++)freq[ifreq] = freq_lo + delta_freq*ifreq;
  return freq;
}
