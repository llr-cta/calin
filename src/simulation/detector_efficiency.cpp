/*

   calin/simulation/detector_efficiency.cpp -- Stephen Fegan -- 2016-10-15

   Classes to do 1D linear or exponential interpolation over detector efficiency
   and Cherenkov bandwidth tables.

   Originally from EGS5 ACT code, see header below.

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


// DetectorEfficiency.hpp - Classes to handle efficiency of detector
// - inculding absoption of Cherenkov light along the path
// Stephen Fegan - sfegan@llr.in2p3.fr - November 2012
// $Id: DetectorEfficiency.cpp 5422 2013-06-26 14:01:03Z sfegan $

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include <calin_global_definitions.hpp>
#include <util/log.hpp>
#include <util/file.hpp>
#include <provenance/chronicle.hpp>
#include <math/special.hpp>
#include <util/string.hpp>
#include <simulation/detector_efficiency.hpp>

using calin::math::special::SQR;
using calin::util::file::expand_filename;
using calin::util::string::chomp;
using calin::util::string::from_string;
using namespace calin::util::log;
using namespace calin::simulation::detector_efficiency;
using namespace calin::math::interpolation_1d;

#define EV_NM 1239.84193009239 // gunits: c/(ev/h) -> nm

// ----------------------------------------------------------------------------
// AtmosphericAbsorption
// ----------------------------------------------------------------------------

AtmosphericAbsorption::AtmosphericAbsorption()
{
  // nothing to see here
}

AtmosphericAbsorption::
AtmosphericAbsorption(const std::string& filename, OldStyleAtmObsFlag flag,
    double ground_level_km, double spacing_km)
{
  std::string fn = expand_filename(filename);
  std::ifstream stream(fn.c_str());
  if(!stream.good())
    throw std::runtime_error("Could not open: "+fn);

  auto* file_record = calin::provenance::chronicle::register_file_open(fn,
    calin::ix::provenance::chronicle::AT_READ, __PRETTY_FUNCTION__);
  std::string line;
  std::getline(stream,line); // skip first line
  file_record->set_comment(line);
  double e = 0;
  double z = 0;
  InterpLinear1D abs;
  std::getline(stream,line);
  while(stream)
  {
    if(line[1] != ' ')
    {
      std::istringstream lstream(line);
      if(e)
      {
        e_ev_.push_back(e);
        absorption_.push_back(abs);
      }
      double lambda;
      lstream >> lambda;
      e = EV_NM/lambda;
      z = 0;
      abs.clear();
    }
    else
    {
      std::istringstream lstream(line);
      double a;
      lstream >> a;
      while(lstream)
      {
        abs.insert(z*100000.0, a);
        z += spacing_km;
        lstream >> a;
      }
    }
    std::getline(stream,line);
  }
  if(e)
  {
    e_ev_.push_back(e);
    absorption_.push_back(abs);
  }
  calin::provenance::chronicle::register_file_close(file_record);
}

AtmosphericAbsorption::AtmosphericAbsorption(const std::string& filename,
  std::vector<double> levels_cm)
{
  std::string fn = expand_filename(filename);
  std::ifstream stream(fn.c_str());
  if(!stream.good())
    throw std::runtime_error("Could not open: "+fn);

  auto* file_record = calin::provenance::chronicle::register_file_open(fn,
    calin::ix::provenance::chronicle::AT_READ, __PRETTY_FUNCTION__);
  std::string line;
  std::string comment;
  std::getline(stream,line);

  while(levels_cm.empty() and stream)
  {
    std::istringstream line_stream(line);
    std::string token;
    line_stream >> token;
    if(token.empty()) {
      comment += '\n';
      // Empty line : move on to next line
      goto next_header_line;
    } else if(token[0]=='#') {
      // Comment line : look for pattern of H2= and H1= and import levels
      comment += line;
      comment += '\n';
      line_stream >> token;
      if(token != "H2=")goto next_header_line;
      double h_km;
      if(not (line_stream >> h_km))goto next_header_line;
      levels_cm.emplace_back(h_km * 1.0e5);
      line_stream >> token;
      if(token == ",")line_stream >> token;
      if(token != "H1=")goto next_header_line;
      while(line_stream >> h_km)
        levels_cm.emplace_back(h_km * 1.0e5);
    } else {
      // Start of data : break and process data below
      break;
    }
next_header_line:
    std::getline(stream,line);
  }

  if(levels_cm.size() < 2) {
    calin::provenance::chronicle::register_file_close(file_record);
    throw std::runtime_error("Must have 2 or more levels in absorption file: "+
      fn);
  }

  while(stream)
  {
    std::istringstream line_stream(line);
    std::string token;
    line_stream >> token;
    if(not token.empty())
    {
      double lambda;
      if(!from_string(token, lambda)) {
        throw std::runtime_error("Could not extract wavelength from line : " +
          line);
      }
      double e;
      e = EV_NM/lambda;

      auto ilevel = levels_cm.begin();
      InterpLinear1D abs;
      abs.insert(*ilevel, 0.0);
      for(++ilevel; ilevel!=levels_cm.end(); ++ilevel)
      {
        double a;
        if(not(line_stream >> a))
          throw std::runtime_error("Insuffient entries for lambda="
            + std::to_string(lambda));
        abs.insert(*ilevel, a);
      }

      e_ev_.emplace_back(e);
      absorption_.emplace_back(abs);
    }
    std::getline(stream,line);
  }
  file_record->set_comment(comment);
  calin::provenance::chronicle::register_file_close(file_record);
}

void AtmosphericAbsorption::set_zref(double zref)
{
  for(unsigned ie=0;ie<e_ev_.size();ie++) {
    double tau0 = absorption_[ie](zref);
    absorption_[ie] -= tau0;
  }
}

InterpLinear1D AtmosphericAbsorption::optical_depth_for_altitude(double h) const
{
  if(h < absorption_.front().xi(0))
    throw std::out_of_range("Altitude " + std::to_string(h) +
      " below lowest level available (" +
      std::to_string(absorption_.front().xi(0)) + ")");
  InterpLinear1D abs;
  for(unsigned ie=0;ie<e_ev_.size();ie++)
    abs.insert(e_ev_[ie], absorption_[ie](h));
  return abs;
}

double AtmosphericAbsorption::optical_depth_for_altitude_and_energy(double h, double e) const
{
  if(h < absorption_.front().xi(0))
    throw std::out_of_range("Altitude " + std::to_string(h) +
      " below lowest level available (" +
      std::to_string(absorption_.front().xi(0)) + ")");
  if(e > e_ev_[0])
    throw std::out_of_range("Energy " + std::to_string(e) +
      " below highest energy available (" +
      std::to_string(e_ev_[0]) + ")");
  auto e_below = std::upper_bound(e_ev_.begin(), e_ev_.end(), e, std::greater<double>());
  if(e < *e_below)
    throw std::out_of_range("Energy " + std::to_string(e) +
      " below lowest energy available (" +
      std::to_string(*e_below) + ")");

  int ie_below = e_below - e_ev_.begin();
  double tau_u = absorption_[ie_below-1](h);
  double tau_l = absorption_[ie_below](h);

  // std::cout << e << ' ' << *e_below << ' ' << *(e_below-1) << ' ' << tau_l << ' ' << tau_u << '\n';

  return (tau_u*(e-*e_below) + tau_l*(*(e_below-1)-e))/(*(e_below-1)-*e_below);
}

ACTEffectiveBandwidth AtmosphericAbsorption::
integrateBandwidth(double h0, double w0, const DetectionEfficiency& eff) const
{
  InterpLinear1D abs0 = optical_depth_for_altitude(h0);
  ACTEffectiveBandwidth bandwidth(w0);

#if 1
  bool obslevel = false;
  for(unsigned ih=0;ih<absorption_.front().nXY();ih++)
  {
    double h = absorption_.front().xi(ih);
    if(h==h0 && !obslevel)
    {
      obslevel = true;
    }
    else if(h>h0 && !obslevel)
    {
      ih--;
      h = h0;
      obslevel = true;
    }
#else
  for(double h=absorption_.front().xmin();
      h<absorption_.front().xmax(); h+=10000)
  {
#endif
    InterpLinear1D Y0;
    InterpLinear1D Y1;
    InterpLinear1D Y2;
    for(unsigned ie=0;ie<e_ev_.size();ie++)
  	{
  	  double e = e_ev_[ie];
  	  double mfp = std::fabs(absorption_[ie](h)-abs0(e));
  	  double abs = std::exp(-mfp/w0);
  	  Y0.insert(e, abs);
  	  Y1.insert(e, mfp/(w0*w0)*abs);
  	  Y2.insert(e, mfp*(0.5*mfp/w0-1.0)/(w0*w0*w0)*abs);
  	}
    Y0 *= eff;
    Y1 *= eff;
    Y2 *= eff;

    bandwidth_t y(Y0.integrate(), Y1.integrate(), Y2.integrate());
#ifdef ACT_LIGHT_YIELD_TAYLOR_SERIES_IN_LOG
    y.dn_dw   = y.dn_dw/y.n;
    y.d2n_dw2 = y.d2n_dw2/y.n - 0.5*SQR(y.dn_dw);
    y.n       = std::log(y.n);
#endif
    bandwidth.insert(h, y);
  }
  return bandwidth;
}

ACTEffectiveBandwidth AtmosphericAbsorption::
integrateBandwidth(double h0, double w0, const DetectionEfficiency& eff,
  double emin, double emax) const
{
  InterpLinear1D abs0 = optical_depth_for_altitude(h0);
  ACTEffectiveBandwidth bandwidth(w0);

  bool obslevel = false;
  for(unsigned ih=0;ih<absorption_.front().nXY();ih++)
  {
    double h = absorption_.front().xi(ih);
    if(h==h0 && !obslevel)
    {
      obslevel = true;
    }
    else if(h>h0 && !obslevel)
    {
      ih--;
      h = h0;
      obslevel = true;
    }
    InterpLinear1D Y0;
    InterpLinear1D Y1;
    InterpLinear1D Y2;
    for(unsigned ie=0;ie<e_ev_.size();ie++)
  	{
  	  double e = e_ev_[ie];
  	  double mfp = std::fabs(absorption_[ie](h)-abs0(e));
  	  double abs = std::exp(-mfp/w0);
  	  Y0.insert(e, abs);
  	  Y1.insert(e, mfp/(w0*w0)*abs);
  	  Y2.insert(e, mfp*(0.5*mfp/w0-1.0)/(w0*w0*w0)*abs);
  	}
    Y0 *= eff;
    Y1 *= eff;
    Y2 *= eff;

    bandwidth_t y(Y0.integrate(emin,emax), Y1.integrate(emin,emax), Y2.integrate(emin,emax));
#ifdef ACT_LIGHT_YIELD_TAYLOR_SERIES_IN_LOG
    y.dn_dw   = y.dn_dw/y.n;
    y.d2n_dw2 = y.d2n_dw2/y.n - 0.5*SQR(y.dn_dw);
    y.n       = std::log(y.n);
#endif
    bandwidth.insert(h, y);
  }
  return bandwidth;
}

calin::math::spline_interpolation::TwoDimensionalCubicSpline*
AtmosphericAbsorption::integrate_bandwidth_to_spline(
  double h0, const DetectionEfficiency& eff,
  std::vector<double> h, std::vector<double> w, double emin, double emax) const
{
  if(emin == 0 and emax == 0) {
    emin = std::numeric_limits<double>::infinity();
    for(auto e : eff.all_xi()) {
      emin = std::min(emin, e);
      emax = std::max(emax, e);
    }
  }

  if(h.empty()) {
    for(double ih=h0; ih<=levels_cm().back(); ih+=1e5) {
      h.push_back(ih);
    }
  }

  Eigen::VectorXd hh = calin::std_to_eigenvec(h);
  Eigen::VectorXd ww = calin::std_to_eigenvec(w);
  Eigen::MatrixXd bw(ww.size(), hh.size());

  InterpLinear1D abs0 = optical_depth_for_altitude(h0);

  for(unsigned ih=0;ih<hh.size();ih++)
  {
    double h = hh(ih);
    std::vector<InterpLinear1D> Y0(w.size());
    for(unsigned ie=0;ie<e_ev_.size();ie++)
  	{
  	  double e = e_ev_[ie];
  	  double od = std::max(absorption_[ie](h)-abs0(e), 0.0);
      for(unsigned iw=0;iw<w.size();iw++) {
	      double abs = std::exp(-od/w[iw]);
	      Y0[iw].insert(e, abs);
      }
  	}

    for(unsigned iw=0;iw<w.size();iw++) {
      Y0[iw] *= eff;
      bw(iw,ih) = Y0[iw].integrate(emin,emax);
    }
  }

  return new calin::math::spline_interpolation::TwoDimensionalCubicSpline(hh,ww,bw);
}

std::vector<double> AtmosphericAbsorption::levels_cm() const
{
  return absorption_.front().all_xi();
}

// ----------------------------------------------------------------------------
// DetectionEfficiency
// ----------------------------------------------------------------------------

DetectionEfficiency::DetectionEfficiency(double const_eff):
  InterpLinear1D(const_eff)
{
  // nothing to see here
}

DetectionEfficiency::DetectionEfficiency(const std::string& filename):
  InterpLinear1D(1.0)
{
  scaleEffFromFile(filename);
}

void DetectionEfficiency::scaleEff(const InterpLinear1D& eff)
{
  *static_cast<InterpLinear1D*>(this) *= eff;
}

void DetectionEfficiency::scaleEffByConst(double c)
{
  *static_cast<InterpLinear1D*>(this) *= c;
}

void DetectionEfficiency::
scaleEffFromFile(const std::string& filename)
{
  InterpLinear1D eff_fn;
  eff_fn.insert_from_2column_file_with_filter(filename,
    [](double& lambda_in_e_out, double& eff) {
      lambda_in_e_out = EV_NM / lambda_in_e_out; return true; });
  scaleEff(eff_fn);
}

void DetectionEfficiency::
scaleEffFromOldStyleFile(const std::string& filename,
		 double lambda0_nm, double dlambda_nm)
{
  std::ifstream stream(filename.c_str());
  auto* file_record = calin::provenance::chronicle::register_file_open(filename,
    calin::ix::provenance::chronicle::AT_READ, __PRETTY_FUNCTION__);
  InterpLinear1D eff_fn;
  std::string line;
  std::getline(stream,line); // skip first line
  line += '\n';
  file_record->set_comment(line);
  double lambda = lambda0_nm;
  double eff;
  stream >> eff;
  while(stream)
  {
    double e = EV_NM / lambda;
    eff_fn.insert(e, eff);
    lambda += dlambda_nm;
    stream >> eff;
  }
  scaleEff(eff_fn);
  calin::provenance::chronicle::register_file_close(file_record);
}

// ----------------------------------------------------------------------------
// AngularEfficiency
// ----------------------------------------------------------------------------

AngularEfficiency::AngularEfficiency(double const_eff):
  InterpLinear1D(const_eff)
{
  // nothing to see here
}

AngularEfficiency::AngularEfficiency(const std::string& filename):
  InterpLinear1D(1.0)
{
  this->insert_from_2column_file_with_filter(filename,
    [](double& theta_in_w_out, double& eff) {
      theta_in_w_out = std::cos(theta_in_w_out/180.0*M_PI); return true; });
}

void AngularEfficiency::scaleEff(const InterpLinear1D& eff)
{
  *static_cast<InterpLinear1D*>(this) *= eff;
}

void AngularEfficiency::scaleEffByConst(double c)
{
  *static_cast<InterpLinear1D*>(this) *= c;
}

void AngularEfficiency::
scaleEffFromFile(const std::string& filename)
{
  AngularEfficiency eff(filename);
  scaleEff(eff);
}

// ----------------------------------------------------------------------------
// ACTEffectiveBandwidth
// ----------------------------------------------------------------------------

ACTEffectiveBandwidth::ACTEffectiveBandwidth(double w0)
  : Interpolation1D(), m_w0(w0)
{
  // nothing to see here
}

double ACTEffectiveBandwidth::bandwidth(double h, double w) const
{
  bandwidth_t _y = y(h);
  double dw = w-m_w0;
#ifdef ACT_LIGHT_YIELD_TAYLOR_SERIES_IN_LOG
  return std::exp(_y.n + dw*(_y.dn_dw + _y.d2n_dw2*dw));
#else
  return _y.n + dw*(_y.dn_dw + _y.d2n_dw2*dw);
#endif
}

// ----------------------------------------------------------------------------
// PEAmplitudeGenerator
// ----------------------------------------------------------------------------

PEAmplitudeGenerator::~PEAmplitudeGenerator()
{
  // nothing to see here
}

SplinePEAmplitudeGenerator::SplinePEAmplitudeGenerator(
    const Eigen::VectorXd& q, const Eigen::VectorXd& dp_dq,
    SplineMode spline_mode, calin::math::rng::RNG* rng, bool adopt_rng):
  PEAmplitudeGenerator(),
  spline_(make_spline(q, dp_dq, spline_mode)), spline_mode_(spline_mode),
  rng_(rng==nullptr ? new calin::math::rng::RNG(__PRETTY_FUNCTION__, "Amplitude generation") : rng),
  adopt_rng_(rng==nullptr ? true : adopt_rng)
{
  // nothing to see here
}

calin::math::spline_interpolation::CubicSpline*
SplinePEAmplitudeGenerator::make_spline(const Eigen::VectorXd& q,
    const Eigen::VectorXd& dp_dq, SplineMode spline_mode,
    bool regularize_spline, bool extend_linear_rhs)
{
  std::vector<double> x;
  std::vector<double> y = calin::eigen_to_stdvec(q);
  x.push_back(0);
  double P = 0;
  for(unsigned iq=0;iq<q.size()-1;++iq) {
    P += 0.5*(dp_dq[iq] + dp_dq[iq+1]) * (q[iq+1] - q[iq]);
    x.push_back(P);
  }
  if(P>1.0) {
    LOG(WARNING) << "Normalization of probability distribution exceeds unity";
  }
  switch(spline_mode) {
  case SM_LINEAR:
    std::transform(x.begin(), x.end(), x.begin(), [](double x) { return x; });
    break;
  case SM_LOG:
    std::transform(x.begin(), x.end(), x.begin(), [](double x) { return -std::log(1-x); });
    break;
  case SM_SQRT_LOG:
    std::transform(x.begin(), x.end(), x.begin(), [](double x) { return std::sqrt(-std::log(1-x)); });
    break;
  }
  auto* spline = new calin::math::spline_interpolation::CubicSpline(x, y,
    calin::math::spline_interpolation::BC_NATURAL, 0,
    calin::math::spline_interpolation::BC_NATURAL, 0);
  if(regularize_spline) {
    auto even_spline =
      spline->new_regularized_spline((spline->xmax()-spline->xmin())/(x.size() - 1));
    delete spline;
    spline = even_spline;
  }
  if(extend_linear_rhs) {
    spline->extend_linear_rhs();
  }
  return spline;
}

SplinePEAmplitudeGenerator::SplinePEAmplitudeGenerator(
    const calin::math::spline_interpolation::CubicSpline& spline,
    SplineMode spline_mode, calin::math::rng::RNG* rng, bool adopt_rng):
  PEAmplitudeGenerator(),
  spline_(new calin::math::spline_interpolation::CubicSpline(spline)),
  spline_mode_(spline_mode),
  rng_(rng==nullptr ? new calin::math::rng::RNG(__PRETTY_FUNCTION__, "Amplitude generation") : rng),
  adopt_rng_(rng==nullptr ? true : adopt_rng)
{
  // nothing to see here
}

SplinePEAmplitudeGenerator::~SplinePEAmplitudeGenerator()
{
  if(adopt_rng_) { delete rng_; }
}

double SplinePEAmplitudeGenerator::generate_amplitude()
{
  double x = rng_->uniform_double();
  switch(spline_mode_) {
  case SM_LINEAR:
    break;
  case SM_LOG:
    x = -std::log(x);
    break;
  case SM_SQRT_LOG:
    x = std::sqrt(-std::log(x));
    break;
  }
  return spline_->value(x);
}
