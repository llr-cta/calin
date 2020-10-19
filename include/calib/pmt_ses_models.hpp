/*

   calin/calib/pmt_ses_models.hpp -- Stephen Fegan -- 2017-04-24

   PMT single-electron spectrum models

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

#pragma once

#include <Eigen/Dense>

#include <math/function.hpp>
#include <math/pdf_1d.hpp>
#include <math/fftw_util.hpp>
#include <calib/pmt_ses_models.pb.h>

namespace calin { namespace calib { namespace pmt_ses_models {

class TwoGaussianSES:
  public calin::math::function::ReducedSpaceParameterizableSingleAxisFunction
{
public:
  TwoGaussianSES(double dx = 0);
  virtual ~TwoGaussianSES();

  double gain();
  double resolution();
  void moments(double& mx, double& mxx);
};

class ReparameterizedTwoGaussianSES:
  public calin::math::pdf_1d::Parameterizable1DPDF
{
public:
  ReparameterizedTwoGaussianSES(double dx = 0);
  virtual ~ReparameterizedTwoGaussianSES();

  unsigned num_parameters() override;
  std::vector<calin::math::function::ParameterAxis> parameters() override;
  Eigen::VectorXd parameter_values() override;
  void set_parameter_values(ConstVecRef values) override;

  calin::math::function::DomainAxis domain_axis() override;

  bool can_calculate_gradient() override;
  bool can_calculate_hessian() override;
  bool can_calculate_parameter_gradient() override;
  bool can_calculate_parameter_hessian() override;

  double value_1d(double x) override;
  double value_and_gradient_1d(double x,  double& dfdx) override;
  double value_gradient_and_hessian_1d(double x, double& dfdx,
                            double& d2fdx2) override;
  double value_and_parameter_gradient_1d(double x,  VecRef gradient) override;
  double value_parameter_gradient_and_hessian_1d(double x, VecRef gradient,
                                                 MatRef hessian) override;

  TwoGaussianSES* underlying_ses() { return &ses_; }
  Eigen::MatrixXd jabobian() const { return jacobian_lu_.inverse(); }
  Eigen::MatrixXd jabobian_inv() const { return jacobian_inv_; }

  Eigen::MatrixXd prob_lo_curv() const { return prob_lo_curv_; }
  Eigen::MatrixXd beta_lo_curv() const { return beta_lo_curv_; }
  Eigen::MatrixXd gain_curv() const { return gain_curv_; }
  Eigen::MatrixXd beta_curv() const { return beta_curv_; }
  Eigen::MatrixXd mxx_curv() const { return mxx_curv_; }

private:
  void update_cached_values();

  TwoGaussianSES ses_;
  double prob_lo_ = 0.176559;
  double beta_lo_ = 0.362319;
  double gain_    = 100;
  double beta_    = 0.451781;
  Eigen::MatrixXd jacobian_inv_ = Eigen::MatrixXd::Identity(4,4);
  Eigen::FullPivLU<Eigen::MatrixXd> jacobian_lu_;
  Eigen::MatrixXd prob_lo_curv_ = Eigen::MatrixXd::Zero(4,4);
  Eigen::MatrixXd beta_lo_curv_ = Eigen::MatrixXd::Zero(4,4);
  Eigen::MatrixXd gain_curv_ = Eigen::MatrixXd::Zero(4,4);
  Eigen::MatrixXd beta_curv_ = Eigen::MatrixXd::Zero(4,4);
  Eigen::MatrixXd mxx_curv_ = Eigen::MatrixXd::Zero(4,4);

  unsigned num_warnings_ = 5;
};

class LombardMartinPrescottPMTModel {
public:
  struct Tableau {
    Tableau(unsigned npoint_, calin::ix::math::fftw_util::FFTWPlanningRigor fftw_rigor = calin::ix::math::fftw_util::ESTIMATE):
      npoint(npoint_),
      pmf(fftw_alloc_real(npoint), fftw_free),
      dft(fftw_alloc_real(npoint), fftw_free),
      basis_dft(fftw_alloc_real(npoint), fftw_free),
      dft_to_pmf_bwd_plan(fftw_plan_r2r_1d(npoint, dft.get(), pmf.get(), FFTW_HC2R,
          calin::math::fftw_util::proto_planning_enum_to_fftw_flag(fftw_rigor)),
        fftw_destroy_plan)
    {
      calin::math::fftw_util::hcvec_delta_dft(basis_dft.get(), 1.0, npoint);
    }
    Eigen::VectorXd pmf_as_vec() const {
      Eigen::VectorXd outvec(npoint);
      std::copy(pmf.get(), pmf.get()+npoint, outvec.data());
      return outvec;
    }
    Eigen::VectorXd dft_as_vec() const {
      Eigen::VectorXd outvec(npoint);
      std::copy(dft.get(), dft.get()+npoint, outvec.data());
      return outvec;
    }

    unsigned npoint;
    calin::math::fftw_util::uptr_fftw_data pmf;
    calin::math::fftw_util::uptr_fftw_data dft;
    calin::math::fftw_util::uptr_fftw_data basis_dft;
    calin::math::fftw_util::uptr_fftw_plan dft_to_pmf_bwd_plan;
  };

  LombardMartinPrescottPMTModel(
    const calin::ix::calib::pmt_ses_models::LombardMartinPrescottPMTModelConfig& config,
    double precision = 1e-10,
    calin::ix::math::fftw_util::FFTWPlanningRigor fftw_rigor = calin::ix::math::fftw_util::ESTIMATE);

  std::vector<double> stage_0_pmf() const { return stage_0_pmf_; }
  std::vector<double> stage_n_pmf() const { return stage_n_pmf_; }

  std::vector<double> stage_0_pmf_zsa() const { return stage_0_pmf_zsa_; }
  std::vector<double> stage_n_pmf_downsampled(Tableau& tableau) const;

  double p0() const { return p0_; }
  double total_gain() const { return total_gain_; }
  double resolution() const { return resolution_; }

  static std::vector<double> polya_pmf(double mean, double rms_frac, double precision = 1e-10);
  static std::vector<double> multi_stage_polya_pmf(unsigned nstage, double mean, double rms_frac, unsigned rebinning = 0, double precision = 1e-10);
#ifndef SWIG
  static std::vector<double> multi_stage_pmf(Tableau& tableau, unsigned nstage,
    const std::vector<double>& pmf, unsigned rebinning = 0, double precision = 1e-10,
    bool suppress_wraparound_warning = false, unsigned* wraparound_warning_count = nullptr);
  static unsigned rebin_pmf(double* pmf, unsigned npmf, unsigned binning);
#endif
  static std::vector<double> multi_stage_pmf(unsigned npoint, unsigned nstage,
    const std::vector<double>& pmf, unsigned rebinning = 0, double precision = 1e-10,
    calin::ix::math::fftw_util::FFTWPlanningRigor fftw_rigor = calin::ix::math::fftw_util::ESTIMATE);
  static std::vector<double> half_gaussian_pmf(double mean, double precision = 1e-10);
  static void rebin_pmf(std::vector<double>& pmf, unsigned binning);

  double stage_0_gain() const { return stage_0_gain_; }
  double stage_n_gain() const { return stage_n_gain_; }

  double stage_0_Exx() const { return stage_0_Exx_; }
  double stage_n_Exx() const { return stage_n_Exx_; }

  void calc_ses(Tableau& tableau);
  Eigen::VectorXd calc_ses(unsigned npoint = 0);

#ifndef SWIG
  void calc_mes(Tableau& tableau, double mean, double rms_frac,
    const double* ped, bool ped_is_fft=false);
  void calc_mes(Tableau& tableau, const std::vector<double>& pe_pmf,
    const double* ped, bool ped_is_fft=false);

  Eigen::VectorXd calc_mes(double mean, double rms_frac, unsigned npoint,
    const double* ped, bool ped_is_fft=false);
  Eigen::VectorXd calc_mes(const std::vector<double>& pe_pmf, unsigned npoint,
    const double* ped, bool ped_is_fft=false);
#endif

  Eigen::VectorXd calc_mes(double mean, double rms_frac, unsigned npoint = 0);
  Eigen::VectorXd calc_mes(const std::vector<double>& pe_pmf, unsigned npoint = 0);

  Eigen::VectorXd calc_mes(double mean, double rms_frac, const Eigen::VectorXd& ped,
    bool ped_is_fft=false);
  Eigen::VectorXd calc_mes(const std::vector<double>& pe_pmf, const Eigen::VectorXd& ped,
    bool ped_is_fft=false);

  static calin::ix::calib::pmt_ses_models::LombardMartinPrescottPMTModelConfig nectarcam_config();
private:
  void calc_spectrum(Tableau& tableau,
    const std::vector<double>* pe_spec = nullptr,
    const double* ped = nullptr, bool ped_is_fft = false);
  void set_stage_n_gain(double stage_n_gain);

  calin::ix::calib::pmt_ses_models::LombardMartinPrescottPMTModelConfig config_;
  double precision_;
  calin::ix::math::fftw_util::FFTWPlanningRigor fftw_rigor_;

  double p0_ = 0;
  double total_gain_ = 0;
  double resolution_ = 0;

  double stage_0_gain_ = 0;
  double stage_n_gain_ = 0;

  double stage_0_Exx_ = 0;
  double stage_n_Exx_ = 0;

  std::vector<double> stage_0_pmf_;
  std::vector<double> stage_0_pmf_zsa_;
  std::vector<double> stage_n_pmf_;
};

// ********************************** OBSOLETE *********************************

class TwoGaussianSESConstrained:
  public calin::math::function::ReducedSpaceParameterizableSingleAxisFunction
{
public:
  TwoGaussianSESConstrained(double dx = 0);
  virtual ~TwoGaussianSESConstrained();
};

class TwoGaussianSESConstrained_Fast:
  public calin::math::function::ReducedSpaceParameterizableSingleAxisFunction
{
public:
  TwoGaussianSESConstrained_Fast(double dx = 0);
  virtual ~TwoGaussianSESConstrained_Fast();
};

} } } // namespace calin::calib::pmt_ses_models
