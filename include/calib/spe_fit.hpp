/*

   calin/calib/spe_fit.hpp -- Stephen Fegan -- 2015-03-01

   Functions to do fit to multi-electron spectrum in the "single PE"
   domain.

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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
#include <vector>
#include <cmath>

#include <fftw3.h>

#include "math/accumulator.hpp"
#include "math/function.hpp"
#include "math/pdf_1d.hpp"
#include "math/histogram.hpp"
#include "math/m_estimate.hpp"
#include "math/data_modeling.hpp"
#include "calib/spe_fit.pb.h"
#include "math/fftw_util.hpp"
#include <calib/pmt_ses_models.hpp>

namespace calin { namespace calib { namespace spe_fit {

class MultiElectronSpectrum: public calin::math::function::Parameterizable
{
 public:
  virtual ~MultiElectronSpectrum();
  virtual double pdf_ped(double x) = 0;
  virtual double pdf_gradient_ped(double x, VecRef gradient) = 0;
  virtual double pdf_gradient_hessian_ped(double x, VecRef gradient,
                                        MatRef hessian) = 0;
  virtual double pdf_mes(double x) = 0;
  virtual double pdf_gradient_mes(double x, VecRef gradient) = 0;
  virtual double pdf_gradient_hessian_mes(double x, VecRef gradient,
                                        MatRef hessian) = 0;

  virtual double intensity_pe() = 0;
  virtual double ped_rms_dc() = 0;
  virtual double ped_zero_dc() = 0;
  virtual double ses_mean_dc() = 0;
  virtual double ses_rms_pe() = 0;
};

class PoissonGaussianMES: public MultiElectronSpectrum
{
 public:
  PoissonGaussianMES(unsigned nmax = 10, bool force_calc_hessian = false);
  virtual ~PoissonGaussianMES();

  unsigned num_parameters() override;
  std::vector<calin::math::function::ParameterAxis> parameters() override;
  Eigen::VectorXd parameter_values() override;
  void set_parameter_values(ConstVecRef values) override;
  bool can_calculate_parameter_gradient() override;
  bool can_calculate_parameter_hessian() override;

  double pdf_ped(double x) override;
  double pdf_gradient_ped(double x, VecRef gradient) override;
  double pdf_gradient_hessian_ped(double x, VecRef gradient,
                                MatRef hessian) override;

  double pdf_mes(double x) override;
  double pdf_gradient_mes(double x, VecRef gradient) override;
  double pdf_gradient_hessian_mes(double x, VecRef gradient,
                                MatRef hessian) override;

  double intensity_pe() override { return intensity_pe_; };
  double ped_rms_dc() override { return ped_rms_dc_; }
  double ped_zero_dc() override { return ped_zero_dc_; }
  double ses_mean_dc() override { return gain_dc_pe_; }
  double ses_rms_pe() override { return ses_rms_pe_; }

 protected:
  PoissonGaussianMES(unsigned nmax, unsigned nmin,
                     bool force_calc_hessian = false);

  CALIN_TYPEALIAS(accumulator, calin::math::accumulator::LikelihoodAccumulator);

  void calc_cached_vars(bool calc_hessian = false);

  unsigned nmax_       = 10;
  unsigned nmin_       = 0;
  double intensity_pe_ = 1.0;
  double gain_dc_pe_   = 1.0;
  double ses_rms_pe_   = 0.4;
  double ped_rms_dc_   = 0.2;
  double ped_zero_dc_  = 0.0;

  double b2_;
  double s2_;
  double g2_;
  double log_s_;

  std::vector<double> C1_;
  std::vector<double> C2_;
  std::vector<double> C3_;

  std::vector<double> dC1_dF_;
  std::vector<double> dC1_ds_;
  std::vector<double> dC1_dg_;
  std::vector<double> dC1_db_;
  std::vector<double> dC2_ds_;
  std::vector<double> dC2_dg_;
  std::vector<double> dC2_db_;
  std::vector<double> dC3_dg_;

  bool force_calc_hessian_ = false;
  bool hessian_elements_good_ = false;

  std::vector<double> d2C1_dF2_;
  std::vector<double> d2C1_ds2_;
  std::vector<double> d2C1_dg2_;
  std::vector<double> d2C1_db2_;
  std::vector<double> d2C1_dsdg_;
  std::vector<double> d2C1_dsdb_;
  std::vector<double> d2C1_dgdb_;

  std::vector<double> d2C2_ds2_;
  std::vector<double> d2C2_dg2_;
  std::vector<double> d2C2_db2_;
  std::vector<double> d2C2_dsdg_;
  std::vector<double> d2C2_dsdb_;
  std::vector<double> d2C2_dgdb_;

  std::vector<double> d2C3_dg2_;
  std::vector<double> d2C3_dsdg_;
  std::vector<double> d2C3_dgdb_;
};

class PoissonGaussianMES_HighAccuracy: public MultiElectronSpectrum
{
 public:
  PoissonGaussianMES_HighAccuracy(double tol = 1e-200);
  virtual ~PoissonGaussianMES_HighAccuracy();

  unsigned num_parameters() override;
  std::vector<calin::math::function::ParameterAxis> parameters() override;
  Eigen::VectorXd parameter_values() override;
  void set_parameter_values(ConstVecRef values) override;
  bool can_calculate_parameter_gradient() override;
  bool can_calculate_parameter_hessian() override;

  double pdf_mes(double x) override;
  double pdf_ped(double x) override;
  double pdf_gradient_mes(double x, VecRef gradient) override;
  double pdf_gradient_ped(double x, VecRef gradient) override;
  double pdf_gradient_hessian_mes(double x, VecRef gradient,
                                MatRef hessian) override;
  double pdf_gradient_hessian_ped(double x, VecRef gradient,
                                MatRef hessian) override;

  double intensity_pe() override { return intensity_pe_; };
  double ped_rms_dc() override { return ped_rms_dc_; }
  double ped_zero_dc() override { return ped_zero_dc_; }
  double ses_mean_dc() override { return gain_dc_pe_; }
  double ses_rms_pe() override { return ses_rms_pe_; }

 private:
  double tol_;
  double intensity_pe_ = 1.0;
  double gain_dc_pe_   = 1.0;
  double ses_rms_pe_   = 0.4;
  double ped_rms_dc_   = 0.2;
  double ped_zero_dc_  = 0.0;
};

class FastSingleValueGeneralPoissonMES;

class GeneralPoissonMES: public MultiElectronSpectrum
{
public:
   CALIN_TYPEALIAS(config_type, calin::ix::calib::
     spe_fit::GeneralPoissonMESConfig);

  GeneralPoissonMES(double x0, double dx, unsigned npoint,
                    calin::math::function::ParameterizableSingleAxisFunction* ses,
                    calin::math::function::ParameterizableSingleAxisFunction* ped,
                    calin::ix::calib::spe_fit::GeneralPoissonMESConfig config =
                      default_config(),
                    bool adopt_ses = false, bool adopt_ped = false);

  virtual ~GeneralPoissonMES();

  unsigned num_parameters() override;
  unsigned num_intrinsic_parameters();
  std::vector<calin::math::function::ParameterAxis> parameters() override;
  Eigen::VectorXd parameter_values() override;
  void set_parameter_values(ConstVecRef values) override;
  bool can_calculate_parameter_gradient() override;
  bool can_calculate_parameter_hessian() override;

  double pdf_ped(double x) override;
  double pdf_gradient_ped(double x, VecRef gradient) override;
  double pdf_gradient_hessian_ped(double x, VecRef gradient,
                                MatRef hessian) override;

  double pdf_mes(double x) override;
  double pdf_gradient_mes(double x, VecRef gradient) override;
  double pdf_gradient_hessian_mes(double x, VecRef gradient,
                                MatRef hessian) override;

  double intensity_pe() override { return intensity_pe_; };
  double ped_rms_dc() override;
  double ped_zero_dc() override;
  double ses_mean_dc() override;
  double ses_rms_pe() override;

  double x0() const { return x0_; }
  double dx() const { return dx_; }
  unsigned num_electrons_in_model() const { return config_.num_pe_convolutions(); }

  double ses_x(unsigned isample) const { return (0.5+double(isample))*dx_; }
  double ped_x(unsigned isample) const {
    return x0_ + (0.5+double(isample))*dx_; }
  double mes_x(unsigned isample) const { return ped_x(isample); }

  unsigned nsample() const { return nsample_; }

  Eigen::VectorXd all_ses_x() const { Eigen::VectorXd x(nsample_);
    for(unsigned i=0;i<nsample_;i++) x[i] = ses_x(i);
    return x; }
  Eigen::VectorXd all_ped_x() const { Eigen::VectorXd x(nsample_);
    for(unsigned i=0;i<nsample_;i++) x[i] = ped_x(i);
    return x; }
  Eigen::VectorXd all_mes_x() const { return all_ped_x(); }

  Eigen::VectorXd multi_electron_spectrum() const;
  Eigen::VectorXd pedestal_spectrum() const;
  Eigen::VectorXd off_pedestal_spectrum() const;
  Eigen::VectorXd n_electron_spectrum(unsigned n) const; // 1<=n<=config_.num_pe_convolutions()
  Eigen::VectorXd n_electron_spectrum_with_pedestal(unsigned n) const; // 0<=n<=config_.num_pe_convolutions()
  Eigen::VectorXd single_electron_spectrum() const {
    return  n_electron_spectrum(1); }
  Eigen::VectorXd mes_n_electron_cpt(unsigned n) const; // 0<=n<=config_.num_pe_convolutions()

  Eigen::VectorXd multi_electron_spectrum_gradient(unsigned iparam) const;
  Eigen::VectorXd pedestal_spectrum_gradient(unsigned iparam) const;
  Eigen::VectorXd off_pedestal_spectrum_gradient(unsigned iparam) const;
  Eigen::VectorXd single_electron_spectrum_gradient(unsigned iparam) const;

  Eigen::VectorXd extract_ped_gradient_values(ConstVecRef gradient);
  Eigen::VectorXd extract_ses_gradient_values(ConstVecRef gradient);
  Eigen::MatrixXd extract_ped_hessian_values(ConstMatRef hessian);
  Eigen::MatrixXd extract_ses_hessian_values(ConstMatRef hessian);

  unsigned iparam_light_intensity() const { return 0; }
  unsigned iparam_off_ped_shift() const {
    assert(config_.include_on_off_ped_shift());
    return 1;
  }
  unsigned iparam_ped() const {
    if(config_.include_on_off_ped_shift())return 2;
    else return 1;
  }
  unsigned iparam_ses() const {
    return iparam_ped() + ped_pdf_->num_parameters();
  }

  calin::ix::calib::spe_fit::GeneralPoissonMESConfig config() { return config_; }
  static calin::ix::calib::spe_fit::GeneralPoissonMESConfig default_config();

protected:
  friend class FastSingleValueGeneralPoissonMES;

  double ses_x_oversampled(unsigned isample) const {
    return (0.5+double(isample))*dx_oversample_; }
  double ped_x_oversampled(unsigned isample) const {
    return x0_ + (0.5+double(isample))*dx_oversample_; }

  int ibin(double x) const;
  void set_cache(bool force = false);

  double downsample_one(double* data, unsigned ibin, double rescale = 1) const
  {
    double val = data[ibin*oversample_factor_];
    for(unsigned ioversample=1; ioversample<oversample_factor_; ioversample++)
      val += data[ibin*oversample_factor_ + ioversample];
    return val*inv_oversample_factor_*rescale;
  }

  void downsample_all(double* data_to, double* data_from, double rescale = 1) const
  {
    if(oversample_factor_ <= 1) {
      std::copy(data_from, data_from+nsample_, data_to);
    } else {
      for(unsigned isample=0; isample<nsample_; isample++) {
        data_to[isample] = downsample_one(data_from, isample, rescale);
      }
    }
  }

  calin::math::function::ParameterizableSingleAxisFunction* ses_pdf_;
  calin::math::function::ParameterizableSingleAxisFunction* ped_pdf_;
  bool adopt_ses_pdf_ = false;
  bool adopt_ped_pdf_ = false;
  calin::ix::calib::spe_fit::GeneralPoissonMESConfig config_ = default_config();
  double intensity_pe_ = 1.0;
  double off_ped_shift_dc_ = 0.0;

  double x0_ = 0;
  double dx_ = 0;
  unsigned nsample_ = 0;
  unsigned oversample_factor_ = 1;
  double inv_oversample_factor_ = 1;
  unsigned noversample_ = 0;
  double dx_oversample_ = 0;
  double* ped_spec_ = nullptr;
  double* off_spec_ = nullptr;
  double* off_dfdx_ = nullptr;
  double* ped_fft_ = nullptr;
  std::vector<double*> nes_fft_;
  double* mes_spec_ = nullptr;
  fftw_plan ses_plan_fwd_; // Forward FFT of SES in place in nes_fft_[0]
  fftw_plan ped_plan_fwd_; // Forward FFT of PED from ped_spec to ped_fft_
  fftw_plan mes_plan_rev_; // Reverse FFT of MES in place in mes_spec_
  std::vector<double*> ped_grad_fft_;
  std::vector<double*> ses_grad_fft_;
  std::vector<fftw_plan> ped_grad_plan_fwd_;
  std::vector<fftw_plan> ses_grad_plan_fwd_;
  std::vector<double*> ped_grad_;
  std::vector<double*> off_grad_;
  std::vector<double*> mes_grad_;
  std::vector<fftw_plan> mes_grad_plan_rev_;
  unsigned n_ses_norm_warning_ = 0;
};

class FastSingleValueGeneralPoissonMES: public MultiElectronSpectrum
{
public:
  FastSingleValueGeneralPoissonMES(GeneralPoissonMES* mes, bool adopt_mes = false);
  virtual ~FastSingleValueGeneralPoissonMES();

  unsigned num_parameters() override;
  std::vector<calin::math::function::ParameterAxis> parameters() override;
  Eigen::VectorXd parameter_values() override;
  void set_parameter_values(ConstVecRef values) override;
  bool can_calculate_parameter_gradient() override;
  bool can_calculate_parameter_hessian() override;

  double pdf_ped(double x) override;
  double pdf_gradient_ped(double x, VecRef gradient) override;
  double pdf_gradient_hessian_ped(double x, VecRef gradient,
                                MatRef hessian) override;

  double pdf_mes(double x) override;
  double pdf_gradient_mes(double x, VecRef gradient) override;
  double pdf_gradient_hessian_mes(double x, VecRef gradient,
                                MatRef hessian) override;

  double intensity_pe() override { return intensity_pe_; };
  double ped_rms_dc() override;
  double ped_zero_dc() override;
  double ses_mean_dc() override;
  double ses_rms_pe() override;

  double x0() const { return mes_->x0(); }
  double dx() const { return mes_->dx(); }
  unsigned num_electrons_in_model() const { return nes_pmf_.cols(); }

  double ses_x(unsigned isample) const { return mes_->ses_x(isample); }
  double ped_x(unsigned isample) const { return mes_->ped_x(isample); }
  double mes_x(unsigned isample) const { return mes_->mes_x(isample); }

  unsigned nsample() const { return nes_pmf_.rows(); }

  Eigen::VectorXd all_ses_x() const { return mes_->all_ses_x(); }
  Eigen::VectorXd all_ped_x() const { return mes_->all_ped_x(); }
  Eigen::VectorXd all_mes_x() const { return mes_->all_mes_x(); }

  Eigen::VectorXd multi_electron_spectrum() const { return nes_pmf_*nes_weight_; }
  Eigen::VectorXd multi_electron_spectrum_gradient() const { return nes_pmf_*nes_weight_deriv_; }
  Eigen::VectorXd pedestal_spectrum() const { return nes_pmf_.col(0); };
  Eigen::VectorXd off_pedestal_spectrum() const { return ped_pmf_; }
  Eigen::VectorXd n_electron_spectrum_with_pedestal(unsigned n) const {
    return nes_pmf_.col(n); }

  void update_from_general_mes();
  GeneralPoissonMES* general_mes() const { return mes_; }

  const Eigen::MatrixXd& nes_pmf_matrix() const { return nes_pmf_; }
  const Eigen::VectorXd& ped_pmf_vector() const { return ped_pmf_; }
  const Eigen::VectorXd& nes_poisson_weight_vector() const { return nes_weight_; }
protected:
  GeneralPoissonMES* mes_ = nullptr;
  bool adopt_mes_ = false;
  double intensity_pe_ = 1.0;
  Eigen::MatrixXd nes_pmf_;
  Eigen::VectorXd ped_pmf_;
  Eigen::VectorXd nes_weight_;
  Eigen::VectorXd nes_weight_deriv_;
};

class LombardMartinPrescottMES: public MultiElectronSpectrum
{
public:
  LombardMartinPrescottMES(double x0, unsigned npoint,
    calin::math::function::ParameterizableSingleAxisFunction* ped_pdf,
    const calin::ix::calib::spe_fit::LombardMartinPrescottMESConfig& config = default_config(),
    bool adopt_ped_pdf = false);
  LombardMartinPrescottMES(double x0, unsigned npoint,
      const calin::ix::calib::spe_fit::LombardMartinPrescottMESConfig& config,
      calin::math::function::ParameterizableSingleAxisFunction* ped_pdf = nullptr,
      bool adopt_ped_ped = false):
    LombardMartinPrescottMES(x0, npoint, ped_pdf, config, adopt_ped_ped)
  { /* nothing to see here */ }

  virtual ~LombardMartinPrescottMES();

  unsigned num_parameters() override;
  std::vector<calin::math::function::ParameterAxis> parameters() override;
  Eigen::VectorXd parameter_values() override;
  void set_parameter_values(ConstVecRef values) override;
  bool can_calculate_parameter_gradient() override;
  bool can_calculate_parameter_hessian() override;

  double pdf_ped(double x) override;
  double pdf_gradient_ped(double x, VecRef gradient) override;
  double pdf_gradient_hessian_ped(double x, VecRef gradient, MatRef hessian) override;

  double pdf_mes(double x) override;
  double pdf_gradient_mes(double x, VecRef gradient) override;
  double pdf_gradient_hessian_mes(double x, VecRef gradient, MatRef hessian) override;

  double intensity_pe() override;
  double ped_rms_dc() override;
  double ped_zero_dc() override;
  double ses_mean_dc() override;
  double ses_rms_pe() override;

  calin::ix::calib::spe_fit::LombardMartinPrescottMESConfig config() const { return config_; }
  calin::ix::calib::pmt_ses_models::LombardMartinPrescottPMTModelConfig pmt() const { return config_.pmt(); }

  double x0() const { return x0_; }
  unsigned npoint() const { return npoint_; }
  unsigned mes_npoint() const { return mes_npoint_; }

  Eigen::VectorXd pe_pmf() const;

  Eigen::VectorXd mes_pmf() const { return mes_pmf_; }
  Eigen::VectorXd off_pmf() const { return off_pmf_; }
  Eigen::VectorXd ses_pmf() const;
  Eigen::VectorXd ses_pmf_full_resolution() const;

  static calin::ix::calib::spe_fit::LombardMartinPrescottMESConfig default_config();
private:
  void calculate_mes();
  double mes_x(int i) const { return x0_ + double(i)*config_.sensitivity() - 0.5*config_.dx(); }
  double off_x(int i) const { return mes_x(i) + config_.on_off_ped_shift(); }
  int ibin(double x) { return std::round((x-x0_)*dx_inv_); }
  void rebin_spectrum(Eigen::VectorXd& pmf_out, const double* mes_in, unsigned nmes) const;

  calin::ix::calib::spe_fit::LombardMartinPrescottMESConfig config_;
  double x0_;
  double dx_inv_ = 0.0;

  unsigned npoint_;
  unsigned mes_npoint_;

  mutable calin::calib::pmt_ses_models::LombardMartinPrescottPMTModel::Tableau tableau_;
  calin::math::fftw_util::uptr_fftw_data ped_ { nullptr, fftw_free };

  double intensity_pe_ = 1.0;
  calin::math::function::ParameterizableSingleAxisFunction* ped_pdf_ = nullptr;
  bool adopt_ped_pdf_ = false;
  Eigen::VectorXd mes_pmf_;
  Eigen::VectorXd off_pmf_;

  double ped_sum_p_ = 0;
  double ped_sum_px_ = 0;
  double ped_sum_pxx_ = 0;

  double pmt_total_gain_ = 0;
  double pmt_resolution_ = 0;
};

class SPELikelihood: public calin::math::function::MultiAxisFunction
{
 public:
  SPELikelihood(MultiElectronSpectrum& mes_model,
                const calin::math::histogram::SimpleHist& mes_data);
  SPELikelihood(MultiElectronSpectrum& mes_model,
                const calin::math::histogram::SimpleHist& mes_data,
                const calin::math::histogram::SimpleHist& ped_data);
  virtual ~SPELikelihood();

  unsigned num_domain_axes() override;
  std::vector<calin::math::function::DomainAxis> domain_axes() override;
  double value(ConstVecRef x) override;
  bool can_calculate_gradient() override;
  double value_and_gradient(ConstVecRef x, VecRef gradient) override;
  bool can_calculate_hessian() override;
  double value_gradient_and_hessian(ConstVecRef x, VecRef gradient,
                                    MatRef hessian) override;
  double error_up() override { return 0.5; }

private:
  void calculate_mes();

  MultiElectronSpectrum* mes_model_;
  unsigned npar_;
  const calin::math::histogram::SimpleHist mes_data_;
  bool has_ped_data_ { false };
  const calin::math::histogram::SimpleHist ped_data_;
};

class SPERobust: public calin::math::function::MultiAxisFunction
{
 public:
  SPERobust(MultiElectronSpectrum& mes_model,
    const calin::math::histogram::SimpleHist& mes_data,
    calin::math::m_estimate::LikelihoodRhoFunction* rho = nullptr, bool adopt_rho = false);
  SPERobust(MultiElectronSpectrum& mes_model,
    const calin::math::histogram::SimpleHist& mes_data,
    const calin::math::histogram::SimpleHist& ped_data,
    calin::math::m_estimate::LikelihoodRhoFunction* rho = nullptr, bool adopt_rho = false);
  virtual ~SPERobust();

  unsigned num_domain_axes() override;
  std::vector<calin::math::function::DomainAxis> domain_axes() override;
  double value(ConstVecRef x) override;
  bool can_calculate_gradient() override;
  double value_and_gradient(ConstVecRef x, VecRef gradient) override;
  bool can_calculate_hessian() override;
  double value_gradient_and_hessian(ConstVecRef x, VecRef gradient,
                                    MatRef hessian) override;
  double error_up() override { return 0.5; }

 private:
  MultiElectronSpectrum* mes_model_;
  unsigned npar_;
  calin::math::m_estimate::LikelihoodRhoFunction* rho_ = nullptr;
  bool adopt_rho_ = false;
  calin::math::data_modeling::IID1DDataMEstimateLikelihoodFunction* mes_cost_ = nullptr;
  calin::math::data_modeling::IID1DDataMEstimateLikelihoodFunction* ped_cost_ = nullptr;
};



} } } // namespace calin::calib::spe_fit
