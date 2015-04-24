/* 

   calin/calib/spe_fit.hpp -- Stephen Fegan -- 2015-03-01

   Functions to do fit to multi-electron spectrum in the "single PE"
   domain.

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

namespace calin { namespace calib { namespace spe_fit {

CALIN_TYPEALIAS(SingleElectronSpectrum, calin::math::pdf_1d::Parameterizable1DPDF);

class MultiElectronSpectrum: public calin::math::function::Parameterizable
{
 public:
  virtual ~MultiElectronSpectrum();
  virtual double pdf_ped(double x) = 0;
  virtual double pdf_gradient_ped(double x, VecRef gradient) = 0;
  virtual double pdf_gradient_hessian_ped(double x, VecRef gradient,
                                        MatRef hessian);
  virtual double pdf_mes(double x) = 0;
  virtual double pdf_gradient_mes(double x, VecRef gradient) = 0;
  virtual double pdf_gradient_hessian_mes(double x, VecRef gradient,
                                        MatRef hessian);

  virtual double intensity_pe() = 0;
  virtual double ped_rms_dc() = 0;
  virtual double ped_zero_dc() = 0;
  virtual double ses_mean_dc() = 0;
  virtual double ses_rms_pe() = 0;
  
  bool can_calculate_parameter_hessian() override;
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

  double pdf_mes(double x) override;
  double pdf_ped(double x) override;
  double pdf_gradient_mes(double x, VecRef gradient) override;
  double pdf_gradient_ped(double x, VecRef gradient) override;

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

class GeneralPoissonMES: public MultiElectronSpectrum
{
 public:
  GeneralPoissonMES(double x0, double dx, unsigned npoint,
                    SingleElectronSpectrum* ses,
                    calin::math::pdf_1d::Parameterizable1DPDF* ped,
                    unsigned nmax = 10,
                    bool adopt_ses = false, bool adopt_ped = false);

  virtual ~GeneralPoissonMES();

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

  double x0() const { return x0_; }
  double dx() const { return dx_; }
  double num_electrons_in_model() const { return nmax_; }
  std::vector<double> multi_electron_spectrum() const;
  std::vector<double> pedestal_spectrum() const;
  std::vector<double> n_electron_spectrum(unsigned n) const; // 1<=n<=nmax_
  std::vector<double> mes_n_electron_cpt(unsigned n) const; // 0<=n<=nmax_
  
 protected:  
  int ibin(double x) const;
  void set_cache();
  void hcvec_scale_and_multiply(double* ovec, const double* ivec1,
                                const double* ivec2, double scale = 1.0) const;
  void hcvec_scale_and_add(double* ovec, const double* ivec, double scale)
      const;
  
  SingleElectronSpectrum* ses_pdf_;
  calin::math::pdf_1d::Parameterizable1DPDF* ped_pdf_;
  bool adopt_ses_pdf_ = false;
  bool adopt_ped_pdf_ = false;
  unsigned nmax_ = 10;
  double intensity_pe_ = 1.0;

  double x0_ = 0;
  double dx_ = 0;
  unsigned nsample_ = 0;
  double* ped_spec_ = nullptr;
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
  std::vector<double*> mes_grad_;
  std::vector<fftw_plan> mes_grad_plan_rev_;
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
  MultiElectronSpectrum* mes_model_;
  unsigned npar_;
  const calin::math::histogram::SimpleHist mes_data_;
  bool has_ped_data_ { false };
  const calin::math::histogram::SimpleHist ped_data_;
};

} } } // namespace calin::calib::spe_fit
