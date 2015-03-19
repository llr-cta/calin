#pragma once

#include <string>
#include <vector>

#include "math/accumulator.hpp"
#include "math/function.hpp"
#include "math/histogram.hpp"

namespace calin { namespace calib { namespace spe_fit {

class SingleElectronSpectrum: public math::Parameterizable
{
 public:
  virtual ~SingleElectronSpectrum();
};

class MultiElectronSpectrum: public math::Parameterizable
{
 public:
  using ConstVecRef = math::function::ConstVecRef;
  using VecRef = math::function::VecRef;
  using MatRef = math::function::MatRef;

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
  using ConstVecRef = math::function::ConstVecRef;
  using VecRef = math::function::VecRef;
  using MatRef = math::function::MatRef;

  PoissonGaussianMES(unsigned nmax = 10, bool force_calc_hessian = false);
  virtual ~PoissonGaussianMES();

  std::vector<math::ParameterAxis> parameters() override;
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

  using accumulator = math::LikelihoodAccumulator;

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
  using ConstVecRef = math::function::ConstVecRef;
  using VecRef = math::function::VecRef;
  using MatRef = math::function::MatRef;

  PoissonGaussianMES_HighAccuracy(double tol = 1e-200);
  virtual ~PoissonGaussianMES_HighAccuracy();

  std::vector<math::ParameterAxis> parameters() override;
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

class SPELikelihood: public math::MultiAxisFunction
{
 public:
  using ConstVecRef = math::function::ConstVecRef;
  using VecRef = math::function::VecRef;
  using MatRef = math::function::MatRef;

  SPELikelihood(MultiElectronSpectrum& mes_model,
                const math::SimpleHist& mes_data);
  SPELikelihood(MultiElectronSpectrum& mes_model,
                const math::SimpleHist& mes_data,
                const math::SimpleHist& ped_data);
  virtual ~SPELikelihood();
  std::vector<math::DomainAxis> domain_axes() override;
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
  const math::SimpleHist mes_data_;
  bool has_ped_data_ { false };
  const math::SimpleHist ped_data_;
};

} // namespace spe_fit

using spe_fit::SingleElectronSpectrum;
using spe_fit::MultiElectronSpectrum;
using spe_fit::PoissonGaussianMES;
using spe_fit::PoissonGaussianMES_HighAccuracy;
using spe_fit::SPELikelihood;

} } // namespace calin::calib
