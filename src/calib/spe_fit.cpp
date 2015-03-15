#include <cmath>
#include <limits>
#include <cassert>

#include "calib/spe_fit.hpp"
#include "calib/pmt_model_pg.hpp"

using namespace calin::math;
using namespace calin::calib::spe_fit;
using namespace calin::calib::pmt_model_pg;

using calin::math::function::assign_parameters;

namespace {

static double SQR(double x)
{
  return x*x;
}

constexpr double c_gauss_norm = 0.5*M_2_SQRTPI*M_SQRT1_2;

} // anonymous namespace

SingleElectronSpectrum::~SingleElectronSpectrum()
{
  // nothing to see here
}

MultiElectronSpectrum::~MultiElectronSpectrum()
{
  // nothing to see here
}


double MultiElectronSpectrum::pdf_derivs_hessian_ped(double x, double* derivs,
                                                     double* hessian)
{
  assert(0);
}

double MultiElectronSpectrum::pdf_derivs_hessian_mes(double x, double* derivs,
                                                     double* hessian)
{
  assert(0);
}

bool MultiElectronSpectrum::can_calculate_parameter_hessian()
{
  return false;
}

// ============================================================================
//
// Fast PG model for single PE regime where maximum number of PEs is specified
//
// ============================================================================

PoissonGaussianMES::PoissonGaussianMES(unsigned nmax, bool force_calc_hessian):
    nmax_(nmax), C1_(nmax), C2_(nmax), C3_(nmax),
    dC1_dF_(nmax), dC1_ds_(nmax), dC1_dg_(nmax), dC1_db_(nmax), 
    dC2_ds_(nmax), dC2_dg_(nmax), dC2_db_(nmax), dC3_dg_(nmax),
    force_calc_hessian_(force_calc_hessian)
{
  calc_cached_vars();
}

PoissonGaussianMES::~PoissonGaussianMES()
{
  // nothing to see here
}

auto PoissonGaussianMES::parameters() -> std::vector<math::ParameterAxis>
{
  constexpr double tiny_val = std::numeric_limits<double>::min();
  constexpr double phuge_val = std::numeric_limits<double>::max();
  constexpr double nhuge_val = std::numeric_limits<double>::lowest();
  return { { "light_intensity", "PE", true, 0, false, 0 },
    { "ped_zero", "DC", false, 0, false, 0 },
    { "ped_width", "DC", true, tiny_val, false, 0 },
    { "gain", "DC/PE", false, 0, false, 0 }, // can be negative
    { "ses_width", "PE", true, tiny_val, false, 0 }
  };
}

std::vector<double> PoissonGaussianMES::parameter_values()
{
  return { intensity_pe_, ped_zero_dc_, ped_rms_dc_,
        gain_dc_pe_, ses_rms_pe_, };
}

void PoissonGaussianMES::set_parameter_values(const double* values)
{
  assign_parameters(values, intensity_pe_, ped_zero_dc_, ped_rms_dc_,
                    gain_dc_pe_, ses_rms_pe_);
  calc_cached_vars();
}

bool PoissonGaussianMES::can_calculate_parameter_derivs()
{
  return true;
}

bool PoissonGaussianMES::can_calculate_parameter_hessian()
{
  return true;
}

double PoissonGaussianMES::pdf_mes(double x)
{
  const double x0 = ped_zero_dc_;
  const double xc = x-x0;
  accumulator apdf;
  for(unsigned n=nmin_; n<nmax_; n++)
  {
    const double sqr_arg = xc-C3_[n];
    const double log_pn = C1_[n] - C2_[n]*SQR(sqr_arg);
    apdf.accumulate(std::exp(log_pn));
  }
  return c_gauss_norm * apdf.total();;
}

double PoissonGaussianMES::pdf_derivs_mes(double x, double* derivs)
{
  const double x0 = ped_zero_dc_;
  const double xc = x-x0;

  accumulator apdf;
  accumulator aderivs[5];

  for(unsigned n=nmin_; n<nmax_; n++)
  {
    const double xn = xc-C3_[n];
    const double xn2 = SQR(xn);
    const double log_pn = C1_[n] - C2_[n]*xn2;
    const double pdf_n = std::exp(log_pn);

    apdf.accumulate(pdf_n);

    const double dlog0 = dC1_dF_[n];
    const double dlog1 = C2_[n]*2.0*xn;
    const double dlog2 = dC1_ds_[n] - dC2_ds_[n]*xn2;
    const double dlog3 = dC1_dg_[n] - dC2_dg_[n]*xn2 + dC3_dg_[n]*xn;
    const double dlog4 = dC1_db_[n] - dC2_db_[n]*xn2;
    
    aderivs[0].accumulate(pdf_n*dlog0);
    aderivs[1].accumulate(pdf_n*dlog1);
    aderivs[2].accumulate(pdf_n*dlog2);
    aderivs[3].accumulate(pdf_n*dlog3);
    aderivs[4].accumulate(pdf_n*dlog4);
  }

  derivs[0] = c_gauss_norm * aderivs[0].total();
  derivs[1] = c_gauss_norm * aderivs[1].total();
  derivs[2] = c_gauss_norm * aderivs[2].total();
  derivs[3] = c_gauss_norm * aderivs[3].total();
  derivs[4] = c_gauss_norm * aderivs[4].total();
  return c_gauss_norm * apdf.total();
}

double PoissonGaussianMES::
pdf_derivs_hessian_mes(double x, double* derivs, double* hessian)
{
  if(!hessian_elements_good_)calc_cached_vars(true);

  const double x0 = ped_zero_dc_;
  const double xc = x-x0;
      
  accumulator apdf;
  accumulator aderivs[5];
  accumulator ahessian[15];

  for(unsigned n=nmin_; n<nmax_; n++)
  {
    const double xn = xc-C3_[n];
    const double xn2 = SQR(xn);
    const double log_pn = C1_[n] - C2_[n]*xn2;
    const double pdf_n = std::exp(log_pn);

    apdf.accumulate(pdf_n);

    const double dln0 = dC1_dF_[n];
    const double dln1 = C2_[n]*2.0*xn;
    const double dln2 = dC1_ds_[n] - dC2_ds_[n]*xn2;
    const double dln3 = dC1_dg_[n] - dC2_dg_[n]*xn2 + dC3_dg_[n]*xn;
    const double dln4 = dC1_db_[n] - dC2_db_[n]*xn2;

    aderivs[0].accumulate(pdf_n*dln0);
    aderivs[1].accumulate(pdf_n*dln1);
    aderivs[2].accumulate(pdf_n*dln2);
    aderivs[3].accumulate(pdf_n*dln3);
    aderivs[4].accumulate(pdf_n*dln4);

    const double d2ln00 = d2C1_dF2_[n];

    const double d2ln11 = C2_[n]*-2.0;
    const double d2ln12 = dC2_ds_[n]*2.0*xn;
    const double d2ln13 = dC2_dg_[n]*2.0*xn - dC3_dg_[n];
    const double d2ln14 = dC2_db_[n]*2.0*xn;

    const double d2ln22 = d2C1_ds2_[n] - d2C2_ds2_[n]*xn2;
    const double d2ln23 = d2C1_dsdg_[n] - d2C2_dsdg_[n]*xn2 + d2C3_dsdg_[n]*xn;
    const double d2ln24 = d2C1_dsdb_[n] - d2C2_dsdb_[n]*xn2;

    // caution how extra terms are included in d2C1_dg2_ and dC3_dg2_
    const double d2ln33 = d2C1_dg2_[n] - d2C2_dg2_[n]*xn2 + d2C3_dg2_[n]*xn;
    const double d2ln34 = d2C1_dgdb_[n] - d2C2_dgdb_[n]*xn2 + d2C3_dgdb_[n]*xn;
    
    const double d2ln44 = d2C1_db2_[n] - d2C2_db2_[n]*xn2;
    
    ahessian[0].accumulate(pdf_n*(dln0*dln0 + d2ln00));
    ahessian[1].accumulate(pdf_n*dln0*dln1);
    ahessian[2].accumulate(pdf_n*dln0*dln2);
    ahessian[3].accumulate(pdf_n*dln0*dln3);
    ahessian[4].accumulate(pdf_n*dln0*dln4);

    ahessian[5].accumulate(pdf_n*(dln1*dln1 + d2ln11));
    ahessian[6].accumulate(pdf_n*(dln1*dln2 + d2ln12));
    ahessian[7].accumulate(pdf_n*(dln1*dln3 + d2ln13));
    ahessian[8].accumulate(pdf_n*(dln1*dln4 + d2ln14));
    
    ahessian[9].accumulate(pdf_n*(dln2*dln2 + d2ln22));
    ahessian[10].accumulate(pdf_n*(dln2*dln3 + d2ln23));
    ahessian[11].accumulate(pdf_n*(dln2*dln4 + d2ln24));
    
    ahessian[12].accumulate(pdf_n*(dln3*dln3 + d2ln33));
    ahessian[13].accumulate(pdf_n*(dln3*dln4 + d2ln34));

    ahessian[14].accumulate(pdf_n*(dln4*dln4 + d2ln44));
  }

  derivs[0] = c_gauss_norm * aderivs[0].total();
  derivs[1] = c_gauss_norm * aderivs[1].total();
  derivs[2] = c_gauss_norm * aderivs[2].total();
  derivs[3] = c_gauss_norm * aderivs[3].total();
  derivs[4] = c_gauss_norm * aderivs[4].total();

  hessian[0] = c_gauss_norm * ahessian[0].total();
  hessian[1] = c_gauss_norm * ahessian[1].total();
  hessian[2] = c_gauss_norm * ahessian[2].total();
  hessian[3] = c_gauss_norm * ahessian[3].total();
  hessian[4] = c_gauss_norm * ahessian[4].total();

  hessian[5] = hessian[1];
  hessian[6] = c_gauss_norm * ahessian[5].total();
  hessian[7] = c_gauss_norm * ahessian[6].total();
  hessian[8] = c_gauss_norm * ahessian[7].total();
  hessian[9] = c_gauss_norm * ahessian[8].total();

  hessian[10] = hessian[2];
  hessian[11] = hessian[7];
  hessian[12] = c_gauss_norm * ahessian[9].total();
  hessian[13] = c_gauss_norm * ahessian[10].total();
  hessian[14] = c_gauss_norm * ahessian[11].total();
  
  hessian[15] = hessian[3];
  hessian[16] = hessian[8];
  hessian[17] = hessian[13];
  hessian[18] = c_gauss_norm * ahessian[12].total();
  hessian[19] = c_gauss_norm * ahessian[13].total();

  hessian[20] = hessian[4];
  hessian[21] = hessian[9];
  hessian[22] = hessian[14];
  hessian[23] = hessian[19];
  hessian[24] = c_gauss_norm * ahessian[14].total();
  
  return c_gauss_norm * apdf.total();  
}

double PoissonGaussianMES::pdf_ped(double x)
{
  const double xc = x - ped_zero_dc_;
  const double log_pn = -(log_s_ + C2_[0]*SQR(xc));
  return c_gauss_norm * std::exp(log_pn);
}

double PoissonGaussianMES::
pdf_derivs_ped(double x, double* derivs)
{
  const double xc = x - ped_zero_dc_;
  const double xc2 = SQR(xc);
  const double log_pn = -(log_s_ + C2_[0]*xc2);
  const double pdf = c_gauss_norm * std::exp(log_pn); 

  const double dln1 = 2.0*C2_[0]*xc;
  const double dln2 = -(1.0/ped_rms_dc_ + dC2_ds_[0]*xc2);

  derivs[0] = 0;
  derivs[1] = pdf * dln1;
  derivs[2] = pdf * dln2;
  derivs[3] = 0;
  derivs[4] = 0;
  
  return pdf;
}

double PoissonGaussianMES::
pdf_derivs_hessian_ped(double x, double* derivs, double* hessian)
{
  if(!hessian_elements_good_)calc_cached_vars(true);

  const double xc = x - ped_zero_dc_;
  const double xc2 = SQR(xc);
  const double log_pn = -(log_s_ + C2_[0]*xc2);
  const double pdf = c_gauss_norm * std::exp(log_pn); 

  const double dln1 = 2.0*C2_[0]*xc;
  const double dln2 = -(1.0/ped_rms_dc_ + dC2_ds_[0]*xc2);

  derivs[0] = 0;
  derivs[1] = pdf * dln1;
  derivs[2] = pdf * dln2;
  derivs[3] = 0;
  derivs[4] = 0;

  hessian[0] = hessian[1] = hessian[2] = hessian[3] = hessian[4] = 0;

  hessian[5] = 0;
  hessian[6] = pdf * (dln1*dln1 - 2.0*C2_[0]);
  hessian[7] = pdf * (dln1*dln2 + 2.0*dC2_ds_[0]*xc);
  hessian[8] = hessian[9] = 0;

  hessian[10] = 0;
  hessian[11] = hessian[7];
  hessian[12] = pdf * (dln2*dln2 + 1.0/SQR(ped_rms_dc_) - d2C2_ds2_[0]*xc2);
  hessian[13] = hessian[14] = 0;
  
  hessian[15] = hessian[16] = hessian[17] = hessian[18] = hessian[19] = 0;
  hessian[20] = hessian[21] = hessian[22] = hessian[23] = hessian[24] = 0;
  
  return pdf;
}

void PoissonGaussianMES::calc_cached_vars(bool calc_hessian)
{
  const double F = intensity_pe_;
  const double g = gain_dc_pe_;
  const double s = ped_rms_dc_;
  const double b = ses_rms_pe_;

  s2_ = SQR(s);
  b2_ = SQR(b);
  g2_ = SQR(g);
  log_s_ = std::log(s);
  
  const double var_n = b2_*g2_;
  const double var_s = s2_;

  calc_hessian |= force_calc_hessian_;
  if(calc_hessian && d2C1_dF2_.empty())
  {
    d2C1_dF2_.resize(nmax_);
    d2C1_ds2_.resize(nmax_);
    d2C1_dg2_.resize(nmax_);
    d2C1_db2_.resize(nmax_);
    d2C1_dsdg_.resize(nmax_);
    d2C1_dsdb_.resize(nmax_);
    d2C1_dgdb_.resize(nmax_);
    d2C2_ds2_.resize(nmax_);
    d2C2_dg2_.resize(nmax_);
    d2C2_db2_.resize(nmax_);
    d2C2_dsdg_.resize(nmax_);
    d2C2_dsdb_.resize(nmax_);
    d2C2_dgdb_.resize(nmax_);
    d2C3_dg2_.resize(nmax_);
    d2C3_dsdg_.resize(nmax_);
    d2C3_dgdb_.resize(nmax_);
  }
  
  for(unsigned n=0; n<nmax_; n++) // always loop from n=0 here
    {
      const double dbl_n = double(n);
      const double var = dbl_n*var_n + var_s;

#if 0
      accumulator acc;
      acc.accumulate(dbl_n*std::log(F));
      acc.accumulate(-lgamma(dbl_n+1.0));
      acc.accumulate(-F);
      acc.accumulate(-0.5*std::log(var));
      C1_[n] = acc.total();
#else 
      C1_[n] = dbl_n*std::log(F) - F - lgamma(dbl_n+1.0) - 0.5*std::log(var);
#endif
      C2_[n] = 0.5/var;
      C3_[n] = dbl_n*g;

      dC1_dF_[n] = dbl_n/F - 1.0;
      dC1_ds_[n] = -s/var;
      dC1_dg_[n] = -dbl_n*b2_*g/var;
      dC1_db_[n] = -dbl_n*g2_*b/var;

      dC2_ds_[n] = dC1_ds_[n]/var;
      dC2_dg_[n] = dC1_dg_[n]/var;
      dC2_db_[n] = dC1_db_[n]/var;

      dC3_dg_[n] = C2_[n]*2.0*dbl_n;

      if(!calc_hessian)
      {
        hessian_elements_good_ = false;
        continue;
      }
      
      d2C1_dF2_[n]  = -dbl_n/SQR(F);
      d2C1_ds2_[n]  = -(1.0/var + 2.0*s*dC2_ds_[n]);
      d2C1_dg2_[n]  = -dbl_n*b2_*(1.0/var + 2.0*g*dC2_dg_[n]);
      d2C1_db2_[n]  = -dbl_n*g2_*(1.0/var + 2.0*b*dC2_db_[n]);
      d2C1_dsdg_[n] = -2.0*s*dC2_dg_[n];
      d2C1_dsdb_[n] = -2.0*s*dC2_db_[n];
      d2C1_dgdb_[n] = -2.0*dbl_n*b*g*(1.0/var + b*dC2_db_[n]);

      d2C2_ds2_[n]  = d2C1_ds2_[n]/var + 2.0*dC1_ds_[n]*dC2_ds_[n];
      d2C2_dg2_[n]  = d2C1_dg2_[n]/var + 2.0*dC1_dg_[n]*dC2_dg_[n];
      d2C2_db2_[n]  = d2C1_db2_[n]/var + 2.0*dC1_db_[n]*dC2_db_[n];
      d2C2_dsdg_[n] = d2C1_dsdg_[n]/var + 2.0*dC1_ds_[n]*dC2_dg_[n];
      d2C2_dsdb_[n] = d2C1_dsdb_[n]/var + 2.0*dC1_ds_[n]*dC2_db_[n];
      d2C2_dgdb_[n] = d2C1_dgdb_[n]/var + 2.0*dC1_dg_[n]*dC2_db_[n];

      d2C3_dg2_[n]  = dC2_dg_[n]*2.0*dbl_n;
      d2C3_dsdg_[n] = dC2_ds_[n]*2.0*dbl_n;
      d2C3_dgdb_[n] = dC2_db_[n]*2.0*dbl_n;

      // Add extra terms needed for "d2ln33"
      d2C3_dg2_[n] *= 2.0;
      d2C1_dg2_[n] -= dC3_dg_[n]*dbl_n;

      hessian_elements_good_ = true;
    }
}

// ============================================================================
//
// High accuracy PG model for low and high intensity regime using legacy code
//
// ============================================================================

PoissonGaussianMES_HighAccuracy::
PoissonGaussianMES_HighAccuracy(double tol):
    tol_(tol)
{
  // nothing to see here
}

PoissonGaussianMES_HighAccuracy::~PoissonGaussianMES_HighAccuracy()
{
  // nothing to see here
}

auto PoissonGaussianMES_HighAccuracy::parameters() ->
    std::vector<math::ParameterAxis>
{
  return { { "light_intensity", "PE", true, 0, false, 0 },
    { "ped_zero", "DC", false, 0, false, 0 },
    { "ped_width", "DC", true, 0, false, 0 },
    { "gain", "DC/PE", false, 0, false, 0 }, // can be negative
    { "ses_width", "PE", true, 0, false, 0 }
  };
}

std::vector<double> PoissonGaussianMES_HighAccuracy::parameter_values()
{
  return { intensity_pe_, ped_zero_dc_, ped_rms_dc_,
        gain_dc_pe_, ses_rms_pe_, };
}

void PoissonGaussianMES_HighAccuracy::
set_parameter_values(const double* values)
{
  assign_parameters(values, intensity_pe_, ped_zero_dc_, ped_rms_dc_,
                    gain_dc_pe_, ses_rms_pe_);
}

bool PoissonGaussianMES_HighAccuracy::can_calculate_parameter_derivs()
{
  return true;
}

double PoissonGaussianMES_HighAccuracy::pdf_mes(double x)
{
  return std::exp(PMTModelPG::logL_multi_pe(x, intensity_pe_, 1.0,
                                            SQR(ses_rms_pe_), gain_dc_pe_,
                                            SQR(ped_rms_dc_), ped_zero_dc_,
                                            tol_));
}

double PoissonGaussianMES_HighAccuracy::pdf_ped(double x)
{
  return std::exp(PMTModelPG::logL_ped(x, SQR(ped_rms_dc_), ped_zero_dc_));
}

double PoissonGaussianMES_HighAccuracy::
pdf_derivs_mes(double x, double* derivs)
{
  PMTModelPGDerivs log_derivs;
  double p =
      std::exp(PMTModelPG::logL_multi_pe_derivs(log_derivs, x,
           intensity_pe_, 1.0, SQR(ses_rms_pe_), gain_dc_pe_,
           SQR(ped_rms_dc_), ped_zero_dc_, tol_));

  derivs[0] = p * log_derivs.dfdF;
  derivs[1] = p * log_derivs.dfdx0;
  derivs[2] = p * 2.0 * ped_rms_dc_ * log_derivs.dfds2;
  derivs[3] = p * log_derivs.dfdg;
  derivs[4] = p * 2.0 * ses_rms_pe_ * log_derivs.dfdb2;

  return p;
}

double PoissonGaussianMES_HighAccuracy::
pdf_derivs_ped(double x, double* derivs)
{
  PMTModelPGDerivs log_derivs;
  double p =
      std::exp(PMTModelPG::logL_ped_derivs(log_derivs, x,
                                           SQR(ped_rms_dc_), ped_zero_dc_));

  derivs[0] = 0; // p * log_derivs.dfdF;
  derivs[1] = p * log_derivs.dfdx0;
  derivs[2] = p * 2.0 * ped_rms_dc_ * log_derivs.dfds2;
  derivs[3] = 0; // p * log_derivs.dfdg;
  derivs[4] = 0; // p * 2.0 * ses_rms_pe_ * log_derivs.dfdb2;

  return p;
}

// ============================================================================
//
// SPELikelihood -- calculate likelihood for optimizer
//
// ============================================================================

SPELikelihood::SPELikelihood(MultiElectronSpectrum& mes_model,
                             const SimpleHist& mes_data):
    MultiAxisFunction(), mes_model_(&mes_model),
    npar_(mes_model.parameters().size()),
    mes_data_(mes_data), has_ped_data_(false), ped_data_(1.0)
{
  // nothing to see here
}

SPELikelihood::SPELikelihood(MultiElectronSpectrum& mes_model,
                             const SimpleHist& mes_data,
                             const SimpleHist& ped_data):
    MultiAxisFunction(), mes_model_(&mes_model),
    npar_(mes_model.parameters().size()),
    mes_data_(mes_data), has_ped_data_(true), ped_data_(ped_data)
{
  // nothing to see here
}

SPELikelihood::~SPELikelihood()
{
  // nothing to see here
}

auto SPELikelihood::domain_axes() -> std::vector<math::DomainAxis>
{
  return mes_model_->parameters();
}

double SPELikelihood::value(const double* x)
{
  mes_model_->set_parameter_values(x);
  math::LikelihoodAccumulator acc;
  for(auto& ibin : mes_data_)
  {
    double pdf = mes_model_->pdf_mes(ibin.xval_center());
    acc.accumulate(std::log(pdf)*ibin.weight());
  }
  
  if(has_ped_data_)
    for(auto& ibin : ped_data_)
    {
      double pdf = mes_model_->pdf_ped(ibin.xval_center());
      acc.accumulate(std::log(pdf)*ibin.weight());
    }
  return -acc.total();
}

bool SPELikelihood::can_calculate_derivs()
{
  return mes_model_->can_calculate_parameter_derivs();
}

bool SPELikelihood::can_calculate_hessian()
{
  return mes_model_->can_calculate_parameter_hessian();
}

double SPELikelihood::value_and_derivs(const double* x, double* derivs)
{
  mes_model_->set_parameter_values(x);
  math::LikelihoodAccumulator acc;
  std::vector<math::LikelihoodAccumulator> derivs_acc(npar_);
  for(auto& ibin : mes_data_)
  {
    double pdf = mes_model_->pdf_derivs_mes(ibin.xval_center(), derivs);
    acc.accumulate(std::log(pdf)*ibin.weight());
    for(unsigned ipar=0;ipar<npar_;ipar++)
      derivs_acc[ipar].accumulate(derivs[ipar]/pdf*ibin.weight());
  }
    
  if(has_ped_data_)
    for(auto& ibin : ped_data_)
    {
      double pdf = mes_model_->pdf_derivs_ped(ibin.xval_center(), derivs);
      acc.accumulate(std::log(pdf)*ibin.weight());
      for(unsigned ipar=0;ipar<npar_;ipar++)
        derivs_acc[ipar].accumulate(derivs[ipar]/pdf*ibin.weight());
    }
  
  for(unsigned ipar=0;ipar<npar_;ipar++)
    derivs[ipar] = -derivs_acc[ipar].total();
  return -acc.total();
}

double SPELikelihood::
value_derivs_and_hessian(const double* x, double* derivs, double* hessian)
{
  mes_model_->set_parameter_values(x);
  math::LikelihoodAccumulator acc;
  std::vector<math::LikelihoodAccumulator> derivs_acc(npar_);
  std::vector<math::LikelihoodAccumulator> hessian_acc(npar_*(npar_+1)/2);
  for(auto& ibin : mes_data_)
  {
    double pdf =
        mes_model_->pdf_derivs_hessian_mes(ibin.xval_center(), derivs, hessian);
    acc.accumulate(std::log(pdf)*ibin.weight());
    for(unsigned ipar=0;ipar<npar_;ipar++)
      derivs_acc[ipar].accumulate(derivs[ipar]/pdf*ibin.weight());
    unsigned itri = 0;
    for(unsigned icol=0;icol<npar_;icol++)
      for(unsigned irow=icol;irow<npar_;irow++)
      {
        unsigned iel = icol*npar_+irow;
        double summand = (hessian[iel] - derivs[icol]*derivs[irow]/pdf)/pdf;
        hessian_acc[itri++].accumulate(summand*ibin.weight());
      }
  }
    
  if(has_ped_data_)
    for(auto& ibin : ped_data_)
    {
      double pdf =
          mes_model_->pdf_derivs_hessian_ped(ibin.xval_center(), derivs, hessian);
      acc.accumulate(std::log(pdf)*ibin.weight());
      for(unsigned ipar=0;ipar<npar_;ipar++)
        derivs_acc[ipar].accumulate(derivs[ipar]/pdf*ibin.weight());
      unsigned iacc = 0;
      for(unsigned icol=0;icol<npar_;icol++)
        for(unsigned irow=icol;irow<npar_;irow++)
        {
          unsigned iel = icol*npar_+irow;
          double summand = (hessian[iel] - derivs[icol]*derivs[irow]/pdf)/pdf;
          hessian_acc[iacc++].accumulate(summand*ibin.weight());
        }
    }
  
  for(unsigned ipar=0;ipar<npar_;ipar++)
    derivs[ipar] = -derivs_acc[ipar].total();

  unsigned iacc = 0;
  for(unsigned icol=0;icol<npar_;icol++)
  {
    for(unsigned irow=0;irow<icol;irow++)
    {
      unsigned iel = icol*npar_+irow;
      unsigned iel_trans = irow*npar_+icol;
      hessian[iel] = hessian[iel_trans];
    }
    for(unsigned irow=icol;irow<npar_;irow++)
    {
      unsigned iel = icol*npar_+irow;
      hessian[iel] = -hessian_acc[iacc++].total();
    }
  }
  
  return -acc.total();
}
