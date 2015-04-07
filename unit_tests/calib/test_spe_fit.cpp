#include <fstream>

#include <gsl/gsl_multimin.h>
#include <gtest/gtest.h>
#include "Eigen/Dense"
#include "nlopt/nlopt.hpp"
#include "calib/spe_fit.hpp"
#include "karkar_data.hpp"
#include "math/optimizer.hpp"
#include "math/nlopt_optimizer.hpp"
#include "math/minuit75_optimizer.hpp"

using namespace calin::math;
using namespace calin::data;
using namespace calin::calib;
using namespace calin::unit_tests;

TEST(TestPoissonGaussianMES, SetAndRecallParameters) {
  PoissonGaussianMES pg_mes1;
  PoissonGaussianMES_HighAccuracy pg_mes2;

  Eigen::VectorXd params(5);
  params <<  1.0, 0.1, 0.2, 1.0, 0.45;
  pg_mes1.set_parameter_values(params);
  pg_mes2.set_parameter_values(params);

  Eigen::VectorXd params1 { pg_mes1.parameter_values() };
  Eigen::VectorXd params2 { pg_mes2.parameter_values() };

  EXPECT_EQ(params, params1);
  EXPECT_EQ(params, params2);
}

TEST(TestPoissonGaussianMES, PDFEqualityWithLegacyCode_Ped) {
  PoissonGaussianMES pg_mes1;
  PoissonGaussianMES_HighAccuracy pg_mes2;

  Eigen::VectorXd params(5);
  params <<  1.0, 0.1, 0.2, 1.0, 0.45;
  pg_mes1.set_parameter_values(params);
  pg_mes2.set_parameter_values(params);
  
  Eigen::VectorXd gradient(5);
  
  for(double x=-1.0;x<1.001;x+=0.01)
  {
    EXPECT_NEAR(pg_mes1.pdf_ped(x),pg_mes2.pdf_ped(x),pg_mes2.pdf_ped(x)*1e-8);
    EXPECT_NEAR(pg_mes1.pdf_ped(x),pg_mes1.pdf_gradient_ped(x,gradient),
                pg_mes1.pdf_ped(x)*1e-8);
    EXPECT_NEAR(pg_mes2.pdf_ped(x),pg_mes2.pdf_gradient_ped(x,gradient),
                pg_mes2.pdf_ped(x)*1e-8);
  }
}

TEST(TestPoissonGaussianMES, PDFEqualityWithLegacyCode_MES) {
  PoissonGaussianMES pg_mes1(20);
  PoissonGaussianMES_HighAccuracy pg_mes2;

  Eigen::VectorXd params(5);
  params <<  1.0, 0.1, 0.2, 1.0, 0.45;
  pg_mes1.set_parameter_values(params);
  pg_mes2.set_parameter_values(params);

  Eigen::VectorXd gradient(5);
  
  for(double x=-1.0;x<10.001;x+=0.01)
  {
    EXPECT_NEAR(pg_mes1.pdf_mes(x),pg_mes2.pdf_mes(x),pg_mes2.pdf_mes(x)*1e-8);
    EXPECT_NEAR(pg_mes1.pdf_mes(x),pg_mes1.pdf_gradient_mes(x,gradient),
                pg_mes1.pdf_mes(x)*1e-8);
    EXPECT_NEAR(pg_mes2.pdf_mes(x),pg_mes2.pdf_gradient_mes(x,gradient),
                pg_mes2.pdf_mes(x)*1e-8);
  }
}

namespace {

using function::ConstVecRef;
using function::VecRef;
using function::MatRef;

void
mes_gradient_test(MultiElectronSpectrum* mes,
                  double(MultiElectronSpectrum::*val_get_f)(double),
                  double(MultiElectronSpectrum::*grad_get_f)(double,VecRef) ,
        double(MultiElectronSpectrum::*hess_get_f)(double,VecRef,MatRef),
                  const std::vector<double> vp, const std::vector<double> vdp,
                  double xlo, double xhi, double dx)
{
  Eigen::VectorXd p(5);
  p << vp[0], vp[1], vp[2], vp[3], vp[4];
  Eigen::VectorXd dp(5);
  dp << vdp[0], vdp[1], vdp[2], vdp[3], vdp[4];
  for(double xval = xlo; xval<xhi; xval+=dx)
  {
    bool check_ok;
    Eigen::VectorXd err(5);
    function::ValGetter<MultiElectronSpectrum> val_get =
        [xval,val_get_f](MultiElectronSpectrum* mes) {
      return (mes->*val_get_f)(xval); };
    function::GradGetter<MultiElectronSpectrum>  grad_get =
        [xval,grad_get_f](MultiElectronSpectrum* mes, VecRef grad) {
      return (mes->*grad_get_f)(xval,grad); };
    function::HessGetter<MultiElectronSpectrum> hess_get =
      [xval,hess_get_f](MultiElectronSpectrum* mes, VecRef grad, MatRef hess) {
      return (mes->*hess_get_f)(xval,grad,hess); };

    check_ok = function::gradient_check_par(*mes, p, dp, err,
                                            val_get, grad_get, hess_get);

    EXPECT_TRUE(check_ok);
    for(unsigned ipar=0;ipar<5;ipar++)EXPECT_LE(err(ipar),0.5);
  }
}
                       
} // anonymous namespace

TEST(TestPoissonGaussianMES, GradientCheck_PED)
{
  double dp1 = 1e-7;
  PoissonGaussianMES mes(40);
  mes_gradient_test(&mes,
                    &MultiElectronSpectrum::pdf_ped,
                    &MultiElectronSpectrum::pdf_gradient_ped,
                    &MultiElectronSpectrum::pdf_gradient_hessian_ped,
                    { 1.123, 0.100000, 0.2, 1.321, 0.45 },
                    { dp1, dp1, dp1, dp1, dp1}, -1.0, 1.0, 0.1);
}

TEST(TestPoissonGaussianMES, GradientCheck_MES)
{
  double dp1 = 1e-7;
  PoissonGaussianMES mes(40);
  mes_gradient_test(&mes,
                    &MultiElectronSpectrum::pdf_mes,
                    &MultiElectronSpectrum::pdf_gradient_mes,
                    &MultiElectronSpectrum::pdf_gradient_hessian_mes,
                    { 1.123, 0.100000, 0.2, 1.321, 0.45 },
                    { dp1, dp1, dp1, dp1, dp1}, -1.0, 10.0, 0.1);
}

TEST(TestPoissonGaussianMES_HighAccuracy, GradientCheck_PED)
{
  double dp1 = 1e-7;
  PoissonGaussianMES_HighAccuracy mes(1e-20);
  mes_gradient_test(&mes,
                    &MultiElectronSpectrum::pdf_ped,
                    &MultiElectronSpectrum::pdf_gradient_ped,
                    &MultiElectronSpectrum::pdf_gradient_hessian_ped,
                    { 1.123, 0.100000, 0.2, 1.321, 0.45 },
                    { dp1, dp1, dp1, dp1, dp1}, -1.0, 1.0, 0.1);
}

TEST(TestPoissonGaussianMES_HighAccuracy, GradientCheck_MES)
{
  double dp1 = 1e-7;
  PoissonGaussianMES_HighAccuracy mes(1e-20);
  mes_gradient_test(&mes,
                    &MultiElectronSpectrum::pdf_mes,
                    &MultiElectronSpectrum::pdf_gradient_mes,
                    &MultiElectronSpectrum::pdf_gradient_hessian_mes,
                    { 1.123, 0.100000, 0.2, 1.321, 0.45 },
                    { dp1, dp1, dp1, dp1, dp1}, -1.0, 10.0, 0.1);
}

namespace {

using function::ConstVecRef;
using function::VecRef;
using function::MatRef;

void
mes_hessian_test(MultiElectronSpectrum* mes,
                 double(MultiElectronSpectrum::*val_get_f)(double),
                 double(MultiElectronSpectrum::*grad_get_f)(double,VecRef) ,
        double(MultiElectronSpectrum::*hess_get_f)(double,VecRef,MatRef),
                 const std::vector<double> vp, const std::vector<double> vdp,
                 double xlo, double xhi, double dx)
{
  Eigen::VectorXd p(5);
  p << vp[0], vp[1], vp[2], vp[3], vp[4];
  Eigen::VectorXd dp(5);
  dp << vdp[0], vdp[1], vdp[2], vdp[3], vdp[4];
  for(double xval = xlo; xval<xhi; xval+=dx)
  {
    bool check_ok;
    Eigen::MatrixXd good(5,5);
    function::ValGetter<MultiElectronSpectrum> val_get =
        [xval,val_get_f](MultiElectronSpectrum* mes) {
      return (mes->*val_get_f)(xval); };
    function::GradGetter<MultiElectronSpectrum>  grad_get =
        [xval,grad_get_f](MultiElectronSpectrum* mes, VecRef grad) {
      return (mes->*grad_get_f)(xval,grad); };
    function::HessGetter<MultiElectronSpectrum> hess_get =
      [xval,hess_get_f](MultiElectronSpectrum* mes, VecRef grad, MatRef hess) {
      return (mes->*hess_get_f)(xval,grad,hess); };

    check_ok = function::hessian_check_par(*mes, p, dp, good,
                                           val_get, grad_get, hess_get);

    EXPECT_TRUE(check_ok);
    for(unsigned ipar=0;ipar<5;ipar++)
      for(unsigned jpar=0;jpar<5;jpar++)EXPECT_LE(good(ipar,jpar),0.5);
  }
}
                       
} // anonymous namespace

TEST(TestPoissonGaussianMES, HessianCheck_PED)
{
  double dp1 = 1e-7;
  PoissonGaussianMES mes(40);
  mes_hessian_test(&mes,
                   &MultiElectronSpectrum::pdf_ped,
                   &MultiElectronSpectrum::pdf_gradient_ped,
                   &MultiElectronSpectrum::pdf_gradient_hessian_ped,
                   { 1.123, 0.100000, 0.2, 1.321, 0.45 },
                   { dp1, dp1, dp1, dp1, dp1}, -1.0, 1.0, 0.1);
}

TEST(TestPoissonGaussianMES, HessianCheck_MES)
{
  double dp1 = 1e-7;
  PoissonGaussianMES mes(40);
  mes_hessian_test(&mes,
                   &MultiElectronSpectrum::pdf_mes,
                   &MultiElectronSpectrum::pdf_gradient_mes,
                   &MultiElectronSpectrum::pdf_gradient_hessian_mes,
                   { 1.123, 0.100000, 0.2, 1.321, 0.45 },
                   { dp1, dp1, dp1, dp1, dp1}, -1.0, 10.0, 0.1);
}

TEST(TestSPELikelihood, KarkarPG_GradientCheck) {
  auto mes_data = karkar_data();
  SimpleHist mes_hist(1.0);
  for(auto idata : mes_data)mes_hist.accumulate(idata);
  PoissonGaussianMES mes_model(20);
  SPELikelihood like(mes_model, mes_hist);
  Eigen::VectorXd x(5);
  x << 0.55349034289601895, 3094.2718624743093,
      19.614139336940855, 89.181964780087668, 0.32388058781378032;
  Eigen::VectorXd dx(5);
  dx << 1e-7, 1e-7, 1e-7, 1e-7, 1e-7;
  Eigen::VectorXd good(5);
  
  bool check_ok = function::gradient_check(like, x, dx, good);
  EXPECT_TRUE(check_ok);
  for(unsigned ipar=0;ipar<5;ipar++)EXPECT_LE(good(ipar),0.5);
}

TEST(TestSPELikelihood, KarkarPG_HessianCheck) {
  auto mes_data = karkar_data();
  SimpleHist mes_hist(1.0);
  for(auto idata : mes_data)mes_hist.accumulate(idata);
  PoissonGaussianMES mes_model(20);
  SPELikelihood like(mes_model, mes_hist);
  Eigen::VectorXd x(5);
  x << 0.55349034289601895, 3094.2718624743093,
      19.614139336940855, 89.181964780087668, 0.32388058781378032;
  Eigen::VectorXd dx(5);
  dx << 1e-7, 1e-7, 1e-7, 1e-7, 1e-7;
  Eigen::MatrixXd good(5,5);
  
  bool check_ok = function::hessian_check(like, x, dx, good);
  EXPECT_TRUE(check_ok);
  for(unsigned ipar=0;ipar<5;ipar++)
    for(unsigned jpar=0;jpar<5;jpar++)EXPECT_LE(good(ipar,jpar),0.5);
}

#if 0

double my_f(const gsl_vector * x, void * params)
{
  SPELikelihood* like = static_cast<SPELikelihood*>(params);
  double value = like->value(x->data);
#if 0
  std::cout << value << ' '
            << x->data[0] << ' '
            << x->data[1] << ' '
            << x->data[2] << ' '
            << x->data[3] << ' '
            << x->data[4] << '\n';
#endif
  return value;
}

void my_df(const gsl_vector * x, void * params, gsl_vector* g)
{
  SPELikelihood* like = static_cast<SPELikelihood*>(params);
  like->value_and_gradient(x->data, g->data);
#if 0
  std::cout << g->data[0] << ' '
            << g->data[1] << ' '
            << g->data[2] << ' '
            << g->data[3] << ' '
            << g->data[4] << ' '
            << x->data[0] << ' '
            << x->data[1] << ' '
            << x->data[2] << ' '
            << x->data[3] << ' '
            << x->data[4] << '\n';
#endif
}

void my_fdf(const gsl_vector * x, void * params, double* value, gsl_vector* g)
{
  SPELikelihood* like = static_cast<SPELikelihood*>(params);
  *value = like->value_and_gradient(x->data, g->data);
#if 0
  std::cout << *value << ' '
            << g->data[0] << ' '
            << g->data[1] << ' '
            << g->data[2] << ' '
            << g->data[3] << ' '
            << g->data[4] << ' '
            << x->data[0] << ' '
            << x->data[1] << ' '
            << x->data[2] << ' '
            << x->data[3] << ' '
            << x->data[4] << '\n';
#endif
}

}

TEST(TestSPELikelihood, Minimize_GSL_Simplex)
{
  auto mes_data = karkar_data();
  SimpleHist mes_hist(1.0);
  for(auto idata : mes_data)mes_hist.accumulate(idata);
  PoissonGaussianMES mes_model(20);
  SPELikelihood like(mes_model, mes_hist);

  const gsl_multimin_fminimizer_type *T = 
    gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss,*x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  /* Starting point */
  x = gsl_vector_alloc (5);
  gsl_vector_set (x, 0, 1.0);
  gsl_vector_set (x, 1, 3100.0);
  gsl_vector_set (x, 2, 20.0);
  gsl_vector_set (x, 3, 100.0);
  gsl_vector_set (x, 4, 0.45);

  /* Step size */
  ss = gsl_vector_alloc (5);
  gsl_vector_set (ss, 0, 0.1);
  gsl_vector_set (ss, 1, 10.0);
  gsl_vector_set (ss, 2, 1.0);
  gsl_vector_set (ss, 3, 10.0);
  gsl_vector_set (ss, 4, 0.05);

  /* Initialize method and iterate */
  minex_func.n = 5;
  minex_func.f = my_f;
  minex_func.params = &like;

  s = gsl_multimin_fminimizer_alloc (T, 5);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      
      if (status) 
        break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-3);

#if 0
      if (status == GSL_SUCCESS)
        {
          printf ("converged to minimum at\n");
        }

      printf ("%5lu %10.8f %10.4f %10.7f %10.7f %10.8f f() = %7.3f size = %.3f\n", 
              iter,
>              gsl_vector_get (s->x, 0), 
              gsl_vector_get (s->x, 1), 
              gsl_vector_get (s->x, 2), 
              gsl_vector_get (s->x, 3), 
              gsl_vector_get (s->x, 4), 
              s->fval, size);
#endif
    }
  while (status == GSL_CONTINUE && iter < 10000);

  EXPECT_EQ(status, GSL_SUCCESS);
  EXPECT_NEAR(s->x->data[0], 0.55349337, 0.0001);
  EXPECT_NEAR(s->x->data[1], 3094.2715, 0.01);
  EXPECT_NEAR(s->x->data[2], 19.6141970, 0.001);
  EXPECT_NEAR(s->x->data[3], 89.1810077, 0.01);
  EXPECT_NEAR(s->x->data[4], 0.32388838, 0.0001);
  
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
}

TEST(TestSPELikelihood, Minimize_GSL_BFGS2)
{
  auto mes_data = karkar_data();
  SimpleHist mes_hist(1.0);
  for(auto idata : mes_data)mes_hist.accumulate(idata);
  PoissonGaussianMES mes_model(20);
  //PoissonGaussianMES_HighAccuracy mes_model;
  SPELikelihood like(mes_model, mes_hist);

  const gsl_multimin_fdfminimizer_type *T = nullptr;
  gsl_multimin_fdfminimizer *s = nullptr;

  gsl_vector *ss, *x;
  gsl_multimin_function_fdf minex_func;

  size_t iter = 0;
  int status;
  double size;

  /* Starting point */
  x = gsl_vector_alloc (5);
  gsl_vector_set (x, 0, 1.0);
  gsl_vector_set (x, 1, 3100.0);
  gsl_vector_set (x, 2, 20.0);
  gsl_vector_set (x, 3, 100.0);
  gsl_vector_set (x, 4, 0.45);

  /* Initialize method and iterate */
  minex_func.n = 5;
  minex_func.f = my_f;
  minex_func.df = my_df;
  minex_func.fdf = my_fdf;
  minex_func.params = &like;

  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  //T = gsl_multimin_fdfminimizer_conjugate_fr;
  //T = gsl_multimin_fdfminimizer_steepest_descent;
  s = gsl_multimin_fdfminimizer_alloc (T, 5);

  gsl_multimin_fdfminimizer_set (s, &minex_func, x, .01, .01);

  do
  {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate (s);
    
    if (status)
      break;
    
    status = gsl_multimin_test_gradient (s->gradient, 0.1);

#if 0
    if (status == GSL_SUCCESS)
      printf ("Minimum found at:\n");
    
    printf ("%5lu %10.8f %10.4f %10.7f %10.7f %10.8f f() = %7.3f\n", 
            iter,
            gsl_vector_get (s->x, 0), 
            gsl_vector_get (s->x, 1), 
            gsl_vector_get (s->x, 2), 
            gsl_vector_get (s->x, 3), 
            gsl_vector_get (s->x, 4), 
            s->f);
#endif
  }
  while (status == GSL_CONTINUE && iter < 10000);

  EXPECT_EQ(status, GSL_SUCCESS);
  EXPECT_NEAR(s->x->data[0], 0.55349337, 0.0001);
  EXPECT_NEAR(s->x->data[1], 3094.2715, 0.01);
  EXPECT_NEAR(s->x->data[2], 19.6141970, 0.001);
  EXPECT_NEAR(s->x->data[3], 89.1810077, 0.01);
  EXPECT_NEAR(s->x->data[4], 0.32388838, 0.0001);
  
  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);
}
#endif

TEST(TestSPELikelihood, Optimize_NLOpt)
{
  auto mes_data = karkar_data();
  SimpleHist mes_hist(1.0);
  for(auto idata : mes_data)mes_hist.accumulate(idata);
  PoissonGaussianMES mes_model(20);
  SPELikelihood like(mes_model, mes_hist);

  //optimizer::NLOptOptimizer opt(nlopt::LN_SBPLX, &like);
  optimizer::NLOptOptimizer opt(nlopt::LD_LBFGS, &like);
  opt.set_scale({0.1,0.1,1.0,1.0,0.05});
  opt.set_verbosity_level(optimizer::OptimizerVerbosityLevel::MAX);
  opt.set_abs_tolerance(0.0001);
  opt.set_initial_values({ 1.0, 3100.0, 20.0, 100.0, 0.45 });
  Eigen::VectorXd x_opt(5);
  double f_val;
  opt.minimize(x_opt, f_val);
  
  //EXPECT_EQ(status, GSL_SUCCESS);
  EXPECT_NEAR(x_opt[0], 0.55349337, 0.0001);
  EXPECT_NEAR(x_opt[1], 3094.2715, 0.01);
  EXPECT_NEAR(x_opt[2], 19.6141970, 0.001);
  EXPECT_NEAR(x_opt[3], 89.1810077, 0.01);
  EXPECT_NEAR(x_opt[4], 0.32388838, 0.0001);  
 
  std::cout << std::fixed << std::setprecision(3);
  std::cout << x_opt[0] << ' ' << x_opt[1] << ' ' << x_opt[2] << ' '
            << x_opt[3] << ' ' << x_opt[4] << '\n';

  Eigen::MatrixXd hessian_mat(5,5);
  Eigen::VectorXd gradient(5);
  like.value_gradient_and_hessian(x_opt, gradient, hessian_mat);
  Eigen::MatrixXd err_mat = hessian_mat.inverse();
  std::cout << std::scientific << std::setprecision(8) << err_mat << '\n';

  std::cout << std::sqrt(err_mat(0,0)) << ' '
            << std::sqrt(err_mat(1,1)) << ' '
            << std::sqrt(err_mat(2,2)) << ' '
            << std::sqrt(err_mat(3,3)) << ' '
            << std::sqrt(err_mat(4,4)) << '\n';

}

#if 0
TEST(TestSPELikelihood, Minimize_Minuit75)
{
  auto mes_data = karkar_data();
  SimpleHist mes_hist(1.0);
  for(auto idata : mes_data)mes_hist.accumulate(idata);
  PoissonGaussianMES mes_model(20);
  SPELikelihood like(mes_model, mes_hist);

  optimizer::Minuit75Optimizer opt(&like);
  Eigen::VectorXd xopt { 1.0, 3100.0, 20.0, 100.0, 0.45 };
  double fopt;
  opt.set_initial_values(xopt);
  opt.set_scale({ 0.001, 0.1, 0.01, 0.1, 0.001 });
  opt.minimize(xopt,fopt);

  #if 0
  //EXPECT_EQ(status, GSL_SUCCESS);
  EXPECT_NEAR(x[0], 0.55349337, 0.0001);
  EXPECT_NEAR(x[1], 3094.2715, 0.01);
  EXPECT_NEAR(x[2], 19.6141970, 0.001);
  EXPECT_NEAR(x[3], 89.1810077, 0.01);
  EXPECT_NEAR(x[4], 0.32388838, 0.0001);  

  std::cout << x[0] << ' ' << x[1] << ' ' << x[2] << ' '
            << x[3] << ' ' << x[4] << '\n';

  Eigen::MatrixXd hessian_mat(5,5);
  std::vector<double> gradient(5);
  like.value_gradient_and_hessian(&x.front(), &gradient.front(), hessian_mat.data());
  Eigen::MatrixXd err_mat = hessian_mat.inverse();
  std::cout << err_mat << '\n';

  std::cout << std::sqrt(err_mat(0,0)) << ' '
            << std::sqrt(err_mat(1,1)) << ' '
            << std::sqrt(err_mat(2,2)) << ' '
            << std::sqrt(err_mat(3,3)) << ' '
            << std::sqrt(err_mat(4,4)) << '\n';
#endif
}

#endif

TEST(TestGeneralPoissonMES_GG, SetAndRecallParameters) {
  pdf_1d::GaussianPDF ped;
  pdf_1d::LimitedGaussianPDF ses(0,std::numeric_limits<double>::infinity());
  GeneralPoissonMES mes(0, 1, 1025, &ses, &ped);
  EXPECT_EQ(mes.num_parameters(), 5);
  Eigen::VectorXd p(5);
  p << 1.0, 100.0, 20.0, 100.0, 35.0;
  mes.set_parameter_values(p);
  EXPECT_EQ(mes.parameter_values(), p);
  std::ofstream file("spec.dat");
  std::vector<double> mes_spec = mes.multi_electron_spectrum();
  std::vector<double> ped_spec = mes.pedestal_spectrum();
  std::vector<double> one_es_spec = mes.n_electron_spectrum(1);
  std::vector<double> two_es_spec = mes.n_electron_spectrum(2);
  for(unsigned i=0;i<mes_spec.size();i++)
    file << mes_spec[i] << ' ' << ped_spec[i] << ' '
         << one_es_spec[i] << ' ' << two_es_spec[i] << '\n';
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
