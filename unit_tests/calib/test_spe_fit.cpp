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

TEST(PoissonGaussianMES, SetAndRecallParameters) {
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

TEST(PoissonGaussianMES, PDFEqualityWithLegacyCode_Ped) {
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

TEST(PoissonGaussianMES, PDFEqualityWithLegacyCode_MES) {
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
                       
} // anonymous namespave



TEST(PoissonGaussianMES, GradientTest_PED)
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

TEST(PoissonGaussianMES, GradientTest_MES)
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

TEST(PoissonGaussianMES_HighAccuracy, GradientTest_PED)
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

TEST(PoissonGaussianMES_HighAccuracy, GradientTest_MES)
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

#if 0

TEST(PoissonGaussianMES_HighAccuracy, Dervative_PedRMS_MES) {
  deriv_test(make_pgha_mes,
             &MultiElectronSpectrum::pdf_mes,
             &MultiElectronSpectrum::pdf_gradient_mes,
             { 1.0, 0.100000, 0.2, 1.0, 0.45 }, 2, 1e-8, 1e-4, -1.0, 20.0);
}

TEST(PoissonGaussianMES_HighAccuracy, Dervative_Gain_MES) {
  deriv_test(make_pgha_mes,
             &MultiElectronSpectrum::pdf_mes,
             &MultiElectronSpectrum::pdf_gradient_mes,
             { 1.0, 0.100000, 0.2, 1.0, 0.45 }, 3, 1e-8, 1e-4, -1.0, 20.0);
}

TEST(PoissonGaussianMES_HighAccuracy, Dervative_SESRMS_MES) {
  deriv_test(make_pgha_mes,
             &MultiElectronSpectrum::pdf_mes,
             &MultiElectronSpectrum::pdf_gradient_mes,
             { 1.0, 0.100000, 0.2, 1.0, 0.45 }, 4, 1e-8, 1e-5, -1.0, 20.0);
}

namespace {

void mes_hessian_test(double (PoissonGaussianMES::*der)(double, double*),
                      double (PoissonGaussianMES::*hes)(double, double*, double *),
                      double xmax = 10.0)
{
  PoissonGaussianMES mes(20);

  std::vector<double> deriv1(5);
  std::vector<double> deriv2(5);
  std::vector<double> derivh(5);
  std::vector<double> hessian(25);

  const double dx = 1e-8;

  for(double x=0; x<xmax; x+=0.1)
    for(unsigned ipar = 0; ipar<5; ipar++)
    {
      std::vector<double> pars { 1.123, 0.234, 0.2, 0.9, 0.33 };

      mes.set_parameter_values(&pars.front());
      double pdf = (mes.*der)(x,&deriv1.front());
      double pdfh = (mes.*hes)(x,&derivh.front(),&hessian.front());

      EXPECT_EQ(pdf,pdfh);
      EXPECT_EQ(deriv1,derivh);

      std::vector<double> pars1 { pars };
      pars1[ipar] -= dx;
      mes.set_parameter_values(&pars1.front());
      (mes.*der)(x,&deriv1.front());      

      pars1 = pars;
      pars1[ipar] += dx;
      mes.set_parameter_values(&pars1.front());
      (mes.*der)(x,&deriv2.front());
      
      for(unsigned jpar=ipar;jpar<5;jpar++)
      {
        double dpdx = (deriv2[jpar]-deriv1[jpar])/(2.0*dx);
        EXPECT_NEAR(dpdx, hessian[ipar*5+jpar], std::abs(dpdx)*1e-5);
      }
    }
}

} // anonymous namespace

TEST(PoissonGaussianMES, HessianTest_MES)
{
  mes_hessian_test(&PoissonGaussianMES::pdf_gradient_mes,
                   &PoissonGaussianMES::pdf_gradient_hessian_mes);
}

TEST(PoissonGaussianMES, HessianTest_PED)
{
  mes_hessian_test(&PoissonGaussianMES::pdf_gradient_ped,
                   &PoissonGaussianMES::pdf_gradient_hessian_ped, 4.0);
}

namespace {

void spe_deriv_test(unsigned ip, double dx, double nearness)
{
  auto mes_data = karkar_data();
  SimpleHist mes_hist(1.0);
  for(auto idata : mes_data)mes_hist.accumulate(idata);
  PoissonGaussianMES mes_model(20);
  SPELikelihood like(mes_model, mes_hist);
  std::vector<double> x { 0.55349034289601895, 3094.2718624743093,
        19.614139336940855, 89.181964780087668, 0.32388058781378032 };
  //std::vector<double> x { 0.6, 3095, 20, 91, 0.4 };
  std::vector<double> gradient(5);
  double lvalue = like.value(&x.front());
  for(unsigned ix=0; ix<100; ix++)
  {
    x[ip] += dx;
    double value = like.value_and_gradient(&x.front(), &gradient.front());
    EXPECT_NEAR(gradient[ip], (value-lvalue)/dx, (lvalue+value)*nearness);
    lvalue = value;
  }
}

}

TEST(SPELikelihood, Deriv_Intensity) {
  spe_deriv_test(0, 1e-4, 1e-5);
  //spe_deriv_test(0, 1e-8, 1e-10);
}

TEST(SPELikelihood, Deriv_PedZero) {
  spe_deriv_test(1, 1e-2, 1e-3);
}

TEST(SPELikelihood, Deriv_PedRMS) {
  spe_deriv_test(2, 1e-3, 1e-4);
}

TEST(SPELikelihood, Deriv_Gain) {
  spe_deriv_test(3, 1e-3, 1e-4);
}

TEST(SPELikelihood, Deriv_SESRMS) {
  spe_deriv_test(4, 1e-4, 1e-5);
}

TEST(SPELikelihood, Hessian_Test)
{
  auto mes_data = karkar_data();
  SimpleHist mes_hist(1.0);
  for(auto idata : mes_data)mes_hist.accumulate(idata);
  PoissonGaussianMES mes_model(20);
  SPELikelihood like(mes_model, mes_hist);
  std::vector<double> x0 { 0.55349034289601895, 3094.2718624743093,
        19.614139336940855, 89.181964780087668, 0.32388058781378032 };
  double dx = 1e-8;
  
  std::vector<double> deriv1(5);
  std::vector<double> deriv2(5);
  std::vector<double> derivh(5);
  std::vector<double> hessian(25);

  for(unsigned ipar = 0; ipar<5; ipar++)
  {
    double pdf = like.value_and_gradient(&x0.front(),&deriv1.front());
    double pdfh = like.value_gradient_and_hessian(&x0.front(),&derivh.front(),
                                                &hessian.front());

    EXPECT_EQ(pdf,pdfh);
    EXPECT_EQ(deriv1,derivh);

    std::vector<double> x1 { x0 };
    x1[ipar] = x0[ipar] - dx;
    like.value_and_gradient(&x1.front(),&deriv1.front());

    x1[ipar] = x0[ipar] + dx;
    like.value_and_gradient(&x1.front(),&deriv2.front());

    for(unsigned jpar=ipar;jpar<5;jpar++)
    {
      double dfdx = (deriv2[jpar]-deriv1[jpar])/(2.0*dx);
      EXPECT_NEAR(dfdx, hessian[ipar*5+jpar], std::abs(dfdx)*2e-5);
    }
  }
}


namespace {

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

TEST(SPELikelihood, Minimize_GSL_Simplex)
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

TEST(SPELikelihood, Minimize_GSL_BFGS2)
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

TEST(SPELikelihood, Minimize_NLOpt_LD_LBFGS)
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
  std::vector<double> x { 1.0, 3100.0, 20.0, 100.0, 0.45 };
  opt.set_initial_values(x);
  double f_val;
  opt.minimize(x, f_val);
  
  //EXPECT_EQ(status, GSL_SUCCESS);
  EXPECT_NEAR(x[0], 0.55349337, 0.0001);
  EXPECT_NEAR(x[1], 3094.2715, 0.01);
  EXPECT_NEAR(x[2], 19.6141970, 0.001);
  EXPECT_NEAR(x[3], 89.1810077, 0.01);
  EXPECT_NEAR(x[4], 0.32388838, 0.0001);  
 
  std::cout << std::fixed << std::setprecision(3);
  std::cout << x[0] << ' ' << x[1] << ' ' << x[2] << ' '
            << x[3] << ' ' << x[4] << '\n';

  Eigen::MatrixXd hessian_mat(5,5);
  std::vector<double> gradient(5);
  like.value_gradient_and_hessian(&x.front(), &gradient.front(), hessian_mat.data());
  Eigen::MatrixXd err_mat = hessian_mat.inverse();
  std::cout << std::scientific << std::setprecision(8) << err_mat << '\n';

  std::cout << std::sqrt(err_mat(0,0)) << ' '
            << std::sqrt(err_mat(1,1)) << ' '
            << std::sqrt(err_mat(2,2)) << ' '
            << std::sqrt(err_mat(3,3)) << ' '
            << std::sqrt(err_mat(4,4)) << '\n';

}

TEST(SPELikelihood, Minimize_Minuit75)
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

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
