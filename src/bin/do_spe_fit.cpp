#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>

#include <Eigen/Dense>
#include <nlopt/nlopt.hpp>
#include <calib/spe_fit.hpp>

using namespace calin::math;
using namespace calin::calib;

unsigned nlopt_f_count(bool increment = true, bool reset = false)
{
  static unsigned iter = 0;
  if(increment)iter++;
  if(reset)iter=0;
  return iter;
}

static double nlopt_f(unsigned n, const double* x, double* grad,
                      void* f_data)
{
  SPELikelihood* like = static_cast<SPELikelihood*>(f_data);
  double value = like->value_and_derivs(x, grad);
#if 1
  std::cout << nlopt_f_count() << ' '
            << std::setprecision(10) << value << ' '
            << x[0] << ' '
            << x[1] << ' '
            << x[2] << ' '
            << x[3] << ' '
  << x[4] << ' '
            << grad[0] << ' '
            << grad[1] << ' '
            << grad[2] << ' '
            << grad[3] << ' '
            << grad[4] << '\n';
#endif
  return value;
}

int main(int argc, char** argv)
{
  std::string progname(*argv);
  argc--, argv++;

  if(argc!=6)
  {
    std::cout
        << "Usage: " << progname << " filename lambda ped0 pedrms gain res\n";
    exit(EXIT_FAILURE);
  }

  std::istream* stream = &std::cin;
  std::ifstream* fstream = nullptr;
  
  if(std::string(*argv) != "-")
    stream = fstream = new std::ifstream(*argv);
  argc--, argv++;

  double xval0;
  double wval0;
  *stream >> xval0 >> wval0;
  double xval1;
  double wval1;
  *stream >> xval1 >> wval1;

  assert(*stream);

  double dx = xval1-xval0;
  double xalign = SimpleHist::xalign_to_center_value(xval0,dx);
  
  SimpleHist hist(dx, xalign);
  hist.accumulate(xval0,wval0);
  hist.accumulate(xval1,wval1);

  *stream >> xval1 >> wval1;
  while(*stream)
  {
    hist.accumulate(xval1,wval1);
    *stream >> xval1 >> wval1;
  }
  
  delete fstream;

  PoissonGaussianMES mes_model(20);
  SPELikelihood like(mes_model, hist);

  nlopt::opt opt(nlopt::LD_LBFGS, 5);
  opt.set_min_objective(nlopt_f, &like);
  std::vector<double> xlim_lo;
  std::vector<double> xlim_hi;
  bool has_xlim_lo = false;
  bool has_xlim_hi = false;
  for(const auto& ipar : like.domain_axes())
  {
    has_xlim_lo |= ipar.has_lo_bound;
    xlim_lo.push_back(ipar.has_lo_bound?ipar.lo_bound:-HUGE_VAL);
    has_xlim_hi |= ipar.has_hi_bound;
    xlim_hi.push_back(ipar.has_hi_bound?ipar.hi_bound:HUGE_VAL);
  }
  if(has_xlim_lo)opt.set_lower_bounds(xlim_lo);
  if(has_xlim_hi)opt.set_upper_bounds(xlim_hi);

  std::vector<double> x { atof(argv[0]),atof(argv[1]),atof(argv[2]),
        atof(argv[3]),atof(argv[4]) };
  std::cout << "START " << like.value(&x.front()) << ' ' << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << x[3] << ' ' << x[4] << '\n';
  double f_val;
  opt.set_ftol_abs(0.001);
  nlopt_f_count(false,/*reset=*/true);
  opt.optimize(x, f_val);

  std::cout << x[0] << ' ' << x[1] << ' ' << x[2] << ' '
            << x[3] << ' ' << x[4] << '\n';
  
  Eigen::MatrixXd hessian_mat;
  std::vector<double> derivs;
  like.value_derivs_and_hessian(x, derivs, hessian_mat);
  Eigen::MatrixXd err_mat = hessian_mat.inverse();
  std::cout << err_mat << '\n';

  std::cout << std::sqrt(0.5*err_mat(0,0)) << ' '
            << std::sqrt(0.5*err_mat(1,1)) << ' '
            << std::sqrt(0.5*err_mat(2,2)) << ' '
            << std::sqrt(0.5*err_mat(3,3)) << ' '
            << std::sqrt(0.5*err_mat(4,4)) << '\n';
  
}
