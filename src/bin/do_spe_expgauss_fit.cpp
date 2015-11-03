/*

   calin/bin/do_spe_expgauss_fit.cpp -- Stephen Fegan -- 2015-04-14

   OBSELETE driver program to do Exp-Gauss fit - use Python now

   Copyright 2015, Stephen Fegan <sfegan@gmail.com>

   This file is part of "calin"
   
   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.
    
   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>

#include <Eigen/Dense>
#include <nlopt/nlopt.hpp>
#include <calib/spe_fit.hpp>
#include <math/nlopt_optimizer.hpp>

using namespace calin::math;
using namespace calin::math::histogram;
using namespace calin::math::optimizer;
using namespace calin::calib;
using namespace calin::calib::spe_fit;

int main(int argc, char** argv)
{
  constexpr double inf = std::numeric_limits<double>::infinity();
  
  std::string progname(*argv);
  argc--, argv++;

  if(argc!=8)
  {
    std::cout
        << "Usage: " << progname << " filename lambda ped_zero ped_rms exp_prob exp_scale gauss_mean gauss_rms\n";
    exit(EXIT_FAILURE);
  }

  std::istream* stream = &std::cin;
  std::ifstream* fstream = nullptr;
  
  if(std::string(*argv) != "-")
    stream = fstream = new std::ifstream(*argv);
  argc--, argv++;

  SimpleHist mes_hist(1.0);
  double xval;
  *stream >> xval;
  while(*stream)
  {
    mes_hist.insert(xval);
    *stream >> xval;
  }
  
  delete fstream;

  pdf_1d::GaussianPDF ped;
  pdf_1d::LimitedExponentialPDF exp_pdf(0,inf,mes_hist.dxval());
  exp_pdf.limit_scale(0.1, inf);
  pdf_1d::LimitedGaussianPDF gauss_pdf(0,inf);
  pdf_1d::TwoComponentPDF ses(&exp_pdf, "exp", &gauss_pdf, "gauss");
  GeneralPoissonMES mes_model(mes_hist.xval_left(0),
                              mes_hist.dxval(),
                              mes_hist.size(), &ses, &ped);
  SPELikelihood like(mes_model, mes_hist);

  NLOptOptimizer opt("LD_LBFGS", &like);
  opt.set_verbosity_level(optimizer::OptimizerVerbosityLevel::MAX);
  Eigen::VectorXd x(7);
  x <<  atof(argv[0]),atof(argv[1]),atof(argv[2]),
      atof(argv[3]),atof(argv[4]),atof(argv[5]),atof(argv[6]);
  std::cout << "START: " << like.value(x) << ' '
            << x(0) << ' ' << x(1) << ' ' << x(2) << ' ' << x(3) << ' ' << x(4)
            << x(5) << ' ' << x(6) << "\n\n";

  double f_val;
  opt.set_abs_tolerance(0.001);
  opt.set_initial_values(x);
  opt.minimize(x, f_val);

  std::cout << "\nFINAL: " << x(0) << ' ' << x(1) << ' ' << x(2) << ' '
            << x(3) << ' ' << x(4) << ' ' << x(5) << ' ' << x(6) << "\n\n";

  Eigen::MatrixXd err_mat(7,7);
  opt.calc_error_matrix(err_mat);

  std::cout << err_mat << "\n\n";

  std::cout << std::sqrt(0.5*err_mat(0,0)) << ' '
            << std::sqrt(0.5*err_mat(1,1)) << ' '
            << std::sqrt(0.5*err_mat(2,2)) << ' '
            << std::sqrt(0.5*err_mat(3,3)) << ' '
            << std::sqrt(0.5*err_mat(4,4)) << ' '
            << std::sqrt(0.5*err_mat(5,5)) << ' '
            << std::sqrt(0.5*err_mat(6,6)) << '\n';

  std::ofstream file("spec.dat");
  std::vector<double> mes_spec = mes_model.multi_electron_spectrum();
  std::vector<double> ped_spec = mes_model.pedestal_spectrum();
  std::vector<double> one_es_spec = mes_model.n_electron_spectrum(1);
  std::vector<double> two_es_spec = mes_model.n_electron_spectrum(2);
  std::vector<double> three_es_spec = mes_model.n_electron_spectrum(3);
  std::vector<double> zero_es_cpt = mes_model.mes_n_electron_cpt(0);
  std::vector<double> one_es_cpt = mes_model.mes_n_electron_cpt(1);
  std::vector<double> two_es_cpt = mes_model.mes_n_electron_cpt(2);
  std::vector<double> three_es_cpt = mes_model.mes_n_electron_cpt(3);

  for(unsigned i=0;i<mes_spec.size();i++)
    file << mes_hist.xval_left(0)+mes_hist.dxval()*(0.5+i) << ' '
         << mes_spec[i] << ' '
         << ped_spec[i] << ' ' << one_es_spec[i] << ' '
         << two_es_spec[i] << ' ' << three_es_spec[i] << ' '
         << zero_es_cpt[i] << ' ' << one_es_cpt[i] << ' '
         << two_es_cpt[i] << ' ' << three_es_cpt[i] << ' ' << '\n';  
}
