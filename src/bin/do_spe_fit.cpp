/*

   calin/bin/do_spe_fit.cpp -- Stephen Fegan -- 2015-03-15

   OBSOLETE driver program to do SPE gaussian fits - use Python now

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

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
  hist.insert(xval0,wval0);
  hist.insert(xval1,wval1);

  *stream >> xval1 >> wval1;
  while(*stream)
  {
    hist.insert(xval1,wval1);
    *stream >> xval1 >> wval1;
  }
  
  delete fstream;

  PoissonGaussianMES mes_model(20);
  SPELikelihood like(mes_model, hist);

  NLOptOptimizer opt("LD_LBFGS", &like);
  Eigen::VectorXd x(5);
  x <<  atof(argv[0]),atof(argv[1]),atof(argv[2]),
      atof(argv[3]),atof(argv[4]);
  std::cout << "START " << like.value(x) << ' '
            << x(0) << ' ' << x(1) << ' ' << x(2) << ' ' << x(3) << ' ' << x(4)
            << '\n';

  double f_val;
  opt.set_abs_tolerance(0.001);
  opt.set_initial_values(x);
  opt.minimize(x, f_val);

  std::cout << x(0) << ' ' << x(1) << ' ' << x(2) << ' '
            << x(3) << ' ' << x(4) << '\n';
  
  Eigen::MatrixXd hessian;
  Eigen::VectorXd gradient;
  like.value_gradient_and_hessian(x, gradient, hessian);
  Eigen::MatrixXd err_mat = hessian.inverse();
  std::cout << err_mat << '\n';

  std::cout << std::sqrt(0.5*err_mat(0,0)) << ' '
            << std::sqrt(0.5*err_mat(1,1)) << ' '
            << std::sqrt(0.5*err_mat(2,2)) << ' '
            << std::sqrt(0.5*err_mat(3,3)) << ' '
            << std::sqrt(0.5*err_mat(4,4)) << '\n';
  
}
