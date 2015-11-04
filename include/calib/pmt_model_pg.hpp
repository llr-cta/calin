// PRE-CALIN LEGACY CODE THAT DOES NOT CONFORM TO CONDING STANDARDS

// PMTModelPG.hpp - Stephen Fegan - sfegan@llr.in2p3.fr - December 2013

/* 

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

   This file is part of "calin"
   
   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.
    
   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include<vector>

namespace calin { namespace calib { namespace pmt_model_pg {

struct PMTModelPGDerivs
{
  PMTModelPGDerivs(): dfdF(), dfdar(), dfdb2(), dfdg(), dfds2(), dfdx0() { }

  double dfdF;
  double dfdar;
  double dfdb2;
  double dfdg;
  double dfds2;
  double dfdx0;

  PMTModelPGDerivs& operator+= (const PMTModelPGDerivs& o)
  {
    dfdF  += o.dfdF;
    dfdar += o.dfdar;
    dfdb2 += o.dfdb2;
    dfdg  += o.dfdg;
    dfds2 += o.dfds2;
    dfdx0 += o.dfdx0;
    return *this;
  }

  PMTModelPGDerivs& operator*= (double x)
  {
    dfdF  *= x;
    dfdar *= x;
    dfdb2 *= x;
    dfdg  *= x;
    dfds2 *= x;
    dfdx0 *= x;
    return *this;
  }
};

class PMTModelPG
{
public:
  static double logL_ped(double x, double s2, double x0 = 0);
  static double logL_ped_derivs(PMTModelPGDerivs &derivs,
				double x, double s2, double x0 = 0);
  static double logL_multi_pe(double x, double F, double a, double b2,
			      double g, double s2, double x0 = 0, 
			      const double tol = 1e-200);
  static double logL_multi_pe_derivs(PMTModelPGDerivs &derivs,
				     double x, double F, double a, double b2,
				     double g, double s2, double x0 = 0,
				     const double tol = 1e-200);
  static double logL_multi_pe_rg(double x, double F, double r, double b2, 
				 double g, double s2, double x0 = 0, 
				 const double tol = 1e-200);
  static double logL_multi_pe_rg_derivs(PMTModelPGDerivs &derivs,
					double x, double F, double r, double b2,
					double g, double s2, double x0 = 0,
					const double tol = 1e-200);

  static double sum_logL_ped(const std::vector<double>& x,
			     double s2, double x0 = 0);
  static double sum_logL_ped_derivs(PMTModelPGDerivs &derivs,
				    const std::vector<double>& x, 
				    double s2, double x0 = 0);
  static double sum_logL_multi_pe(const std::vector<double>& x, 
				  double F, double a, double b2,
				  double g, double s2, double x0 = 0, 
				  const double tol = 1e-200);
  static double sum_logL_multi_pe_derivs(PMTModelPGDerivs &derivs,
					 const std::vector<double>& x, 
					 double F, double a, double b2,
					 double g, double s2, double x0 = 0,
					 const double tol = 1e-200);
  static double sum_logL_multi_pe_rg(const std::vector<double>& x, 
				     double F, double r, double b2, 
				     double g, double s2, double x0 = 0, 
				     const double tol = 1e-200);
  static double sum_logL_multi_pe_rg_derivs(PMTModelPGDerivs &derivs,
					    const std::vector<double>& x, 
					    double F, double r, double b2,
					    double g, double s2, double x0 = 0,
					    const double tol = 1e-200);

  using HistValues = std::pair<double,double>;
  
  static double sum_logL_ped(const std::vector<HistValues>& x,
			     double s2, double x0 = 0);
  static double sum_logL_ped_derivs(PMTModelPGDerivs &derivs,
				    const std::vector<HistValues>& x, 
				    double s2, double x0 = 0);
  static double sum_logL_multi_pe(const std::vector<HistValues>& x, 
				  double F, double a, double b2,
				  double g, double s2, double x0 = 0, 
				  const double tol = 1e-200);
  static double sum_logL_multi_pe_derivs(PMTModelPGDerivs &derivs,
					 const std::vector<HistValues>& x, 
					 double F, double a, double b2,
					 double g, double s2, double x0 = 0,
					 const double tol = 1e-200);
  static double sum_logL_multi_pe_rg(const std::vector<HistValues>& x, 
				     double F, double r, double b2, 
				     double g, double s2, double x0 = 0, 
				     const double tol = 1e-200);
  static double sum_logL_multi_pe_rg_derivs(PMTModelPGDerivs &derivs,
					    const std::vector<HistValues>& x, 
					    double F, double r, double b2,
					    double g, double s2, double x0 = 0,
					    const double tol = 1e-200);


  static double spe_sum_logL_multi_pe(const std::vector<double>& x, 
				      double F, double a, double b2,
				      double g, double s2, double x0 = 0, 
				      const unsigned nmax = 100);

};

} } } // namespace calin::calib::pmt_model_pg
