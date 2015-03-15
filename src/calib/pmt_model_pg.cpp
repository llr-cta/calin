#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>

#include "calib/pmt_model_pg.hpp"

using namespace calin::calib::pmt_model_pg;

//#define DEBUG_LASERCALIB  

namespace {

static inline double SQR(const double x) 
{ 
  return x*x; 
}

static const double scGAUSSNORM = 0.5*std::log(2.0*M_PI);

static inline double l_poisson(double n, double Fa)
{ 
  return n*std::log(Fa)-Fa-lgamma(n+1.0); 
}

static inline double l_poisson_summand(double n, double Fa,
				       double x, double g, 
				       double var_n, double var_s)
{
  const double var = n*var_n + var_s;
  return l_poisson(n,Fa) - 0.5*(SQR(x-n*g)/var + std::log(var));
}

static unsigned find_n_max(double Fa, double x, double g, 
			 double var_n, double var_s)
{
  // Golden search for "n" max between bounds of F*a and x/g
  unsigned na = std::min(std::max(0,abs(int(x/g))),int(Fa));
  unsigned nb = std::max(std::max(0,abs(int(x/g))),int(Fa))+20;
  unsigned nc = na;
  //double xa = l_poisson_summand(na, Fa, x, g, var_n, var_s);
  //double xb = l_poisson_summand(nb, Fa, x, g, var_n, var_s);
  double xc = 0; /*xa; */
  unsigned iloop = 0;
  while(nb-na>3 && iloop<200)
    {
      if(nc==na || nc==nb)
	{
	  nc = na + floor(double(nb-na)/2.6180339887498948482);
	  xc = l_poisson_summand(nc, Fa, x, g, var_n, var_s);
	}
      unsigned nd = na + nb - nc;
      double xd = l_poisson_summand(nd, Fa, x, g, var_n, var_s);
#if 0
      std::cout << iloop << ' '
		<< na << ' ' << nb << ' ' << nc << ' ' << nd << ' '
	/*<< xa << ' ' << xb << ' ' */<< xc << ' ' << xd << '\n' ;
#endif
      if(xc>xd) { 
	if(nd>nc) { /*xb=xd;*/ nb=nd; }
	else { /*xa=xd;*/ na=nd; }
      }
      else { 
	if(nd>nc) { /*xa=xc;*/ na=nc; xc=xd; nc=nd; }
	else { /*xb=xc;*/ nb=nc; xc=xd; nc=nd; }
      }
      iloop++;
    }
  unsigned nmax = na;
  double xmax =  l_poisson_summand(nmax, Fa, x, g, var_n, var_s);
  for(nc=na+1;nc<nb;nc++)
    {
      xc = l_poisson_summand(nc, Fa, x, g, var_n, var_s);
      if(xc > xmax)xmax=xc, nmax=nc;
    }
  return nmax;
}

static inline double 
l_poisson_deriv_summand(const double n, 
			const double F, const double a, const double Fa,
			const double x, const double x2,
			const double g, const double g2,
			const double b2, const double b2g2, const double s2,
			double& dfdF, double& dfdar, double& dfdb2, 
			double& dfdg, double& dfds2, double& dfdx0, 
			double& logIoffset)
{
  const double var   = n*b2g2 + s2;
  const double var2  = var*var;
  const double n2    = n*n;
  const double xn    = x-n*g;
  double f = 
    l_poisson(n,Fa) - 0.5*(SQR(xn)/var + std::log(var));
  if(logIoffset == 0)logIoffset = f, f = 0;
  else f -= logIoffset;
  f = std::exp(f);
  if(std::isnan(f) || std::isinf(f))return f;
  dfdF  += f*(n-Fa)/F;
  dfdar += f*(n-Fa)/a;
  dfdb2 += f*g2*n*(g2*n2 - b2*g2*n - 2*g*n*x + x2 - s2)/(2*var2);
  dfdg  += 
    -f*n*(n*b2*b2g2*g + n*b2g2*x + b2*g*s2 - b2*g*x2 + n*g*s2 - s2*x)/var2;
  dfds2 += f*(g2*n2 - b2g2*n - 2*g*n*x + x2 - s2)/(2*var2);
  dfdx0 += f*xn/var;
  return f;
}

static void scale_and_clamp(double &x, double I)
{
  x = std::min(std::numeric_limits<double>::max(),
	       std::max(-std::numeric_limits<double>::max(),
			x/I));
}

} // anonymous namespace

// ============================================================================
//
// Log likelihood for individual values
//
// ============================================================================

double PMTModelPG::
logL_ped(double x, double s2, double x0)
{
  const double x2 = SQR(x-x0);
  return -0.5*(x2/s2 + std::log(s2)) - scGAUSSNORM;
}

double PMTModelPG::
logL_ped_derivs(PMTModelPGDerivs &derivs, double x, double s2, double x0)
{
  const double x2 = SQR(x-x0);
  derivs.dfdF  = 0;
  derivs.dfdar = 0;
  derivs.dfdb2 = 0;
  derivs.dfdg  = 0;
  derivs.dfds2 = 0.5*(x2/s2-1.0)/s2;
  derivs.dfdx0 = (x-x0)/s2;
  return -0.5*(x2/s2 + std::log(s2)) - scGAUSSNORM;
}

double PMTModelPG::
logL_multi_pe(double x, double F, double a, double b2, double g, double s2,
	      double x0, const double tol)
{
  const double xp    = x-x0;
  const double g2    = g*g;
  const double var_n = b2*g2;
  const double var_s = s2;
  const double Fa    = F*a;
  unsigned nhi       = find_n_max(Fa,x,g,var_n,var_s);
  unsigned nlo       = nhi;

  double logIoffset = l_poisson_summand(nlo, Fa, xp, g, var_n, var_s);
  if(std::isnan(logIoffset) || std::isinf(logIoffset))
    {
#ifdef DEBUG_LASERCALIB  
      std::cout << "-inf" << '\n';
#endif
      return -std::numeric_limits<double>::infinity();
    }
  double I = 1;

  while(true)
    {
      nhi += 1;
      double sp = 
	l_poisson_summand(nhi, Fa, xp, g, var_n, var_s) - logIoffset;
      //     std::cout << nhi << ' ' << exp(l_poisson(nhi,Fa)) << ' ' << sp << ' ' << I << '\n';
      sp = std::exp(sp);
      if(std::isnan(sp) || std::isinf(sp))break;
      I += sp;
      if(nhi-nlo > 1000){ return 0; 
	std::cout << x << ' ' << a << ' ' << b2 << ' ' << g << ' ' << s2 << ' '
		  << nhi << '\n';}
      if(sp/I < tol)break;
    }
  const double Icutoff = I;
  while(nlo>0)
    {
      nlo -= 1;
      double sp = 
	l_poisson_summand(nlo, Fa, xp, g, var_n, var_s) - logIoffset;
      //     std::cout << nlo << ' ' << exp(l_poisson(nlo,Fa)) << ' ' << sp << ' ' << I << '\n';
      sp = std::exp(sp);
      if(std::isnan(sp) || std::isinf(sp))break;
      I += sp;
      if(sp/Icutoff < tol)break;
    }
  I = std::log(I) + logIoffset - scGAUSSNORM;
  return I;
}

double PMTModelPG::
logL_multi_pe_derivs(PMTModelPGDerivs &derivs, 
	   double x, double F, double a, double b2, double g, double s2,
	   double x0, const double tol)
{
  const double xp    = x-x0;
  const double Fa    = F*a;
  const double g2    = g*g;
  const double x2    = xp*xp;
  const double b2g2  = b2*g2;
  unsigned nhi       = find_n_max(Fa,xp,g,b2g2,s2);
  unsigned nlo       = nhi;

#ifdef DEBUG_LASERCALIB  
  std::cout 
    << x << ' ' << F << ' ' << a << ' ' << b2 << ' ' << g << ' ' << s2 << ' ';
#endif

  derivs.dfdF  = 0;
  derivs.dfdar = 0;
  derivs.dfdb2 = 0;
  derivs.dfdg  = 0;
  derivs.dfds2 = 0;
  derivs.dfdx0 = 0;
  
  double logIoffset = 0;
  double I = 
    l_poisson_deriv_summand(nlo, F, a, Fa, xp, x2, g, g2, b2, b2g2, s2,
			    derivs.dfdF, derivs.dfdar, derivs.dfdb2, 
			    derivs.dfdg, derivs.dfds2, derivs.dfdx0,
			    logIoffset);
  if(std::isnan(logIoffset) || std::isinf(logIoffset))
    {
#ifdef DEBUG_LASERCALIB  
      std::cout << "-inf" << '\n';
#endif
      derivs.dfdF  = 0;
      derivs.dfdar = 0;
      derivs.dfdb2 = 0;
      derivs.dfdg  = 0;
      derivs.dfds2 = 0;
      derivs.dfdx0 = 0;
      return -std::numeric_limits<double>::infinity();
    }
  
  while(true)
    {
      nhi += 1;
      double sp = 
	l_poisson_deriv_summand(nhi, F, a, Fa, xp, x2, g, g2, b2, b2g2, s2,
				derivs.dfdF, derivs.dfdar, derivs.dfdb2, 
				derivs.dfdg, derivs.dfds2, derivs.dfdx0,
				logIoffset);
      //std::cout << nhi << ' ' << exp(l_poisson(nhi,Fa)) << ' ' << sp << ' ' << I << ' ' << derivs.dfdg << '\n';
      if(std::isnan(sp) || std::isinf(sp))break;
      I += sp;
      if(sp/I < tol)break;
    }
  const double I_cutoff = I;
  while(nlo>0)
    {
      nlo -= 1;
      double sp = 
	l_poisson_deriv_summand(nlo, F, a, Fa, xp, x2, g, g2, b2, b2g2, s2,
				derivs.dfdF, derivs.dfdar, derivs.dfdb2, 
				derivs.dfdg, derivs.dfds2, derivs.dfdx0,
				logIoffset);
      //std::cout << nlo << ' ' << exp(l_poisson(nhi,Fa)) << ' ' << sp << ' ' << I << ' ' << derivs.dfdg << '\n';
      if(std::isnan(sp) || std::isinf(sp))break;
      I += sp;
      if(sp/I_cutoff < tol)break;
    }

#ifdef DEBUG_LASERCALIB  
  std::cout << x << ' ' << F << ' ' << I << ' '
	    << derivs.dfdF << ' ' << std::log(I) + logIoffset - scGAUSSNORM
	    << '\n';
#endif
  
  scale_and_clamp(derivs.dfdF,  I);
  scale_and_clamp(derivs.dfdar, I);
  scale_and_clamp(derivs.dfdb2, I);
  scale_and_clamp(derivs.dfdg,  I);
  scale_and_clamp(derivs.dfds2, I);
  scale_and_clamp(derivs.dfdx0, I);
  I = std::log(I) + logIoffset - scGAUSSNORM;

#ifdef DEBUG_LASERCALIB  
  std::cout << I << '\n';
#endif
  
  return I;
}

double PMTModelPG::
logL_multi_pe_rg(double x, double F, double r, double b2, double g, double s2,
		 double x0, const double tol)
{
  return logL_multi_pe(x, F, r/g, b2, g, s2, x0, tol);
}

double PMTModelPG::
logL_multi_pe_rg_derivs(PMTModelPGDerivs &derivs, 
	  double x, double F, double r, double b2, double g, double s2,
	  double x0, const double tol)
{
  double f = logL_multi_pe_derivs(derivs, x, F, r/g, b2, g, s2, x0, tol);
  derivs.dfdg  -= r/(g*g)*derivs.dfdar;
  derivs.dfdar /= g;
  return f;
}

// ============================================================================
//
// Log likelihood for list of values
//
// ============================================================================

double PMTModelPG::
sum_logL_ped(const std::vector<double>& x,
	     double s2, double x0)
{
  double logl = 0;
  for(const auto& ix : x)logl += logL_ped(ix, s2, x0);
  return logl;
}

double PMTModelPG::
sum_logL_ped_derivs(PMTModelPGDerivs &derivs,
		    const std::vector<double>& x, 
		    double s2, double x0)
{
  double logl = 0;
  derivs = PMTModelPGDerivs();
  PMTModelPGDerivs ix_derivs;
  for(const auto& ix : x)
    {
      logl += logL_ped_derivs(ix_derivs, ix, s2, x0);
      derivs += ix_derivs;
    }
  return logl;
}

double PMTModelPG::
sum_logL_multi_pe(const std::vector<double>& x, 
		  double F, double a, double b2,
		  double g, double s2, double x0, 
		  const double tol)
{
  double logl = 0;
  for(const auto& ix : x)
    logl += logL_multi_pe(ix, F, a, b2, g, s2, x0, tol);
  return logl;
}

double PMTModelPG::
sum_logL_multi_pe_derivs(PMTModelPGDerivs &derivs,
			 const std::vector<double>& x,
			 double F, double a, double b2,
			 double g, double s2, double x0,
			 const double tol)
{
  double logl = 0;
  derivs = PMTModelPGDerivs();
  PMTModelPGDerivs ix_derivs;
  for(const auto& ix : x)
  {
    logl += logL_multi_pe_derivs(ix_derivs, ix, F, a, b2, g, s2, x0, tol);
    derivs += ix_derivs;
  }
  return logl;
}

double PMTModelPG::
sum_logL_multi_pe_rg(const std::vector<double>& x, 
		     double F, double r, double b2, 
		     double g, double s2, double x0, 
		     const double tol)
{
  double logl = 0;
  for(const auto& ix : x)
    logl += logL_multi_pe_rg(ix, F, r, b2, g, s2, x0, tol);
  return logl;
}

double PMTModelPG::
sum_logL_multi_pe_rg_derivs(PMTModelPGDerivs &derivs,
			    const std::vector<double>& x, 
			    double F, double r, double b2,
			    double g, double s2, double x0,
			    const double tol)
{
  double logl = 0;
  derivs = PMTModelPGDerivs();
  PMTModelPGDerivs ix_derivs;
  for(const auto& ix : x)
  {
    logl += 
	logL_multi_pe_rg_derivs(ix_derivs, ix, F, r, b2, g, s2, x0, tol);
    derivs += ix_derivs;
  }
  return logl;
}

// ============================================================================
//
// Log likelihood for list of values
//
// ============================================================================

double PMTModelPG::
sum_logL_ped(const std::vector<HistValues>& x,
	     double s2, double x0)
{
  double logl = 0;
  for(const auto& ix : x)logl+= ix.second * logL_ped(ix.first, s2, x0);
  return logl;
}

double PMTModelPG::
sum_logL_ped_derivs(PMTModelPGDerivs &derivs,
		    const std::vector<HistValues>& x,
		    double s2, double x0)
{
  double logl = 0;
  derivs = PMTModelPGDerivs();
  PMTModelPGDerivs ix_derivs;
  for(const auto& ix : x)
  {
    logl+= ix.second * logL_ped_derivs(ix_derivs, ix.first, s2, x0);
    derivs += (ix_derivs *= ix.second);
  }
  return logl;
}

double PMTModelPG::
sum_logL_multi_pe(const std::vector<HistValues>& x, 
		  double F, double a, double b2,
		  double g, double s2, double x0, 
		  const double tol)
{
  double logl = 0;
  for(const auto& ix : x)
    logl += ix.second * logL_multi_pe(ix.first, F, a, b2, g, s2, x0, tol);
  return logl;
}

double PMTModelPG::
sum_logL_multi_pe_derivs(PMTModelPGDerivs &derivs,
			 const std::vector<HistValues>& x,
			 double F, double a, double b2,
			 double g, double s2, double x0,
			 const double tol)
{
  double logl = 0;
  derivs = PMTModelPGDerivs();
  PMTModelPGDerivs ix_derivs;
  for(const auto& ix : x)
  {
    logl += ix.second *
      logL_multi_pe_derivs(ix_derivs, ix.first, F, a, b2, g, s2, x0, tol);
    derivs += (ix_derivs *= ix.second);
  }
  return logl;
}

double PMTModelPG::
sum_logL_multi_pe_rg(const std::vector<HistValues>& x, 
		     double F, double r, double b2, 
		     double g, double s2, double x0, 
		     const double tol)
{
  double logl = 0;
  for(const auto& ix : x)
    logl += ix.second * logL_multi_pe_rg(ix.first, F, r, b2, g, s2, x0, tol);
  return logl;
}

double PMTModelPG::
sum_logL_multi_pe_rg_derivs(PMTModelPGDerivs &derivs,
			    const std::vector<HistValues>& x, 
			    double F, double r, double b2,
			    double g, double s2, double x0,
			    const double tol)
{
  double logl = 0;
  derivs = PMTModelPGDerivs();
  PMTModelPGDerivs ix_derivs;
  for(const auto& ix : x)
  {
    logl += ix.second *
        logL_multi_pe_rg_derivs(ix_derivs, ix.first, F, r, b2, g, s2, x0, tol);
    derivs += (ix_derivs *= ix.second);
  }
  return logl;
}

// ============================================================================
//
// Simplified (faster?) code for low intensity
//
// ============================================================================

double PMTModelPG::
spe_sum_logL_multi_pe(const std::vector<double>& x, 
		      double F, double a, double b2,
		      double g, double s2, double x0, 
		      const unsigned nmax)
{
  const double g2    = g*g;
  const double var_n = b2*g2;
  const double var_s = s2;
  const double Fa    = F*a;
  
  double* C1 = (double*)alloca(nmax*sizeof(double));
  double* C2 = (double*)alloca(nmax*sizeof(double));
  double* C3 = (double*)alloca(nmax*sizeof(double));
  for(unsigned n=0; n<nmax; n++)
    {
      const double var = n*var_n + var_s;  
      C1[n] = n*std::log(Fa)-Fa-lgamma(n+1.0) - 0.5*std::log(var);
      C2[n] = 0.5/var;
      C3[n] = n*g;
    }

  double logl = 0;
  for(unsigned ix=0;ix<x.size();ix++)
    {
      const double xix = x[ix]-x0;
      double lix = 0;
      for(unsigned n=0; n<nmax; n++)
	{
	  const double sqr_arg = xix-C3[n];
	  const double log_pn = C1[n] - C2[n]*(sqr_arg*sqr_arg);
	  lix += std::exp(log_pn);
	}
      logl += std::log(lix);
    }
  return logl - x.size()*scGAUSSNORM;
}


