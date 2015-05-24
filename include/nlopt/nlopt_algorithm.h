/* Copyright (c) 2007-2014 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

#ifndef NLOPT_ALGORITHM_H
#define NLOPT_ALGORITHM_H

typedef enum {
     /* Naming conventions:

        NLOPT_{G/L}{D/N}_* 
	    = global/local derivative/no-derivative optimization, 
              respectively 
 
	*_RAND algorithms involve some randomization.

	*_NOSCAL algorithms are *not* scaled to a unit hypercube
	         (i.e. they are sensitive to the units of x)
	*/

     NLOPT_GN_DIRECT = 0,
     NLOPT_GN_DIRECT_L,
     NLOPT_GN_DIRECT_L_RAND,
     NLOPT_GN_DIRECT_NOSCAL,
     NLOPT_GN_DIRECT_L_NOSCAL,
     NLOPT_GN_DIRECT_L_RAND_NOSCAL,

     NLOPT_GN_ORIG_DIRECT,
     NLOPT_GN_ORIG_DIRECT_L,

     NLOPT_GD_STOGO,
     NLOPT_GD_STOGO_RAND,

     NLOPT_LD_LBFGS_NOCEDAL,

     NLOPT_LD_LBFGS,

     NLOPT_LN_PRAXIS,

     NLOPT_LD_VAR1,
     NLOPT_LD_VAR2,

     NLOPT_LD_TNEWTON,
     NLOPT_LD_TNEWTON_RESTART,
     NLOPT_LD_TNEWTON_PRECOND,
     NLOPT_LD_TNEWTON_PRECOND_RESTART,

     NLOPT_GN_CRS2_LM,

     NLOPT_GN_MLSL,
     NLOPT_GD_MLSL,
     NLOPT_GN_MLSL_LDS,
     NLOPT_GD_MLSL_LDS,

     NLOPT_LD_MMA,

     NLOPT_LN_COBYLA,

     NLOPT_LN_NEWUOA,
     NLOPT_LN_NEWUOA_BOUND,

     NLOPT_LN_NELDERMEAD,
     NLOPT_LN_SBPLX,

     NLOPT_LN_AUGLAG,
     NLOPT_LD_AUGLAG,
     NLOPT_LN_AUGLAG_EQ,
     NLOPT_LD_AUGLAG_EQ,

     NLOPT_LN_BOBYQA,

     NLOPT_GN_ISRES,

     /* new variants that require local_optimizer to be set,
	not with older constants for backwards compatibility */
     NLOPT_AUGLAG,
     NLOPT_AUGLAG_EQ,
     NLOPT_G_MLSL,
     NLOPT_G_MLSL_LDS,

     NLOPT_LD_SLSQP,

     NLOPT_LD_CCSAQ,

     NLOPT_GN_ESCH,

     NLOPT_NUM_ALGORITHMS /* not an algorithm, just the number of them */
} nlopt_algorithm;

#endif // NLOPT_ALGORITHM_H
