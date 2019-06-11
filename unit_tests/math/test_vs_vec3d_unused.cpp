/*

   calin/unit_tests/math/test_vs_vec3d.cpp -- Stephen Fegan -- 2015-11-15

   Unit tests for Vec3D

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>
#include <vector>
#include <limits>

#include "math/vs_vec3d.hpp"

using namespace calin::math::vs_physics;

double sqrt_eps = std::sqrt(std::numeric_limits<double>::epsilon());

TEST(TestVec3D, BasicMathematicalOperations) {
  Vec3D a(10,20,30);
  Vec3D b(1,2,3);

  EXPECT_EQ(a.x(), 10);
  EXPECT_EQ(a.y(), 20);
  EXPECT_EQ(a.z(), 30);
  EXPECT_EQ(b.x(), 1);
  EXPECT_EQ(b.y(), 2);
  EXPECT_EQ(b.z(), 3);

  EXPECT_EQ(a, Vec3D(10,20,30));
  EXPECT_NE(a, Vec3D(11,20,30));
  EXPECT_NE(a, Vec3D(10,21,30));
  EXPECT_NE(a, Vec3D(10,20,31));
  EXPECT_EQ(b, Vec3D(1,2,3));

  EXPECT_EQ(a.Norm2(), 1400.0);
  EXPECT_EQ(b.Norm2(), 14.0);

  EXPECT_EQ(a.Norm(), std::sqrt(1400.0));
  EXPECT_EQ(b.Norm(), std::sqrt(14.0));

  EXPECT_EQ(a+=b, Vec3D(11,22,33));
  EXPECT_EQ(a, Vec3D(11,22,33));
  EXPECT_EQ(b, Vec3D(1,2,3));

  EXPECT_EQ(a-=b, Vec3D(10,20,30));
  EXPECT_EQ(a, Vec3D(10,20,30));
  EXPECT_EQ(b, Vec3D(1,2,3));

  EXPECT_EQ(a*=1.5, Vec3D(15,30,45));
  EXPECT_EQ(a, Vec3D(15,30,45));

  EXPECT_EQ(a/=1.5, Vec3D(10,20,30));
  EXPECT_EQ(a, Vec3D(10,20,30));

  EXPECT_EQ(a.Reset(), Vec3D(0,0,0));
  EXPECT_EQ(a, Vec3D(0,0,0));
  EXPECT_EQ(b.Reset(), Vec3D(0,0,0));
  EXPECT_EQ(b, Vec3D(0,0,0));

  EXPECT_EQ(a.Reset(1,0,0), Vec3D(1,0,0));
  EXPECT_EQ(a, Vec3D(1,0,0));
  EXPECT_EQ(a.Reset(0,1,0), Vec3D(0,1,0));
  EXPECT_EQ(a, Vec3D(0,1,0));
  EXPECT_EQ(a.Reset(0,0,1), Vec3D(0,0,1));
  EXPECT_EQ(a, Vec3D(0,0,1));

#if 0
  EXPECT_EQ(a.Reset(1,0,0), Vec3D(1,0,0));
  EXPECT_EQ(b.Reset(0,0,M_PI/2), Vec3D(0,0,M_PI/2));
  a.Rotate(b);
  EXPECT_NEAR(a.x(), 0.0, sqrt_eps);
  EXPECT_NEAR(a.y(), 1.0, sqrt_eps);
  EXPECT_NEAR(a.z(), 0.0, sqrt_eps);
  EXPECT_EQ(b, Vec3D(0,0,M_PI/2));

  EXPECT_EQ(a.Reset(0,1,0), Vec3D(0,1,0));
  EXPECT_EQ(b.Reset(0,0,M_PI/2), Vec3D(0,0,M_PI/2));
  a.Rotate(b);
  EXPECT_NEAR(a.x(), -1.0, sqrt_eps);
  EXPECT_NEAR(a.y(), 0.0, sqrt_eps);
  EXPECT_NEAR(a.z(), 0.0, sqrt_eps);
  EXPECT_EQ(b, Vec3D(0,0,M_PI/2));

  EXPECT_EQ(a.Reset(0,0,1), Vec3D(0,0,1));
  EXPECT_EQ(b.Reset(0,0,M_PI/2), Vec3D(0,0,M_PI/2));
  a.Rotate(b);
  EXPECT_NEAR(a.x(), 0.0, sqrt_eps);
  EXPECT_NEAR(a.y(), 0.0, sqrt_eps);
  EXPECT_NEAR(a.z(), 1.0, sqrt_eps);
  EXPECT_EQ(b, Vec3D(0,0,M_PI/2));
#endif

  EXPECT_EQ(a.Reset(1,-2,3), Vec3D(1,-2,3));
  a.P();
  EXPECT_EQ(a, Vec3D(-1,2,-3));

  EXPECT_EQ(a.Reset(1,0,0), Vec3D(1,0,0));
  EXPECT_EQ(b.Reset(1,0,0), Vec3D(1,0,0));
  EXPECT_EQ(a^=b, Vec3D(0,0,0));
  EXPECT_EQ(a, Vec3D(0,0,0));
  EXPECT_EQ(b, Vec3D(1,0,0));

  EXPECT_EQ(a.Reset(1,0,0), Vec3D(1,0,0));
  EXPECT_EQ(b.Reset(0,1,0), Vec3D(0,1,0));
  EXPECT_EQ(a^=b, Vec3D(0,0,1));
  EXPECT_EQ(a, Vec3D(0,0,1));
  EXPECT_EQ(b, Vec3D(0,1,0));

  EXPECT_EQ(a.Reset(1,0,0), Vec3D(1,0,0));
  EXPECT_EQ(b.Reset(0,0,1), Vec3D(0,0,1));
  EXPECT_EQ(a^=b, Vec3D(0,-1,0));
  EXPECT_EQ(a, Vec3D(0,-1,0));
  EXPECT_EQ(b, Vec3D(0,0,1));

  Vec3D c;

  EXPECT_EQ(a=b=c=Vec3D(), Vec3D());
  EXPECT_EQ(a, Vec3D());
  EXPECT_EQ(b, Vec3D());
  EXPECT_EQ(c, Vec3D());

  EXPECT_EQ(c=b=(a=Vec3D())+Vec3D(1,0,0), Vec3D(1,0,0));
  EXPECT_EQ(a, Vec3D());
  EXPECT_EQ(b, Vec3D(1,0,0));
  EXPECT_EQ(c, Vec3D(1,0,0));

  a=Vec3D(1,0,0);
  b=Vec3D(1,-1,0);
  EXPECT_EQ(c=a+b, Vec3D(2,-1,0));
  EXPECT_EQ(a, Vec3D(1,0,0));
  EXPECT_EQ(b, Vec3D(1,-1,0));
  EXPECT_EQ(c, Vec3D(2,-1,0));

  EXPECT_EQ(c=a-b, Vec3D(0,1,0));
  EXPECT_EQ(a, Vec3D(1,0,0));
  EXPECT_EQ(b, Vec3D(1,-1,0));
  EXPECT_EQ(c, Vec3D(0,1,0));

  EXPECT_EQ(a*b, 1);
  EXPECT_EQ(a, Vec3D(1,0,0));
  EXPECT_EQ(b, Vec3D(1,-1,0));

  EXPECT_EQ(Vec3D(3,4,5)*Vec3D(2,0,0), 6);
  EXPECT_EQ(Vec3D(3,4,5)*Vec3D(0,3,0), 12);
  EXPECT_EQ(Vec3D(3,4,5)*Vec3D(0,0,4), 20);
  EXPECT_EQ(Vec3D(3,4,5)*Vec3D(2,3,0), 18);
  EXPECT_EQ(Vec3D(3,4,5)*Vec3D(2,0,4), 26);
  EXPECT_EQ(Vec3D(3,4,5)*Vec3D(0,3,4), 32);
  EXPECT_EQ(Vec3D(3,4,5)*Vec3D(2,3,4), 38);

  EXPECT_EQ(Vec3D(3,4,5)^Vec3D(2,0,0), Vec3D(0, 10, -8));
  EXPECT_EQ(Vec3D(3,4,5)^Vec3D(0,3,0), Vec3D(-15, 0, 9));
  EXPECT_EQ(Vec3D(3,4,5)^Vec3D(0,0,4), Vec3D(16, -12, 0));
  EXPECT_EQ(Vec3D(3,4,5)^Vec3D(2,3,0), Vec3D(-15, 10, 1));
  EXPECT_EQ(Vec3D(3,4,5)^Vec3D(2,0,4), Vec3D(16, -2, -8));
  EXPECT_EQ(Vec3D(3,4,5)^Vec3D(0,3,4), Vec3D(1, -12, 9));
  EXPECT_EQ(Vec3D(3,4,5)^Vec3D(2,3,4), Vec3D(1, -2, 1));

  EXPECT_EQ(Vec3D(3,4,5)+Vec3D(2,0,0), Vec3D(5, 4, 5));
  EXPECT_EQ(Vec3D(3,4,5)+Vec3D(0,3,0), Vec3D(3, 7, 5));
  EXPECT_EQ(Vec3D(3,4,5)+Vec3D(0,0,4), Vec3D(3, 4, 9));
  EXPECT_EQ(Vec3D(3,4,5)+Vec3D(2,3,0), Vec3D(5, 7, 5));
  EXPECT_EQ(Vec3D(3,4,5)+Vec3D(2,0,4), Vec3D(5, 4, 9));
  EXPECT_EQ(Vec3D(3,4,5)+Vec3D(0,3,4), Vec3D(3, 7, 9));
  EXPECT_EQ(Vec3D(3,4,5)+Vec3D(2,3,4), Vec3D(5, 7, 9));

  EXPECT_EQ(Vec3D(3,4,5)-Vec3D(2,0,0), Vec3D(1, 4, 5));
  EXPECT_EQ(Vec3D(3,4,5)-Vec3D(0,3,0), Vec3D(3, 1, 5));
  EXPECT_EQ(Vec3D(3,4,5)-Vec3D(0,0,4), Vec3D(3, 4, 1));
  EXPECT_EQ(Vec3D(3,4,5)-Vec3D(2,3,0), Vec3D(1, 1, 5));
  EXPECT_EQ(Vec3D(3,4,5)-Vec3D(2,0,4), Vec3D(1, 4, 1));
  EXPECT_EQ(Vec3D(3,4,5)-Vec3D(0,3,4), Vec3D(3, 1, 1));
  EXPECT_EQ(Vec3D(3,4,5)-Vec3D(2,3,4), Vec3D(1, 1, 1));

  EXPECT_EQ(-Vec3D(3,4,5), Vec3D(-3,-4,-5));
  EXPECT_EQ(1.5*Vec3D(3,4,5), Vec3D(4.5,6,7.5));
  EXPECT_EQ(Vec3D(3,4,5)*1.5, Vec3D(4.5,6,7.5));
  EXPECT_EQ(Vec3D(3,4,5)/1.5, Vec3D(2,4/1.5,5/1.5));

#if 0
  std::istringstream s("(10,11,12)");
  s >> a;
  EXPECT_EQ(a, Vec3D(10,11,12));
#endif

#if 0
  for(unsigned nrot = 1; nrot<=20; nrot++)
    {
      double mean_residual = 0;
      double max_residual = 0;

      for(unsigned n=0; n<1000; n++)
	{
	  double phi = rng.Uniform() * M_PI * 2.0;
	  double costheta = rng.Uniform() * 2.0 - 1.0;
	  double sintheta = sin(acos(costheta));
	  Vec3D a(sintheta*cos(phi),sintheta*sin(phi),costheta);
	  Vec3D b(a);
	  Vec3D rc;

	  std::vector<Vec3D> R;
	  std::vector<Vec3D> Rc;

	  for(unsigned i=0; i<nrot; i++)
	    {
	      phi = rng.Uniform() * M_PI * 2.0;
	      costheta = rng.Uniform() * 2.0 - 1.0;
	      sintheta = sin(acos(costheta));
	      Vec3D r(sintheta*cos(phi),sintheta*sin(phi),costheta);
	      r *= rng.Uniform()*2.0*M_PI;

	      b.Rotate(r);
	      rc &= r;

	      R.push_back(r);
	      Rc.push_back(rc);
	    }

	  a.Rotate(rc);

	  Vec3D c(a-b);
	  double residual = c.Norm();
	  if(residual>max_residual)max_residual=residual;
	  mean_residual+=residual;

#if 0
	  if(residual>1e-14)
	    {
	      for(unsigned i=0; i<R.size(); i++)
		std::cout << R[i].Norm()
			  << ' '  << R[i] << ' ' << Rc[i] << std::endl;
	      std::cout << a << ' ' << b << ' ' << c << std::endl;
	    }
#endif
	}

      std::cout << nrot << " rotations: mean residual: "
		<< mean_residual/1000.0
		<< " max residual: " << max_residual << std::endl;
    }
#endif
}

#if 0
TEST(TestVec3D, RotationComposition) {
  Vec3D a;
  Vec3D b;
  Vec3D c;

  Vec3D r1=M_PI/2*Vec3D(1,0,0);
  Vec3D r2=M_PI/2*Vec3D(0,1,0);
  Vec3D rc = r1 & r2;

  a=Vec3D(1,0,0); b=a; b.Rotate(r1); b.Rotate(r2); c=a; c.Rotate(rc);
  EXPECT_NEAR(b.x(), 0, sqrt_eps);
  EXPECT_NEAR(b.y(), 0, sqrt_eps);
  EXPECT_NEAR(b.z(), -1, sqrt_eps);
  EXPECT_NEAR(b.x(), c.x(), sqrt_eps);
  EXPECT_NEAR(b.y(), c.y(), sqrt_eps);
  EXPECT_NEAR(b.z(), c.z(), sqrt_eps);

  a=Vec3D(0,1,0); b=a; b.Rotate(r1); b.Rotate(r2); c=a; c.Rotate(rc);
  EXPECT_NEAR(b.x(), 1, sqrt_eps);
  EXPECT_NEAR(b.y(), 0, sqrt_eps);
  EXPECT_NEAR(b.z(), 0, sqrt_eps);
  EXPECT_NEAR(b.x(), c.x(), sqrt_eps);
  EXPECT_NEAR(b.y(), c.y(), sqrt_eps);
  EXPECT_NEAR(b.z(), c.z(), sqrt_eps);

  a=Vec3D(0,0,2); b=a; b.Rotate(r1); b.Rotate(r2); c=a; c.Rotate(rc);
  EXPECT_NEAR(b.x(), 0, sqrt_eps);
  EXPECT_NEAR(b.y(), -2, sqrt_eps);
  EXPECT_NEAR(b.z(), 0, sqrt_eps);
  EXPECT_NEAR(b.x(), c.x(), sqrt_eps);
  EXPECT_NEAR(b.y(), c.y(), sqrt_eps);
  EXPECT_NEAR(b.z(), c.z(), sqrt_eps);

  a=Vec3D(1,1,0); b=a; b.Rotate(r1); b.Rotate(r2); c=a; c.Rotate(rc);
  EXPECT_NEAR(b.x(), 1, sqrt_eps);
  EXPECT_NEAR(b.y(), 0, sqrt_eps);
  EXPECT_NEAR(b.z(), -1, sqrt_eps);
  EXPECT_NEAR(b.x(), c.x(), sqrt_eps);
  EXPECT_NEAR(b.y(), c.y(), sqrt_eps);
  EXPECT_NEAR(b.z(), c.z(), sqrt_eps);

  a=Vec3D(.1,0,0);
  a &= a & a & a;
  EXPECT_NEAR(a.x(), 0.4, sqrt_eps);
  EXPECT_NEAR(a.y(), 0, sqrt_eps);
  EXPECT_NEAR(a.z(), 0, sqrt_eps);
}
#endif

TEST(TestVec3D, ProtoSaveAndRestore) {
  Vec3D x(1.0,2.0,3.0);
  calin::ix::common_types::Vector3D x_data;
  x.dump_as_proto(&x_data);
  Vec3D y(x_data);
  EXPECT_EQ(x,y);
}
