/* 

   calin/math/vs_vec4d.cpp -- Stephen Fegan -- 2015-11-05

   Class for 4-vector operations. This code is derived from simulation
   code largely implemented by Vladimir Vassiliev at UCLA in 2004. It
   is not complient with the calin coding style.

   Some portions are :

   Copyright 2004, 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

/*! \file vec4D.cpp
  vec4D class implementation file
  
  \author   Stephen Fegan             \n
            UCLA                      \n
	    sfegan@astro.ucla.edu     \n

  \author   Vladimir Vassiliev        \n
            UCLA                      \n
	    vvv@astro.ucla.edu        \n

  \date     11/09/2004
  \version  1.1
  \note
*/

#include<iostream>

#include<math/vs_vec4d.hpp>

using namespace calin::math::vs_physics;

/// Print vector components and vector norm in square
void Vec4D::Dump(std::ostream& stream) const
{
  stream << std::endl;
  stream << " .T:   " << r0 << std::endl;
  stream << " .X:   " << r.x << std::endl;
  stream << " .Y:   " << r.y << std::endl;
  stream << " .Z:   " << r.z << std::endl;
  stream << "Norm2: " << Norm2() << std::endl;
}

void Vec4D::DumpShort(std::ostream& stream) const
{
  stream << "( " << r0 << ", " << r.x << ", " << r.y << ", " << r.z << " )";
}


#ifdef TESTMAIN

// Compile with: g++ -O3 -DTESTMAIN -o test Vec4D.cpp Vec3D.o -lm

#include<sstream>

int main()
{
  std::cout << "sizeof(Vec4D): " << sizeof(Vec4D) << std::endl;
  std::cout << "sizeof(Vec4D*): " << sizeof(Vec4D*) << std::endl << std::endl;

  Vec4D a(1,2,3,4);
  Vec4D b(10,20,30,40);

  std::cout << "a:" << a << " b:" << b << std::endl;
  std::cout << "a.Norm2():" << a.Norm2() << "  "
	    << "b.Norm2():" << b.Norm2() << std::endl << std::endl;
  
  std::cout << "a+=b: " << (a+=b) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a-=b: " << (a-=b) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a*=1.5: " << (a*=1.5) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a/=1.5: " << (a/=1.5) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a.Reset(): " << a.Reset() << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl;
  std::cout << "a.Reset(b): " << a.Reset(b) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  Vec3D r;

  std::cout << "a.Reset(1,1,0,0): " << a.Reset(1,1,0,0) << std::endl;
  std::cout << "r.Reset(0,0,pi/4): " << r.Reset(0,0,M_PI/2) << std::endl;
  std::cout << "a:" << a << " r:" << r << std::endl;
  std::cout << "a.Rotate(r)" << std::endl;
  a.Rotate(r);
  std::cout << "a:" << a << " r:" << r << std::endl << std::endl;

  std::cout << "a.Reset(1,0,1,0): " << a.Reset(1,0,1,0) << std::endl;
  std::cout << "r.Reset(0,0,pi/4): " << r.Reset(0,0,M_PI/2) << std::endl;
  std::cout << "a:" << a << " r:" << r << std::endl;
  std::cout << "a.Rotate(r)" << std::endl;
  a.Rotate(r);
  std::cout << "a:" << a << " r:" << r << std::endl << std::endl;

  std::cout << "a.Reset(1,0,0,1): " << a.Reset(1,0,0,1) << std::endl;
  std::cout << "r.Reset(0,0,pi/4): " << r.Reset(0,0,M_PI/2) << std::endl;
  std::cout << "a:" << a << " r:" << r << std::endl;
  std::cout << "a.Rotate(r)" << std::endl;
  a.Rotate(r);
  std::cout << "a:" << a << " r:" << r << std::endl << std::endl;

  std::cout << "a.Reset(1,-2,3,-4): " << a.Reset(1,-2,3,-4) << std::endl;
  std::cout << "a.P()" << std::endl;
  a.P();
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  Vec4D c;

  std::cout << "c=b=a=():" << (c=b=a=Vec4D()) << std::endl;
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl << std::endl;
  
  std::cout << "c=b=(a=())+(1,2,0,0)):" 
	    << (c=b=(a=Vec4D())+Vec4D(1,2,0,0)) << std::endl;
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl << std::endl;

  std::cout << "a=(1,1,0,0):" << (a=Vec4D(1,1,0,0)) << std::endl;
  std::cout << "b=(1,0,-1,0):" << (b=Vec4D(1,0,-1,0)) << std::endl;
  std::cout << "c=a+b:" << (c=a+b) << std::endl;
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl << std::endl;
  
  std::cout << "c=a-b:" << (c=a-b) << std::endl;
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl << std::endl;
 
  std::cout << "a*b:" << (a*b) << std::endl;
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl << std::endl;

  std::cout << "(3,4,5,6)*(10,0,0,0): " << (Vec4D(3,4,5,6)*Vec4D(10,0,0,0)) << std::endl;
  std::cout << "(3,4,5,6)*(0,20,0,0): " << (Vec4D(3,4,5,6)*Vec4D(10,0,0,0)) << std::endl;
  std::cout << "(3,4,5,6)*(0,0,30,0): " << (Vec4D(3,4,5,6)*Vec4D(10,0,0,0)) << std::endl;
  std::cout << "(3,4,5,6)*(0,0,0,40): " << (Vec4D(3,4,5,6)*Vec4D(10,0,0,0)) << std::endl;
  std::cout << "(3,4,5,6)*(10,20,0,0): " << (Vec4D(3,4,5,6)*Vec4D(10,20,0,0)) << std::endl;
  std::cout << "(3,4,5,6)*(10,0,30,0): " << (Vec4D(3,4,5,6)*Vec4D(10,0,30,0)) << std::endl;
  std::cout << "(3,4,5,6)*(10,0,0,40): " << (Vec4D(3,4,5,6)*Vec4D(10,0,0,40)) << std::endl;
  std::cout << "(3,4,5,6)*(0,20,30,0): " << (Vec4D(3,4,5,6)*Vec4D(0,20,30,0)) << std::endl;
  std::cout << "(3,4,5,6)*(0,20,0,40): " << (Vec4D(3,4,5,6)*Vec4D(0,20,0,40)) << std::endl;
  std::cout << "(3,4,5,6)*(0,0,30,40): " << (Vec4D(3,4,5,6)*Vec4D(0,0,30,40)) << std::endl;
  std::cout << "(3,4,5,6)*(10,20,30,0): " << (Vec4D(3,4,5,6)*Vec4D(10,20,30,0)) << std::endl;
  std::cout << "(3,4,5,6)*(10,20,0,40): " << (Vec4D(3,4,5,6)*Vec4D(10,20,0,40)) << std::endl;
  std::cout << "(3,4,5,6)*(10,0,30,40): " << (Vec4D(3,4,5,6)*Vec4D(10,0,30,40)) << std::endl;
  std::cout << "(3,4,5,6)*(0,20,30,40): " << (Vec4D(3,4,5,6)*Vec4D(0,20,30,40)) << std::endl;
  std::cout << "(3,4,5,6)*(10,20,30,40): " << (Vec4D(3,4,5,6)*Vec4D(10,20,30,40)) << std::endl;
  std::cout << std::endl;

  std::cout << "(3,4,5,6)+(10,0,0,0): " << (Vec4D(3,4,5,6)+Vec4D(10,0,0,0)) << std::endl;
  std::cout << "(3,4,5,6)+(0,20,0,0): " << (Vec4D(3,4,5,6)+Vec4D(10,0,0,0)) << std::endl;
  std::cout << "(3,4,5,6)+(0,0,30,0): " << (Vec4D(3,4,5,6)+Vec4D(10,0,0,0)) << std::endl;
  std::cout << "(3,4,5,6)+(0,0,0,40): " << (Vec4D(3,4,5,6)+Vec4D(10,0,0,0)) << std::endl;
  std::cout << "(3,4,5,6)+(10,20,0,0): " << (Vec4D(3,4,5,6)+Vec4D(10,20,0,0)) << std::endl;
  std::cout << "(3,4,5,6)+(10,0,30,0): " << (Vec4D(3,4,5,6)+Vec4D(10,0,30,0)) << std::endl;
  std::cout << "(3,4,5,6)+(10,0,0,40): " << (Vec4D(3,4,5,6)+Vec4D(10,0,0,40)) << std::endl;
  std::cout << "(3,4,5,6)+(0,20,30,0): " << (Vec4D(3,4,5,6)+Vec4D(0,20,30,0)) << std::endl;
  std::cout << "(3,4,5,6)+(0,20,0,40): " << (Vec4D(3,4,5,6)+Vec4D(0,20,0,40)) << std::endl;
  std::cout << "(3,4,5,6)+(0,0,30,40): " << (Vec4D(3,4,5,6)+Vec4D(0,0,30,40)) << std::endl;
  std::cout << "(3,4,5,6)+(10,20,30,0): " << (Vec4D(3,4,5,6)+Vec4D(10,20,30,0)) << std::endl;
  std::cout << "(3,4,5,6)+(10,20,0,40): " << (Vec4D(3,4,5,6)+Vec4D(10,20,0,40)) << std::endl;
  std::cout << "(3,4,5,6)+(10,0,30,40): " << (Vec4D(3,4,5,6)+Vec4D(10,0,30,40)) << std::endl;
  std::cout << "(3,4,5,6)+(0,20,30,40): " << (Vec4D(3,4,5,6)+Vec4D(0,20,30,40)) << std::endl;
  std::cout << "(3,4,5,6)+(10,20,30,40): " << (Vec4D(3,4,5,6)+Vec4D(10,20,30,40)) << std::endl;
  std::cout << std::endl;

  std::cout << "(3,4,5,6)-(10,0,0,0): " << (Vec4D(3,4,5,6)-Vec4D(10,0,0,0)) << std::endl;
  std::cout << "(3,4,5,6)-(0,20,0,0): " << (Vec4D(3,4,5,6)-Vec4D(10,0,0,0)) << std::endl;
  std::cout << "(3,4,5,6)-(0,0,30,0): " << (Vec4D(3,4,5,6)-Vec4D(10,0,0,0)) << std::endl;
  std::cout << "(3,4,5,6)-(0,0,0,40): " << (Vec4D(3,4,5,6)-Vec4D(10,0,0,0)) << std::endl;
  std::cout << "(3,4,5,6)-(10,20,0,0): " << (Vec4D(3,4,5,6)-Vec4D(10,20,0,0)) << std::endl;
  std::cout << "(3,4,5,6)-(10,0,30,0): " << (Vec4D(3,4,5,6)-Vec4D(10,0,30,0)) << std::endl;
  std::cout << "(3,4,5,6)-(10,0,0,40): " << (Vec4D(3,4,5,6)-Vec4D(10,0,0,40)) << std::endl;
  std::cout << "(3,4,5,6)-(0,20,30,0): " << (Vec4D(3,4,5,6)-Vec4D(0,20,30,0)) << std::endl;
  std::cout << "(3,4,5,6)-(0,20,0,40): " << (Vec4D(3,4,5,6)-Vec4D(0,20,0,40)) << std::endl;
  std::cout << "(3,4,5,6)-(0,0,30,40): " << (Vec4D(3,4,5,6)-Vec4D(0,0,30,40)) << std::endl;
  std::cout << "(3,4,5,6)-(10,20,30,0): " << (Vec4D(3,4,5,6)-Vec4D(10,20,30,0)) << std::endl;
  std::cout << "(3,4,5,6)-(10,20,0,40): " << (Vec4D(3,4,5,6)-Vec4D(10,20,0,40)) << std::endl;
  std::cout << "(3,4,5,6)-(10,0,30,40): " << (Vec4D(3,4,5,6)-Vec4D(10,0,30,40)) << std::endl;
  std::cout << "(3,4,5,6)-(0,20,30,40): " << (Vec4D(3,4,5,6)-Vec4D(0,20,30,40)) << std::endl;
  std::cout << "(3,4,5,6)-(10,20,30,40): " << (Vec4D(3,4,5,6)-Vec4D(10,20,30,40)) << std::endl;
  std::cout << std::endl;

  std::cout << "-(3,-4,5,-6): " << (-Vec4D(3,-4,5,-6)) << std::endl;
  std::cout << "1.5*(3,4,5,6): " << (1.5*Vec4D(3,4,5,6)) << std::endl;
  std::cout << "(3,4,5,6)*1.5: " << (Vec4D(3,4,5,6)*1.5) << std::endl;
  std::cout << "(3,4,5,6)/1.5: " << (Vec4D(3,4,5,6)/1.5) << std::endl << std::endl;

  std::cout << "std::istringstream s(\"(10,11,12,13)\"); s >> a;" << std::endl;
  std::istringstream s("(10,11,12,13)");
  s >> a;
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl << std::endl;
}

#endif // TESTMAIN
