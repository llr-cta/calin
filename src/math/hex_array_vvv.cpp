/* 

   calin/io/hex_array_vvv.cpp -- Stephen Fegan -- 2015-10-21

   Vladimir Vassiliev's original implementation of the hex grid
   functions adapted into calin namespace primarily for testing of new
   implementation.

*/

#include "math/hex_array.hpp"

/* xytohex.c */
/*
 * This set of routines provide conversion between 
 * regular coordinate system (x,y) and coordinate
 * system of hexagonal cells covering all 2D space. 
 * The spacing between cells in x direction is equal 
 * to 1 which is also the side to side cell size.
 * The cell in origin of (x,y) plane is cell number
 * one. The center of cell number two is at (1,0).
 * Then cells of the first hexagonal ring are numbered 
 * in clockwise direction around (x,y) origin. The center 
 * of cell number eight is at (2,0) and cells of the 
 * second hexagonal ring are numbered in clockwise 
 * direction again. This numbering process is continued
 * until all 2D space is covered. There are two routines
 * available for a user. Routine  
 *
 * xy_to_nh(double *x, double *y, int *n);
 *
 * maps (x,y) coordinates onto cell number n and (x,y)
 * pair of coordinates in translated reference frame 
 * which has origin at the center of the cell number n.
 * Routine 
 *
 * nh_to_xy(int *n, double *x, double *y);
 *
 * performs inverse transformation. Using given cell
 * number, n, it returns (x,y) coordinates of the
 * cell center. Because all cells are numbered by 
 * integers >0, n must be >0. If routine is called 
 * with negative or zero n, n value will be re-assigned 
 * to 1 (origin).
 *
 * March 24, 2000
 *
 * vvv
 */ 

#include <cfloat>

namespace {

/* STRUCTURE DECLARATIONS */
/*  vector structure          */
  typedef struct vec {
     double        x;  /* x  */ 
     double        y;  /* y */
    } vec;

/*  hex coordinate structure  */
  typedef struct hex {
     int           m;  /* radius  */ 
     int           j;  /* segment */
     int           l;  /* cell    */
    } hex;


/* DEFINE STATEMENTS */
#define s 1.73205080756887729352744634150587  /* sqrt(3.) */




double 	hex_nghbr(double x, double y, hex rsc[3])
/*
 * converts (x,y) into "hex" coordinates of three neighbors
 * and returns outer radius of hexagon to which (x,y) belongs.
 *
 * vvv - March 24, 2000
 */
{
  static vec     el[6]={{-1./2.,-s/2.},{  -1.,   0.},{-1/2., s/2.},
                        { 1./2., s/2.},{   1.,   0.},{ 1/2.,-s/2.}};
  static vec     es[6]={{    1.,-1./s},{   0.,-2./s},{  -1.,-1./s},
                        {   -1., 1./s},{   0., 2./s},{   1., 1./s}};
  int    m,j=0,l;
  double r=0.,r_;
  int          i;
  
  /* radius, segment, and cell */
  for(i=0;i<6;i++){
	r_=x*es[i].x+y*es[i].y;
    if(r_>r) {
      r=r_;
      j=i;
    }
  }
  m=(int)floor(r);
  l=(int)floor(x*el[j].x+y*el[j].y+(m+1.)/2.);

  /* 3 neighbors */
  rsc[0].m=  m;
  rsc[0].j=  j;
  rsc[0].l=  l;
  rsc[1].m=m+1;
  rsc[1].j=  j;
  rsc[1].l=  l;
  rsc[2].m=m+1;
  rsc[2].j=  j;
  rsc[2].l=l+1;
  for(i=0;i<3;i++){
	if(rsc[i].l == rsc[i].m) {
	   rsc[i].l=0;
	   rsc[i].j++;
	   if(rsc[i].j == 6) rsc[i].j=0;
	}
  }

  return r;
}

void hex_to_xy(hex hx, double *x, double *y)
/*
 * calculates cell (x,y) coordinate using cell "hex" coordinate
 *
 * vvv - March 24, 2000
 */
{
  static vec     em[6]={{    1.,   0.},{ 1/2.,-s/2.},{-1/2.,-s/2.},
                        {   -1.,   0.},{-1/2., s/2.},{ 1/2., s/2.}};
  static vec     el[6]={{-1./2.,-s/2.},{  -1.,   0.},{-1/2., s/2.},
                        { 1./2., s/2.},{   1.,   0.},{ 1/2.,-s/2.}};

  *x=hx.m*em[hx.j].x+hx.l*el[hx.j].x;
  *y=hx.m*em[hx.j].y+hx.l*el[hx.j].y;

  return;
}

void xy_to_hex(double *x, double *y, hex *rsc)
/*
 * converts (x,y) into "hex" coordinates of the cell,
 * and returns (x,y) in the reference frame in which 
 * cell center is placed in origin.
 *
 * vvv - March 24, 2000
 */
{

  vec           d[3];
  hex   hx[3],hxa[3];

  double     r,xp,yp;
  int         i,k=-1;

  r=hex_nghbr(*x, *y, hx);
  for(i=0;i<3;i++){
    hex_to_xy(hx[i], &xp, &yp);
	d[i].x=*x-xp;
	d[i].y=*y-yp;
    r=s*hex_nghbr(d[i].y,d[i].x, hxa);
    if(r<=1.+DBL_EPSILON) k=i;
  }

  if(k!=-1) {
    rsc->m=hx[k].m;
    rsc->j=hx[k].j;
    rsc->l=hx[k].l;
	*x=d[k].x;
	*y=d[k].y;
  } else {
	fprintf(stdout," Error in xy_to_hex: pixel not found %e %e \n",*x,*y);
    exit(1);
  }

  return;
}

void hex_to_nh(hex hx, int *n)
/*
 * calculates cell number using cell "hex" coordinate
 *
 * vvv - March 24, 2000
 */
{
  *n=1+3*hx.m*(hx.m-1)+hx.m*hx.j+hx.l+1;
  if(hx.m==0) *n=1;	

  return;
}

void nh_to_hex(int n, hex *hx)
/*
 * calculates cell "hex" coordinate using cell number
 *
 * vvv - March 24, 2000
 */
{
  int n_=1;

  hx->m=0;
  hx->j=0;
  hx->l=0;
  if(n==1) return;

  /* radius */
  do{
	hx->m++;
	n_+=6*hx->m;
  }while (n_<n);
  n_-=6*hx->m;

  /* segment */
  do{
	hx->j++;
	n_+=hx->m;
  }while (n_<n);
  n_-=hx->m;
  hx->j--;
  
  /* cell */
  hx->l=n-n_-1;

  return;
}

} // unnamed namespace

void calin::math::hex_array::vvv::
xy_to_nh(double *x, double *y, int *nh)
/*
 * converts (x,y) into hexagonal cell number,
 * and returns (x,y) in the reference frame in which 
 * cell center is placed in origin. It is assumed that
 * spacing of cells in x direction is equal to 1. Side
 * to side cell size is therefore equal to 1 too.
 *
 * vvv - March 24, 2000
 */
{
  hex rsc;
  xy_to_hex(x, y, &rsc);
  hex_to_nh(rsc, nh);
  
  return;
}

void calin::math::hex_array::vvv::
nh_to_xy(int *nh, double *x, double *y)
/*
 * Converts hexagonal cell number into (x,y) coordinates.
 * Cell number must always be integer > 0. If it is <= 0,
 * it is assigned to 1.	It is assumed that spacing of cells 
 * in x direction is equal to 1. Side to side hexagonal 
 * cell size is therefore equal to 1 too.
 *
 * vvv - March 24, 2000
 */
{
  hex rsc;

  if(*nh < 1) *nh=1;

  nh_to_hex(*nh, &rsc);
  hex_to_xy(rsc, x, y);

  return;
}
