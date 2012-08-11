/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "PmeBase.h"
#include <math.h>

void compute_b_spline(double *frac, double *M, double *dM, int order) {
  int j, n;
  double x,y,z,x1,y1,z1, div;
  double *Mx, *My, *Mz, *dMx, *dMy, *dMz;
  Mx=M-1; My=M+order-1; Mz=M+2*order-1;
  dMx=dM-1; dMy =dM+order-1; dMz=dM+2*order-1;
  x=frac[0];
  y=frac[1];
  z=frac[2];
  x1=1.0-x; y1=1.0-y; z1=1.0-z;
  /* Do n=3 case first */
  Mx[1]=.5*x1*x1;
  Mx[2]=x1*x + .5;
  Mx[3]=0.5*x*x;
  Mx[order]=0.0;
  My[1]=.5*y1*y1;
  My[2]=y1*y + .5;
  My[3]=0.5*y*y;
  My[order]=0.0;
  Mz[1]=.5*z1*z1;
  Mz[2]=z1*z + .5;
  Mz[3]=0.5*z*z;
  Mz[order]=0.0;

  /* Recursively fill in the rest.  */
  for (n=4; n<=order-1; n++) {
    double div=1.0/(n-1);
    int j;
    Mx[n] = x*div*Mx[n-1];
    My[n] = y*div*My[n-1];
    Mz[n] = z*div*Mz[n-1];
    for (j=1; j<=n-2; j++) {
      Mx[n-j] = ((x+j)*Mx[n-j-1] + (n-x-j)*Mx[n-j])*div;
      My[n-j] = ((y+j)*My[n-j-1] + (n-y-j)*My[n-j])*div;
      Mz[n-j] = ((z+j)*Mz[n-j-1] + (n-z-j)*Mz[n-j])*div;
    }
    Mx[1] *= (1.0-x)*div;
    My[1] *= (1.0-y)*div;
    Mz[1] *= (1.0-z)*div;
  }
  /* Now do the derivatives  */
  dMx[1]=-Mx[1]; dMy[1]=-My[1]; dMz[1]=-Mz[1];
  for (j=2; j <= order; j++) {
    dMx[j] = Mx[j-1] - Mx[j];
    dMy[j] = My[j-1] - My[j];
    dMz[j] = Mz[j-1] - Mz[j];
  }
  /* Now finish the job!    */
  div=1.0/(order-1);
  Mx[order] = x*div*Mx[order-1];
  My[order] = y*div*My[order-1];
  Mz[order] = z*div*Mz[order-1];
  for (j=1; j<=order-2; j++) {
    Mx[order-j] = ((x+j)*Mx[order-j-1] + (order-x-j)*Mx[order-j])*div;
    My[order-j] = ((y+j)*My[order-j-1] + (order-y-j)*My[order-j])*div;
    Mz[order-j] = ((z+j)*Mz[order-j-1] + (order-z-j)*Mz[order-j])*div;
  }
  Mx[1] *= (1.0-x)*div;
  My[1] *= (1.0-y)*div;
  Mz[1] *= (1.0-z)*div;
}


