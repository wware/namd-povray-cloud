/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PME_BASE_H__
#define PME_BASE_H__

#include <math.h>
#include <charm++.h>
#include "MathArray.h"
#include "Vector.h"
#include "Lattice.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct PmeGrid {
  int K1, K2, K3;
  int dim2, dim3;
  int order;
  int block1, block2, block3;
  int xBlocks, yBlocks, zBlocks;
};

struct PmeParticle {
  double x, y, z;
  double cg;
};

void compute_b_spline(double frac[3], double *M, double *dM, int order);

#define PME_MAX_EVALS 255 
typedef MathArray<double,7> PmeReduction;

#ifndef SQRT_PI
#define SQRT_PI 1.7724538509055160273 /* mathematica 15 digits*/
#endif


static inline void scale_coordinates(PmeParticle p[], int N, Lattice lattice, PmeGrid grid) {
  Vector origin = lattice.origin();
  Vector recip1 = lattice.a_r();
  Vector recip2 = lattice.b_r();
  Vector recip3 = lattice.c_r();
  double ox = origin.x;
  double oy = origin.y;
  double oz = origin.z;
  double r1x = recip1.x;
  double r1y = recip1.y;
  double r1z = recip1.z;
  double r2x = recip2.x;
  double r2y = recip2.y;
  double r2z = recip2.z;
  double r3x = recip3.x;
  double r3y = recip3.y;
  double r3z = recip3.z;
  int K1 = grid.K1;
  int K2 = grid.K2;
  int K3 = grid.K3;

  for (int i=0; i<N; i++) {
    double px = p[i].x - ox;
    double py = p[i].y - oy;
    double pz = p[i].z - oz;
    double sx = px*r1x + py*r1y + pz*r1z;
    double sy = px*r2x + py*r2y + pz*r2z;
    double sz = px*r3x + py*r3y + pz*r3z;
    p[i].x = K1 * ( sx - floor(sx) );
    p[i].y = K2 * ( sy - floor(sy) );
    p[i].z = K3 * ( sz - floor(sz) );
    //  Check for rare rounding condition where K * ( 1 - epsilon ) == K
    //  which was observed with g++ on Intel x86 architecture.
    if ( p[i].x == K1 ) p[i].x = 0;
    if ( p[i].y == K2 ) p[i].y = 0;
    if ( p[i].z == K3 ) p[i].z = 0;
  }
}


static inline void scale_forces(Vector f[], int N, Lattice &lattice) {
  Vector recip1 = lattice.a_r();
  Vector recip2 = lattice.b_r();
  Vector recip3 = lattice.c_r();
  double r1x = recip1.x;
  double r1y = recip1.y;
  double r1z = recip1.z;
  double r2x = recip2.x;
  double r2y = recip2.y;
  double r2z = recip2.z;
  double r3x = recip3.x;
  double r3y = recip3.y;
  double r3z = recip3.z;

  for (int i=0; i<N; i++) {
    double f1 = f[i].x;
    double f2 = f[i].y;
    double f3 = f[i].z;
    f[i].x = f1*r1x + f2*r2x + f3*r3x;
    f[i].y = f1*r1y + f2*r2y + f3*r3y;
    f[i].z = f1*r1z + f2*r2z + f3*r3z;
  }
}

#endif
