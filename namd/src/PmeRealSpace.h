/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PME_REAL_SPACE_H__
#define PME_REAL_SPACE_H__

#include "PmeBase.h"
#include "Vector.h"

class PmeRealSpace {
  
public:
  PmeRealSpace(PmeGrid grid, int natoms);
  ~PmeRealSpace();

  void fill_charges(double **q_arr, double **q_arr_list, int &q_arr_count,
                       int &stray_count, char *f_arr, char *fz_arr, PmeParticle p[]); 
  void compute_forces(const double * const *q_arr, const PmeParticle p[], 
                      Vector f[]);
                      
  void compute_forces_order4_partial(int first, int last, const double * const *q_arr, const PmeParticle p[], 
                      Vector f[]);                      
private:
  void fill_charges_order4(double **q_arr, double **q_arr_list, int &q_arr_count,
                       int &stray_count, char *f_arr, char *fz_arr, PmeParticle p[]); 
  void compute_forces_order4(const double * const *q_arr, const PmeParticle p[], 
                      Vector f[]);
  void fill_b_spline(PmeParticle p[]);

  const int N;
  const PmeGrid myGrid;
  double *M, *dM;
};


#endif

