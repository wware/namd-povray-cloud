/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEPME_H
#define COMPUTEPME_H

#include "ComputeHomePatches.h"
#include "PmeBase.h"
#include "NamdTypes.h"

class PmeRealSpace;
class ComputeMgr;
class SubmitReduction;
class PmeGridMsg;
class ComputePmeMgr;

class ComputePme : public ComputeHomePatches {
public:
  ComputePme(ComputeID c);
  virtual ~ComputePme();
  void doWork();
  void sendData(int, int*, int*, int*);
  void sendPencils();
  void copyResults(PmeGridMsg *);
  void copyPencils(PmeGridMsg *);
  void ungridForces();
  void setMgr(ComputePmeMgr *mgr) { myMgr = mgr; }
#if CMK_PERSISTENT_COMM
  void setup_recvgrid_persistent();
#endif
 private:
  PmeGrid myGrid;
  int qsize, fsize, bsize;
  int alchFepOn, alchThermIntOn, lesOn, lesFactor, pairOn, selfOn, numGrids;
  int alchDecouple;
  BigReal alchElecLambdaStart;
  BigReal elecLambdaUp;
  BigReal elecLambdaDown;
  
  double **q_arr;
  double **q_list;
  int q_count;
  char *f_arr;
  char *fz_arr;
  PmeReduction evir[PME_MAX_EVALS];
  SubmitReduction *reduction;
  SubmitReduction *amd_reduction;
  int strayChargeErrors;
  int resultsRemaining;
  PmeRealSpace *myRealSpace[PME_MAX_EVALS];
  int numLocalAtoms;
  PmeParticle *localData;
  unsigned char *localPartition;
  int numGridAtoms[PME_MAX_EVALS];
  PmeParticle *localGridData[PME_MAX_EVALS];
  ComputePmeMgr *myMgr;

#if CMK_PERSISTENT_COMM
  PersistentHandle   *recvGrid_handle;
#endif

};

#endif

