#ifndef Declaration_H
#define Declaration_H

#include <fstream>
#include <vector>

extern std::fstream Inputs, Outputs;

extern long nD,nV,nN;

extern int *Sv, *Ev;
extern double *Pj, *Qj, *Qv, *LIMv;

extern double *Cij, *CPij;

extern int* ISTAR;

extern double ALPHA, BETA, GAMMA, W_LOAD, W_DIS;

extern int nInit_Col_Rand_VRPs;
extern long Max_ESSPRC_nROUTE; 
extern bool HypoTenuseCutKey;
extern double TimeLim, Time;
#define SUP_FASHION  0 // 0 = all vehicles together/ 1 = one by one
#define ExtraDualBound  1 // Only used in SUP_FASHION = 0 to calculate dual bound for each node
#define nThreads 1
#define OptEPS .000001
#define EPS .0000001
#define INF 1e16
#define EqualityEPS .000001

#endif
