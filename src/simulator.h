#ifndef _SIMULATOR_H_
#define _SIMULATOR_H_

#include "tools.h"
#include <limits.h>

#define NM_MAX numeric_limits<double>::max()

NTL_CLIENT
using namespace std;

// qmean[0] = everage q of bkz_10
// qmean from BKZ_10 to BKZ_44
extern double qmean[35];

typedef struct CostAll{
	double rounds;
	double oneround;
    double preproc;
}costAll;

costAll DP_NelderMead(const mat_ZZ& B0, const long & n, const ZZ& det, const double &target_d, long &M_op,  int& bkz_op, \
    int DP_SIM_FLAG=0, int k=6, bool NM_VERBOSE = true);

costAll DP_simulator_warp(const mat_ZZ& B_instance, const long& n, const ZZ& det, const double& target_d, \
    long M, int bkz,
    int DP_SIM_FLAG=0, int TOURS=6);

costAll DP_sim_pure(const ZZ& det, const double& target_d, const long& M, const long& n, int k, int bkz);

costAll DP_sim_hybrid(const mat_ZZ& B0, const double& target_d, const long& M, const long& n, int k, int bkz, \
    FloatType float_type = FT_DEFAULT);

costAll DP_sim_basis(const mat_ZZ& B0, const double& target_d, const long& M, const long& n, int k, int bkz, FloatType float_type = FT_DEFAULT);

void DP_ENUM_Wrap_sim(long & counter, const Vec<int>& pi, double Bound, const Vec<double>& r, long n, \
    Vec<long>& sol, Vec<int>& sol_z, const Mat<long>& B, const Vec<double>& B1norm2_d, \
    const Mat<double>& U, const double target2);

// costAll findMinCost_pure(const ZZ& det, const double& target_d, const long& n, \
//     long& M, int &bkz, int& k);

#endif