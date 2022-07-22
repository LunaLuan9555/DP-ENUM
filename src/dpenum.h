#ifndef _DPENUM_H_
#define _DPENUM_H_

NTL_CLIENT
using namespace std;

Vec<long> DP_enumeration_wrap(const long& dim, const mat_ZZ& B_init, const double & target2_d, \
	const int& blocksize, const int& TOURS,\
	const long& M, int insertion=0);


costAll DP_optimization_warp(const mat_ZZ& B_instance, const long& n, const ZZ& det, const double& target_d, \
	long& M_op,  int& bkz_op,\
	int DP_SIM_FLAG=0, int TOURS=6);

#endif