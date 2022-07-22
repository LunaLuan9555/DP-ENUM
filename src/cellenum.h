#ifndef CELLENUM__H
#define CELLENUM__H

#include "tools.h"
#include "reduction.h"

NTL_CLIENT
using namespace std;
using namespace fplll;


// ------------------------ cellEnum.cpp --------------------------
typedef struct CandidateNode{
	Vec<int> vec_z;
	Vec<long> vec_v;
	long v_norm;
	long h_v;
}cnode;

bool lastOdd(const Vec<int>& vout, long n);

void updateBasis(mat_ZZ &B, const cnode& node, long n);

void cellENUM(Vec<Vec<int>>& List, long & counter, const Vec<int>& pi, double Bound, const Vec<double>& r, long n);

void ENUM_discrete_strategy(Vec<long>& sol, Vec<int>& sol_z, Vec<cnode>& Candidate, \
	const Vec<Vec<int>> List, const Mat<long>& B, const Vec<double>& B1norm2_d, \
	const Mat<double>& U, const double& target2, const long & n);

void ENUM_discrete_postprocessing(Vec<long>& sol, Vec<int>& sol_z, const Vec<Vec<int>> List, const Mat<long>& B, \
	const Vec<double>& B1norm2_d, const Mat<double>& U, const double& target2, const long & n);

// --- what's NEW ---
int DP_ENUM_ins_decode(Vec<long>& sol, Vec<int>& sol_z, Vec<cnode>& Candidate, const Vec<int>& z, const Mat<long>& B, \
	const Vec<double>& B1norm2_d, const Vec<double>& B1norm2_d_relax,\
	const Mat<double>& U, const double & target2, const double & target2_relax, const long & n);

void DP_ENUM_ins_Wrap(long & counter, const Vec<int>& pi, double Bound, const Vec<double>& r, long n, \
	Vec<long>& sol, Vec<int>& sol_z, Vec<cnode>& Candidate, const Mat<long>& B, const Vec<double>& B1norm2_d, \
	const Mat<double>& U, const double& target2);

int DP_ENUM_decode(Vec<long>& sol, Vec<int>& sol_z, const Vec<int>& z, const Mat<long>& B, const Vec<double>& B1norm2_d, \
	const Mat<double>& U, const double& target2, const long & n);

void DP_ENUM_Wrap(long & counter, const Vec<int>& pi, double Bound, const Vec<double>& r, long n, \
	Vec<long>& sol, Vec<int>& sol_z, const Mat<long>& B, const Vec<double>& B1norm2_d, \
	const Mat<double>& U, const double target2);

long cellENUM_probe(const Vec<int>& pi, double Bound, const Vec<double>& r, long n, long M);

double computeBound(double Bleft, double Bright, const Vec<int>& pi, const Vec<double>& r, long n, long M);

void computePrunPara(Vec<int>& pi, Vec<double>& r, const Vec<double>& Bnorm2, long n);

void ENUM_discrete_postprocessing_M(long M, Vec<long>& sol, Vec<int>& sol_z, const Vec<Vec<int>> List, \
	const Mat<long>& B, const Vec<double>& B1norm2_d, \
	const Mat<double>& U, const double target2, long n);

void ENUM_discrete_strategy_M(long M, Vec<long>& sol, Vec<int>& sol_z, Vec<cnode>& Candidate, const Vec<Vec<int>> List, \
	const Mat<long>& B, const Vec<double>& B1norm2_d, \
	const Mat<double>& U, const double target2, long n);

#endif