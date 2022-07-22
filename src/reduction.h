#ifndef REDUCTION__H
#define REDUCTION__H

#include "tools.h"

NTL_CLIENT
using namespace std;
using namespace fplll;


// BKZ as preprocessing, namely only run one tour of BKZ
int BKZ_one_tour(mat_ZZ &B, int blocksize, FloatType float_type = FT_DEFAULT, int MPFR_PREC=150);

int BKZ_k_tours(mat_ZZ &B, int blocksize, int tours=1, FloatType float_type = FT_DEFAULT, int MPFR_PREC=150);

// BKZ 2.0 reduction routine
int BKZ_default(mat_ZZ& B, int blocksize, FloatType float_type = FT_MPFR, int MPFR_PREC = 150);

// BKZ with early-abort
int BKZ_quick(mat_ZZ& B, int blocksize, FloatType float_type = FT_DEFAULT, int MPFR_PREC = 150);

// BKZ with linear pruning
int BKZ_linearPrune(mat_ZZ &B, int blocksize, FloatType float_type = FT_MPFR, int MPFR_PREC = 150);

// BKZ with full enum
int BKZ_fullENUM(mat_ZZ &B, int blocksize, FloatType float_type = FT_MPFR, int MPFR_PREC = 150);

int LLL_default(mat_ZZ &B, double delta, FloatType float_type = FT_MPFR, int MPFR_PREC = 150);

int LLL_quick(mat_ZZ &B, double delta);

// implementations
template <class FT>  //FT flote type for GSO 
int BKZ_one_tour_f(mat_ZZ &B0, int blocksize, int sel_ft) 
{
	//gamma = 1.05 by default
	ZZ_mat<mpz_t> B;
	convtpye_A2B(B0, B);
	int status = 0;
	long n, m, i, j;
	n = B0.NumRows();  //dimension of lattice
	m = B0.NumCols();  //dimension of Euclidean space
	int numrows = (int)n;
	// setting enum param
	vector<Strategy> strategies = load_strategies_json(PRUN_PATH);
	BKZParam param(blocksize, strategies); //blocksize, strategy
	param.flags = BKZ_GH_BND;

	// generate MatGSO object
	int gso_flags = 0;
	if (sel_ft == FT_DOUBLE || sel_ft == FT_LONG_DOUBLE)
		gso_flags |= GSO_ROW_EXPO;

	ZZ_mat<long> Bl;
	ZZ_mat<mpz_t> empty_mat;
	ZZ_mat<mpz_t> &u = empty_mat;
	ZZ_mat<mpz_t> &u_inv = empty_mat;

	if (convert<long, mpz_t>(Bl, B, 10))
	{
		ZZ_mat<long> ul;
		convert<long, mpz_t>(ul, u, 0);
		ZZ_mat<long> ul_inv;
		convert<long, mpz_t>(ul_inv, u_inv, 0);
		MatGSO<Z_NR<long>, FT> m_gso(Bl, ul, ul_inv, gso_flags);
		m_gso.update_gso();
		LLLReduction<Z_NR<long>, FT> lll_obj(m_gso, LLL_DEF_DELTA, LLL_DEF_ETA, LLL_DEFAULT);
		lll_obj.lll();
		BKZReduction<Z_NR<long>, FT> bkz_obj(m_gso, lll_obj, param);
		bkz_obj.tour(0, numrows, param, 0, numrows);

		convert<mpz_t, long>(B, Bl, 0);
		status = bkz_obj.status;
	}
	else
	{
		MatGSO<Z_NR<mpz_t>, FT> m_gso(B, u, u_inv, gso_flags);
		LLLReduction<Z_NR<mpz_t>, FT> lll_obj(m_gso, LLL_DEF_DELTA, LLL_DEF_ETA, LLL_DEFAULT);
		lll_obj.lll();
		BKZReduction<Z_NR<mpz_t>, FT> bkz_obj(m_gso, lll_obj, param);
		bkz_obj.tour(0, numrows, param, 0, numrows);
		status = bkz_obj.status;
	}
	convtpye_A2B(B, B0);
	return status;
}

template <class FT>  //FT flote type for GSO 
int BKZ_k_tours_f(mat_ZZ &B0, int blocksize, int tours, int sel_ft) 
{
	//gamma = 1.05 by default
	ZZ_mat<mpz_t> B;
	convtpye_A2B(B0, B);
	int status = 0;
	long n, m;
	n = B0.NumRows();  //dimension of lattice
	m = B0.NumCols();  //dimension of Euclidean space
	int numrows = (int)n;
	// setting enum param
	vector<Strategy> strategies = load_strategies_json(PRUN_PATH);
	BKZParam param(blocksize, strategies); //blocksize, strategy
	param.flags = BKZ_GH_BND;

	// generate MatGSO object
	int gso_flags = 0;
	if (sel_ft == FT_DOUBLE || sel_ft == FT_LONG_DOUBLE)
		gso_flags |= GSO_ROW_EXPO;

	ZZ_mat<long> Bl;
	ZZ_mat<mpz_t> empty_mat;
	ZZ_mat<mpz_t> &u = empty_mat;
	ZZ_mat<mpz_t> &u_inv = empty_mat;

	if (convert<long, mpz_t>(Bl, B, 10))
	{
		ZZ_mat<long> ul;
		convert<long, mpz_t>(ul, u, 0);
		ZZ_mat<long> ul_inv;
		convert<long, mpz_t>(ul_inv, u_inv, 0);
		MatGSO<Z_NR<long>, FT> m_gso(Bl, ul, ul_inv, gso_flags);
		m_gso.update_gso();
		LLLReduction<Z_NR<long>, FT> lll_obj(m_gso, LLL_DEF_DELTA, LLL_DEF_ETA, LLL_DEFAULT);
		lll_obj.lll();
		BKZReduction<Z_NR<long>, FT> bkz_obj(m_gso, lll_obj, param);
		for(int i = 0; i<tours; i++)
		{
			bkz_obj.tour(0, numrows, param, 0, numrows);
		}

		convert<mpz_t, long>(B, Bl, 0);
		status = bkz_obj.status;
	}
	else
	{
		MatGSO<Z_NR<mpz_t>, FT> m_gso(B, u, u_inv, gso_flags);
		LLLReduction<Z_NR<mpz_t>, FT> lll_obj(m_gso, LLL_DEF_DELTA, LLL_DEF_ETA, LLL_DEFAULT);
		lll_obj.lll();
		BKZReduction<Z_NR<mpz_t>, FT> bkz_obj(m_gso, lll_obj, param);
		for(int i = 0; i<tours; i++)
		{
			bkz_obj.tour(0, numrows, param, 0, numrows);
		}
		status = bkz_obj.status;
	}
	convtpye_A2B(B, B0);
	return status;
}

#endif