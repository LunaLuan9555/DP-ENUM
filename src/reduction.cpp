#include "tools.h"
#include "reduction.h"

NTL_CLIENT
using namespace std;
using namespace fplll;

// ======================== BKZ & LLL ========================

int BKZ_one_tour(mat_ZZ &B, int blocksize, FloatType float_type, int MPFR_PREC)
{
	int status = 0;
	FloatType sel_ft = (float_type != FT_DEFAULT) ? float_type : FT_DOUBLE;
	if(sel_ft == FT_DOUBLE)
		status = BKZ_one_tour_f<FP_NR<double>>(B, blocksize, sel_ft);
	else if (sel_ft == FT_LONG_DOUBLE)
		status = BKZ_one_tour_f<FP_NR<long double>>(B, blocksize, sel_ft);
	else if (sel_ft == FT_MPFR)
	{
		int old_prec = FP_NR<mpfr_t>::set_prec(MPFR_PREC);
		status = BKZ_one_tour_f<FP_NR<mpfr_t>>(B, blocksize, sel_ft);
		FP_NR<mpfr_t>::set_prec(old_prec);
	}
	return status;
}

int BKZ_k_tours(mat_ZZ &B, int blocksize, int tours, FloatType float_type, int MPFR_PREC)
{
	int status = 0;
	FloatType sel_ft = (float_type != FT_DEFAULT) ? float_type : FT_DOUBLE;
	if(sel_ft == FT_DOUBLE)
		status = BKZ_k_tours_f<FP_NR<double>>(B, blocksize, tours, sel_ft);
	else if (sel_ft == FT_LONG_DOUBLE)
		status = BKZ_k_tours_f<FP_NR<long double>>(B, blocksize, tours, sel_ft);
	else if (sel_ft == FT_MPFR)
	{
		int old_prec = FP_NR<mpfr_t>::set_prec(MPFR_PREC);
		status = BKZ_k_tours_f<FP_NR<mpfr_t>>(B, blocksize, tours, sel_ft);
		FP_NR<mpfr_t>::set_prec(old_prec);
	}
	return status;
}


int BKZ_quick(mat_ZZ &B, int blocksize, FloatType float_type, int MPFR_PREC)
{
	ZZ_mat<mpz_t> A;
	convtpye_A2B(B, A);
	//use bkz to solve and record the time 
	clock_t start_time = clock();

	int status = 0;
	// set bkz parameters
	vector<Strategy> strategies = load_strategies_json(PRUN_PATH);
	BKZParam bkz_param(blocksize, strategies); //blocksize, strategy
	// bkz_param.flags = BKZ_VERBOSE;
	bkz_param.flags = BKZ_AUTO_ABORT;
	// bkz reduction
	status = bkz_reduction(&A, NULL, bkz_param, float_type, MPFR_PREC);
	if(status)
		cout<<"!!! BKZ ERROR"<<endl;
	convtpye_A2B(A, B);
	return status;
}

int BKZ_default(mat_ZZ &B, int blocksize, FloatType float_type, int MPFR_PREC)
{
	ZZ_mat<mpz_t> A;
	convtpye_A2B(B, A);
	//use bkz to solve and record the time 
	clock_t start_time = clock();

	int status = 0;
	// set bkz parameters
	vector<Strategy> strategies = load_strategies_json(PRUN_PATH);
	BKZParam bkz_param(blocksize, strategies); //blocksize, strategy
	// bkz_param.flags = BKZ_VERBOSE;
	// bkz reduction
	status = bkz_reduction(&A, NULL, bkz_param, float_type, MPFR_PREC);
	if(status)
		cout<<"!!! BKZ ERROR"<<endl;
	convtpye_A2B(A, B);
	return status;
}

int BKZ_fullENUM(mat_ZZ &B, int blocksize, FloatType float_type, int MPFR_PREC)
{
	ZZ_mat<mpz_t> A;
	convtpye_A2B(B, A);
	//use bkz to solve and record the time 
	clock_t start_time = clock();

	int status = 0;
	// set bkz parameters
	vector<Strategy> strategies;
	BKZParam bkz_param(blocksize, strategies); //blocksize, strategy
	// bkz_param.flags = BKZ_VERBOSE;
	// bkz reduction
	status = bkz_reduction(&A, NULL, bkz_param, float_type, MPFR_PREC);
	if(status)
		cout<<"!!! BKZ ERROR"<<endl;
	convtpye_A2B(A, B);
	return status;
}

int BKZ_linearPrune(mat_ZZ &B, int blocksize, FloatType float_type, int MPFR_PREC)
{
	ZZ_mat<mpz_t> A;
	convtpye_A2B(B, A);
	//use bkz to solve and record the time 
	clock_t start_time = clock();

	int status = 0;
	// set bkz parameters
	vector<Strategy> strategies;
	for (long b = 0; b < blocksize; b++)
	{
		Strategy strategy = Strategy::EmptyStrategy(b);
		if (b == 10)
		{
			strategy.preprocessing_block_sizes.emplace_back(5);
		}
		else if (b == 20)
		{
			strategy.preprocessing_block_sizes.emplace_back(10);
		}
		else if (b == 30)
		{
			strategy.preprocessing_block_sizes.emplace_back(15);
		}
		strategies.emplace_back(std::move(strategy));
	}
	Strategy strategy;
	strategy.pruning_parameters.emplace_back(
		PruningParams::LinearPruningParams(blocksize, blocksize / 2));
		strategies.emplace_back(std::move(strategy));

	BKZParam bkz_param(blocksize, strategies);
	// bkz_param.flags = BKZ_DEFAULT;
	// bkz_param.flags = BKZ_VERBOSE;
	// bkz reduction
	status = bkz_reduction(&A, NULL, bkz_param, float_type, MPFR_PREC);
	if(status)
		cout<<"!!! BKZ ERROR"<<endl;
	convtpye_A2B(A, B);
	return status;
}

int LLL_default(mat_ZZ &B, double delta, FloatType float_type, int MPFR_PREC)
{
	int status = 0;
	ZZ_mat<mpz_t> A;
	convtpye_A2B(B, A);
	status = lll_reduction(A, delta, LLL_DEF_ETA, LM_PROVED, float_type, MPFR_PREC);
	convtpye_A2B(A, B);
	return status;
}

int LLL_quick(mat_ZZ &B, double delta)
// use FT_DEFAULT = double precision
{
	int status = 0;
	ZZ_mat<mpz_t> A;
	convtpye_A2B(B, A);
	status = lll_reduction(A, delta);
	convtpye_A2B(A, B);
	return status;
}

