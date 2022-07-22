#include "gsa.h"
#include "simulator.h"
#include "tools.h"
#include "cellenum.h"
#include "psucc.h"
#include "reduction.h"
#include "dpenum.h"
#include <iostream>
#include <cmath>

NTL_CLIENT
using namespace std;

Vec<long> DP_enumeration_wrap(const long& dim, const mat_ZZ& B_init, const double & target2_d, const int& blocksize, \
	const int& TOURS, const long& M, int insertion)
{
	int PREC = RR::precision();
	long i, j;
	long m, n, seed;
	mat_ZZ B, B0;
	Mat<long> B_l; 
	mat_RR U;
	Mat<double> U_d;
	Vec<int> pi; //coordinate permutation
	Vec<double> r, B1norm2_d;
	vec_RR B1norm2;
	RR GSS, GSS0;
	double GSS_d;
	// dicrete pruning settings
	double NNRBound; // NNR expected value's upperbound
	double minBound, maxBound;
	Vec<Vec<int>> List;
	bool rerandom = true;
	long DP_rounds = 0;
	// time counter
	struct timespec starttime, endtime, time1, time2;
	// solution
	Vec<long> sol;
	Vec<int> sol_z;
	sol.SetLength(dim);
	sol_z.SetLength(dim);
	// for insertion
	Vec<cnode> Candidate;
	// status for fplll functions;
	int status = 0;

	n = dim;
	B = B_init;
	B0 = B_init;
	// preprocess
	status = LLL_default(B, 0.99, FT_DEFAULT, PREC);
	m = B.NumRows();
	n = B.NumCols();
	B1norm2.SetLength(m);
	B1norm2_d.SetLength(m);
	rePermutation(B);
	status |= LLL_default(B, 0.99, FT_DEFAULT, PREC);
	status |= BKZ_default(B, blocksize, FT_DEFAULT);
	// compute GS information
	ComputeGS(B, U, B1norm2);
	conv(B1norm2_d, B1norm2);
	conv(U_d, U);
	conv(B_l, B);
	GSS = SUM(B1norm2);
	conv(GSS_d, GSS);
	GSS0 = GSS;
	rerandom = true;
	DP_rounds = 0;
	clear(sol);
	// start main loops
	clock_gettime(CLOCK_MONOTONIC, &starttime);
	while(rerandom)
	{
		conv(B1norm2_d, B1norm2);
		DP_rounds ++;
		// === define discrete pruning param ===
		minBound = 0; 							//when tag is zero vector, f(Cell(t)) = 0
		maxBound = 0;
	    for(i=0; i<n; i++)
	    {
	        double zi = pow(M, 1/(double)n);
	        zi = ceil(zi);
	        maxBound += (zi+1) * zi * B1norm2_d[i];
	    }
		computePrunPara(pi, r, B1norm2_d, n); 	//r[i] = Bnorm2[pi[i]]
		NNRBound = computeBound(minBound, maxBound, pi, r, n , M);
		if(NNRBound <= 0)
		{
			cout<<"Upperbound is too small or there is something wrong."<<endl;
			return sol;
		}
		// === cell ENUM and Decode ===
		long counter;
		if(insertion){
			// List.SetLength(0);
			// cellENUM(List, counter, pi, NNRBound, r, n);
			// ENUM_discrete_strategy(sol, sol_z, Candidate, List, B_l, B1norm2_d, U_d, target2_d, n);
			DP_ENUM_ins_Wrap(counter, pi, NNRBound, r, n, sol, sol_z, Candidate, B_l, B1norm2_d, U_d, target2_d);
		}
		else {
			// List.SetLength(0);
			// cellENUM(List, counter, pi, NNRBound, r, n);
			// ENUM_discrete_postprocessing(sol, sol_z, List, B_l, B1norm2_d, U_d, target2_d, n);
			DP_ENUM_Wrap(counter, pi, NNRBound, r, n, sol, sol_z, B_l, B1norm2_d, U_d, target2_d);
		}
		// === check if solution is found ===
		if(!isZero(sol))
		{
			long len2;
			MUL(len2, sol, sol);
			cout << "Find solution: " << sol <<" with len^2 = "<<len2<<endl;

			rerandom = false;
			break;
		}
		// === Insertion and(or) Rerandomize ===
		long Clength = Candidate.length();
		if(Clength && insertion)
		{
			long min_h=n+1;
			cnode min_node = Candidate[0];
			for(i=1; i<Clength; i++)
			{
				if(Candidate[i].h_v < min_h)
				{
					min_h = Candidate[i].h_v;
					min_node = Candidate[i];
				}
				else if(Candidate[i].h_v == min_h)
				{
					if(Candidate[i].v_norm < min_node.v_norm)
					{
						min_node = Candidate[i];
					}
				}
			}
			// save old basis & GSO
			B0 = B;
			GSS0 = GSS;
			if(min_node.h_v > 0 && min_node.h_v < n/2) 
			 // good insertion vector existed, (min_node.h_v < n/2) can be modified to smaller thredshold
			{
				updateBasis(B, min_node, n);
				if(blocksize > 0)
					status |= BKZ_k_tours(B, blocksize, TOURS, FT_DEFAULT);
				else if(blocksize < 0)
					status |= BKZ_quick(B, -blocksize, FT_DEFAULT);
				// else LLL done in update
			}
			ComputeGS(B, U, B1norm2);
			GSS = SUM(B1norm2);
			if(GSS > conv<RR>("0.999")*GSS0) // insertion and update failed, rerandomize old basis
			{
				B = B0;
				rePermutation_keephead(B, n/5);

				if(blocksize > 0)
					status |= BKZ_k_tours(B, blocksize, TOURS, FT_DEFAULT);
				else
					status |= LLL_default(B, 0.99, FT_DEFAULT, PREC);

				ComputeGS(B, U, B1norm2);
				conv(B1norm2_d, B1norm2);
				conv(U_d, U);
				conv(B_l, B);
				GSS = SUM(B1norm2);
				conv(GSS_d, GSS);
			}
			// else: update success, get a basis with smaller GSS, continue;
			else
			{
				conv(B1norm2_d, B1norm2);
				conv(U_d, U);
				conv(B_l, B);
				conv(GSS_d, GSS);
			}
		}
		else
		{
			rePermutation_keephead(B, 5);

			if(blocksize > 0)
				status |= BKZ_k_tours(B, blocksize, TOURS, FT_DEFAULT);
			else
				status |= LLL_default(B, 0.99, FT_DEFAULT, PREC);

			ComputeGS(B, U, B1norm2);
			conv(B1norm2_d, B1norm2);
			conv(U_d, U);
			conv(B_l, B);
			GSS = SUM(B1norm2);
			conv(GSS_d, GSS);
		} // === END if(Clength && insertion)
	}
	clock_gettime(CLOCK_MONOTONIC, &endtime);
	double totalTime = timeDiff(starttime, endtime);
	cout<<"Total time for solving SVP: "<<totalTime<<endl;
	return sol;
}

costAll DP_optimization_warp(const mat_ZZ& B_instance, const long& n, const ZZ& det, const double& target_d, \
	long& M_op,  int& bkz_op,\
	int DP_SIM_FLAG, int TOURS)
{
	/*
	There are three methods provided when call NM method: 
	DP_SIM_FLAG = 0 -> only use simulated GS sequnce to estimate the running time;
	DP_SIM_FLAG = 1 -> use the BKZ_beta-reduced(2.0 ver) input basis `B_instance` to estimate p_succ;
	DP_SIM_FLAG = 2 -> execute one round of DP enum to record the time, and use the reduced `B_instance` to estimate p_succ;
	*/
	costAll cost;
	cost = DP_NelderMead(B_instance, n, det, target_d, M_op, bkz_op, DP_SIM_FLAG, TOURS, true);
	return cost;
}