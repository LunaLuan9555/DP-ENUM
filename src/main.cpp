#include <unistd.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "tools.h"
#include "reduction.h"
#include "psucc.h"
#include "cellenum.h"
#include "gsa.h"
#include "simulator.h"
#include "dpenum.h"

// #define RERANDOM 30
#define RADIUS_RELAXATION_FACTOR2	1.2
#define INSERTION_FACTOR_DELTA		0.99

#define CPUFREQ 2.1e9

NTL_CLIENT
using namespace std;
using namespace fplll;


int main(int argc,char **argv)
{
	// RR::SetOutputPrecision(40);
	int PREC = RR::precision();
	
	long i, j, m, n, status;
	// ============ read  parameters =============.
	n = 60;
	long seed = -1;
	int blocksize = 16;
	int TOURS = 6;
	long REPEAT = 1;
	int insertion = 0;
	double alpha = 1.05;
	int SIMTYPE=-1;
	RR GHLen, target2;
	double target2_d, target_d;
	ZZ det;
	long M_op;
	int bkz_op;
	costAll cost_est;
	mat_ZZ B;

	PARSE_MAIN_ARGS{
		MATCH_MAIN_ARGID("--dim", n);
		MATCH_MAIN_ARGID("--seed", seed);
		MATCH_MAIN_ARGID("--bkz", blocksize);
		MATCH_MAIN_ARGID("--M", M);
		MATCH_MAIN_ARGID("--tours", TOURS);
		MATCH_MAIN_ARGID("--insert", insertion);
		MATCH_MAIN_ARGID("--approx", alpha);
		MATCH_MAIN_ARGID("--simtype", SIMTYPE);
		SYNTAX();
	}

	char *fp = (char *)calloc(150, sizeof(char));
	char *outinfo = (char *)calloc(500, sizeof(char));
	
	if(seed<0)
		seed = rand()%8000 + 1000;
	GenerateBasis(B, n, seed); //generate SVP challenge by default (take care of the psedo-random number generator of NTL!)
	status = LLL_default(B, 0.99, FT_DEFAULT, PREC);
	if(n>100)
	{
		BKZ_default(B, 21, FT_DEFAULT);
	}

	m = B.NumRows();
	n = B.NumCols();
	GHLen = GHL(B); 
	target2 = conv<RR>(alpha) * conv<RR>(alpha) * GHLen * GHLen;
	conv(target2_d, target2);
	conv(target_d, SqrRoot(target2));
	det = SqrRoot(determinant( B*transpose(B)));
	if(SIMTYPE >= 0)
	{
		// ===== DP simulator =====
		cost_est = DP_optimization_warp(B, n, det, target_d, M_op, bkz_op, SIMTYPE);
		sprintf(outinfo, "The optimal parameters is : M = %ld, beta = %d. Expected time =%lf",\
			n, M_op, bkz_op, cost_est.preproc+cost_est.oneround*cost_est.rounds);
		cout << outinfo << endl;
	}
	else
	{
		// ===== DP solver ======
		DP_enumeration_wrap(n, B, target2_d, blocksize, TOURS, M); 
	}
	
	// =====================
	return 0;
}

