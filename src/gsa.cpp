#include "gsa.h"
#include "tools.h"
#include "cellenum.h"
#include "psucc.h"
#include "simulator.h"
#include <iostream>
#include <cmath>


NTL_CLIENT
using namespace std;

// ============= Simulator ==================
double SimulateCost_q(const ZZ& det, const int bkz, const long n, const double target_d, const long M)
{
    Vec<double> GS2_fitted;
    double q, logvol, logb0;
    GS2_fitted.SetLength(n);
    if(bkz>=10 && bkz <45)
        q = qmean[bkz-10];
    else
        q = 1-exp(-0.008668967071869243*bkz-3.405597365733596);
 // cout<< "Test in [SimulateCost_q]: q = "<<q<<endl;
    logvol = conv<double>(log(det));
    logb0 = (logvol - (double)n*(n-1)/2*log(q)) / (double)n;
    // === padding GS2 ======
    GS2_fitted[0] = exp(2*logb0);
    for(long i=1; i<n; i++)
    {
        GS2_fitted[i] = exp(2*(log(q)*i + logb0));
    }
    // === calculate p_succ ==============
    double NNRBound; // NNR expected value's upperbound
    double minBound, maxBound;
    Vec<Vec<int>> List;
    Vec<int> pi; //coordinate permutation
    Vec<double> r;
    double GSS_d;
    GSS_d= SUM(GS2_fitted);
    // 
    minBound = 0;                           //when tag is zero vector, f(Cell(t)) = 0
    maxBound = 10.0 * GSS_d;                //heuristically set. I donnot know how to set the upperbound.
    computePrunPara(pi, r, GS2_fitted, n);   //r[i] = GS2[pi[i]]
    NNRBound = computeBound(minBound, maxBound, pi, r, n , M);
    List.SetLength(0);
    long counter;
    cellENUM(List, counter, pi, NNRBound, r, n);
 // cout<< "Test in [SimulateCost_q]: enum vector(s) = "<<counter<<endl;
    return P_succ_wrap(List, n, GS2_fitted, target_d, counter);
}

// double SimulateCost_q_even(const ZZ& det, const int bkz, const long n, const double target_d, const long M)
// {
//     Vec<double> GS2_fitted;
//     double q, logvol, logb0;
//     GS2_fitted.SetLength(n);
//     q = 1-exp(-0.008668967071869243*bkz-3.405597365733596);
//  // cout<< "Test in [SimulateCost_q]: q = "<<q<<endl;
//     logvol = conv<double>(log(det));
//     logb0 = (logvol - (double)n*(n-1)/2*log(q)) / (double)n;
//     // === padding GS2 ======
//     GS2_fitted[0] = exp(2*logb0);
//     for(long i=1; i<n; i++)
//     {
//         GS2_fitted[i] = exp(2*(log(q)*i + logb0));
//     }
//     // === calculate p_succ ==============
//     double NNRBound; // NNR expected value's upperbound
//     double minBound, maxBound;
//     Vec<Vec<int>> List;
//     Vec<int> pi; //coordinate permutation
//     Vec<double> r;
//     double GSS_d;
//     GSS_d= SUM(GS2_fitted);
//     // 
//     minBound = 0;                           //when tag is zero vector, f(Cell(t)) = 0
//     maxBound = 10.0 * GSS_d;                //heuristically set. I donnot know how to set the upperbound.
//     computePrunPara(pi, r, GS2_fitted, n);   //r[i] = GS2[pi[i]]
//     NNRBound = computeBound(minBound, maxBound, pi, r, n , M);
//     List.SetLength(0);
//     long counter;
//     cellENUM(List, counter, pi, NNRBound, r, n);
//     makeListEven(List, n);
//  // cout<< "Test in [SimulateCost_q]: enum vector(s) = "<<counter<<endl;
//     return P_succ_wrap(List, n, GS2_fitted, target_d, counter);
// }


double SimulateCost_GSA(const Vec<double> & B1norm2_d, const double target_d, const long M)
{
    double NNRBound; // NNR expected value's upperbound
    double minBound, maxBound;
    Vec<Vec<int>> List;
    Vec<int> pi; //coordinate permutation
    Vec<double> r;
    long n;
    double GSS_d;
    // === fit GSnorm ==========
    Vec<double> GS2_fitted;
    GSA_fit(GS2_fitted, B1norm2_d);
    GSS_d= SUM(GS2_fitted);
    // === define discrete pruning param ===
    n = B1norm2_d.length();
    minBound = 0;                           //when tag is zero vector, f(Cell(t)) = 0
    maxBound = 10.0 * GSS_d;                //heuristically set. I donnot know how to set the upperbound.
    computePrunPara(pi, r, GS2_fitted, n);   //r[i] = GS2[pi[i]]
    NNRBound = computeBound(minBound, maxBound, pi, r, n , M);
    List.SetLength(0);
    long counter;
    cellENUM(List, counter, pi, NNRBound, r, n);
 // cout<< "Test in [SimulateCost_GSA]: enum vector(s) = "<<counter<<endl;
    return P_succ_wrap(List, n, GS2_fitted, target_d, counter);
}

double SimulateCost_GSA_even(const Vec<double> & B1norm2_d, const double target_d, const long M)
{
    double NNRBound; // NNR expected value's upperbound
    double minBound, maxBound;
    Vec<Vec<int>> List;
    Vec<int> pi; //coordinate permutation
    Vec<double> r;
    long n;
    double GSS_d;
    // === fit GSnorm ==========
    Vec<double> GS2_fitted;
    GSA_fit(GS2_fitted, B1norm2_d);
    GSS_d= SUM(GS2_fitted);
    // === define discrete pruning param ===
    n = B1norm2_d.length();
    minBound = 0;                           //when tag is zero vector, f(Cell(t)) = 0
    maxBound = 10.0 * GSS_d;                //heuristically set. I donnot know how to set the upperbound.
    computePrunPara(pi, r, GS2_fitted, n);   //r[i] = GS2[pi[i]]
    NNRBound = computeBound(minBound, maxBound, pi, r, n , M);
    List.SetLength(0);
    long counter;
    cellENUM(List, counter, pi, NNRBound, r, n);
    makeListEven(List, n);
 // cout<< "Test in [SimulateCost_GSA]: enum vector(s) = "<<counter<<endl;
    return P_succ_wrap(List, n, GS2_fitted, target_d, counter);
}
// =================== GSA ==================
void GSA_fit(Vec<double> & GS_fit, const Vec<double> & GS2)
{
	long n, i;
	n = GS2.length();
	GS_fit.SetLength(n);
	Vec<double> X, Y, coeff;
	X.SetLength(n);
	Y.SetLength(n);
    coeff.SetLength(5);
    clear(coeff);
	for(i=0; i<n; i++)
	{
		X[i] = (double)i;
		Y[i] = log(GS2[i])/2.0; //ln(||B_i^*||)
	}
    EMatrix(X, Y, 2, coeff); //linear fit
	for(i=0; i<n; i++)
	{
		GS_fit[i] = exp(2*(coeff[1]+coeff[2]*X[i]));
	}
}

// ============= linear regressiong ===================
double Em[6][4];

double RelatePow(const Vec<double>& Vx, const int ex)
{
    int n;
    n = Vx.length();
    double ReSum = 0;
    for (int i = 0; i < n; i++)
    {
        ReSum += pow(Vx[i], (double)ex);
    }
    return ReSum;
}

double RelateMutiXY(const Vec<double>& Vx, const Vec<double>& Vy, const int ex)
{
    double dReMultiSum = 0;
    int n; 
    n = Vx.length();
    if(Vy.length()!=n)
    {
        cerr<<"X and Y dimension dismatch! Will return 0"<<endl;
        return 0;
    }
    for (long i = 0; i < n; i++)
    {
        dReMultiSum += pow(Vx[i], (double)ex)*Vy[i];
    }
    return dReMultiSum;
}

void EMatrix(const Vec<double>& Vx, const Vec<double>& Vy, int ex, Vec<double>& coefficient)
{
    int n = Vx.length();
    if(Vy.length()!=n)
    {
        cerr<<"X and Y dimension dismatch! Will truncate input."<<endl;
        n = Vx.length()<Vy.length() ? Vx.length():Vy.length();
    }
    for (int i = 1; i <= ex; i++)
    {
        for (int j = 1; j <= ex; j++)
        {
            Em[i][j] = RelatePow(Vx, i + j - 2);
        }
        Em[i][ex + 1] = RelateMutiXY(Vx, Vy, i - 1);
    }
    Em[1][1] = n;
    CalEquation(ex, coefficient);
}

void CalEquation(int exp, Vec<double>& coefficient)
{
    for (int k = 1; k < exp; k++) //消元过程
    {
        for (int i = k + 1; i < exp + 1; i++)
        {
            double p1 = 0;

            if (Em[k][k] != 0)
                p1 = Em[i][k] / Em[k][k];

            for (int j = k; j < exp + 2; j++)
                Em[i][j] = Em[i][j] - Em[k][j] * p1;
        }
    }
    coefficient[exp] = Em[exp][exp + 1] / Em[exp][exp];
    for (int l = exp - 1; l >= 1; l--)   //回代求解
        coefficient[l] = (Em[l][exp + 1] - F(coefficient, l + 1, (int)exp)) / Em[l][l];
}

double F(const Vec<double>& c, int l, int m)
{
    double sum = 0;
    for (int i = l; i <= m; i++)
        sum += Em[l - 1][i] * c[i];
    return sum;
}
