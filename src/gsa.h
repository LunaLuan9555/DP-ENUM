#ifndef GSA__H
#define GSA__H

#include <complex>
#include <cmath>
#include <math.h>
#include <assert.h>
#include <vector>
#include "tools.h"

NTL_CLIENT
using namespace std;


void makeListEven(Vec<Vec<int>> & List, const long & n);

// ============ GSA ================
double SimulateCost_q(const ZZ& det, const int bkz, const long n, const double target_d, const long M);

double SimulateCost_GSA(const Vec<double> & B1norm2_d, const double target_d, const long M);

double SimulateCost_q_even(const ZZ& det, const int bkz, const long n, const double target_d, const long M);

double SimulateCost_GSA_even(const Vec<double> & B1norm2_d, const double target_d, const long M);

void GSA_fit(Vec<double> & GS_fit, const Vec<double> & GS2);

// ==================== linear regression(least squares method) =====================

double RelatePow(const Vec<double>& Vx, const int ex);

double RelateMutiXY(const Vec<double>& Vx, const Vec<double>& Vy, const int ex);

void EMatrix(const Vec<double>& Vx, const Vec<double>& Vy, int ex, Vec<double>& coefficient);

void CalEquation(int exp, Vec<double>& coefficient);

double F(const Vec<double>& c, int l, int m); 

// #define N 1e-13


#endif