#pragma once

#include <complex>
#include <cmath>
#include <math.h>
#include <assert.h>
#include "tools.h"

NTL_CLIENT
using namespace std;


const double PI = 3.1415926535897932384626433832795028841971693993751058209749446;
const double PI2 = 9.86960440108935861883449099987615113531369940724079062641335; //PI^2
const double PIroot = 1.77245385090551602729816748334114518279754945612238712821380779; //sqrt(PI)


void makeListOdd(Vec<Vec<int>> & List, const long & n);
void makeListEven(Vec<Vec<int>> & List, const long & n);


double P_succ_wrap(Vec<Vec<int>> & List, const long & n, const Vec<double> & B1norm2_d, \
    const double & target_d, const long & M, const int sound = 0);

//  ======== quick estimation =======
void findMaxCell(Vec<int>& cell_max, const Vec<Vec<int>>& List, const long & n);

void ERF_buildTable(Vec<Vec<Vec<complex<double>> > >& ERFtable, const Vec<double>& GSnorm, const double& R, \
    const Vec<int>& max_t, \
    const double& a, const long & max_m, const long& n);

double DPE_ImF_sm(const Vec<Vec<Vec<complex<double>> > >& ERFtable, const Vec<double>& GSnorm, \
    const double& R, const Vec<int>& t, const double& a, const long& m, const long& n);

double DPE_psucc(Vec<Vec<Vec<complex<double>> > >& ERFtable, const Vec<double>& GSnorm, \
    const double& R, const Vec<int>& t, const double& a, const long& n, \
    const long& K, const long& J);

double estimatePsuccEnsemble(Vec<Vec<Vec<complex<double>> > >& ERFtable, const Vec<Vec<int>>& List, \
    const Vec<double>& GSnorm, const Vec<double>& GSnorm2_d,\
    const long& M, const double& R, const double& a, const long& n, const long& K, const long& J);


// ==== concrete psucc calculation ======
double DPE_ImF_sm_sound(const Vec<double>& GSnorm, const double& R, const Vec<int>& t, \
    const double& a, const long& m, const long& n);


double DPE_psucc_sound(const Vec<double>& GSnorm, const double& R, const Vec<int>& t, \
    const double& a, const long& n, const long& K, const long& J);


double estimatePsuccEnsemble_sound(const Vec<Vec<int>>& List, const Vec<double>& GSnorm, \
    const Vec<double>& GSnorm2_d,const long& M, const double& R, \
    const double& a, const long& n, const long& K, const long& J);

//  =========== functional tools ================
complex<double> EQ1_cplx_Fm(const double& a, const long& m, const long& n);

complex<double> sqrtz_sm(double a, long m);

void quickSort(Vec<Vec<int>>& List, Vec<double>& Ev, long low, long high);

inline complex<double> mul_db_cplx(const double& a, const complex<double>& z)
// complex * double
{
	return complex<double>(a*z.real(), a*z.imag());
}

inline complex<double> div_cplx_db(const complex<double>& z, const double& a)
// complex / double
{
	return complex<double>(z.real()/a, z.imag()/a);
}

inline double arg(const complex<double>& z)
// this is exactly the same with atan2() function in <math.h>
{
	if(z.real()>0)
		return atan(z.imag()/z.real());
	if(z.real()==0 && z.imag()>0)
		return PI/2;
	if(z.real()<0 && z.imag()>=0)
		return PI + atan(z.imag()/z.real());
	if(z.real()<0 && z.imag()<0)
		return atan(z.imag()/z.real()) - PI;
	if(z.real()==0 && z.imag()<0)
		return -1*PI/2;
}


inline double sinc(double x)
{
    if (fabs(x) < 1E-10) return 1.;
    else return sin(x) / x;
}

inline double exp_minus(double a, double b)
{
    assert(b >= 0);
    if (b > 0.5) {
        return exp(a + b) - exp(a - b);
    } else {
        return exp(a - b) * expm1(2. * b);
    }
}

// This implementation is following:
// https://granite.phys.s.u-tokyo.ac.jp/svn/LCGT/trunk/sensitivity/Matlab/official/@double/erfz.pdf
// implementation: lihawl AT hotmail DOT com
inline complex<double> erfz(std::complex<double> const &c)
{
    // OK! I think one can stil improve the speed of the code by at least 10%, by some naive optimization.
    // But before you do some `smart` optimizations, think about if your optimization will cause an overflow!
    // And after any modification, run the test_erfz.cpp to check results!
    typedef complex<double> Complex;
    double const Pi = 3.141592653589793238462643383279502884197169399375;
    double const Re = c.real();
    double const R2 = Re * Re;

    double Im_ = c.imag();
    if (Im_ == 0) return complex<double>(erf(Re), 0);
    bool const conj = Im_ < 0;
    double const Im = conj ? (-Im_) : Im_;

    double const erfR = erf(Re);


    Complex const E_ = exp(-R2) / (2 * Pi) * 2 * Im *
        Complex(sin(Re * Im) * sinc(Re * Im), sinc(2 * Re * Im));


    // Note: exp(-0.25 * 12 ^ 2) / (0.25 * 12 ^ 2) = 6E-18
    int const N = 15;
    double F_R = 0;
    for (int n = N; n >= 1; --n) {
        F_R += exp(-0.25 * n * n) / (0.25 * n * n + R2);
    }
    F_R *= exp(-R2) * Re / Pi;


    int const M = (int)(2 * Im);
    int const M_N = max(M - N, 1);

    Complex HG(0, 0);
    Complex H(0, 0);
    Complex G(0, 0);

    // for H
    for (int n = min(M_N - 1, N); n >= 1; --n) {
        H += exp(-0.25 * n * n - n * Im - R2) / (0.25 * n * n + R2) * Complex(Re, 0.5 * n);
    }

    // overlap of H & G
    for (int n = N; n > M_N - 1; --n) {

        double HG_R = (exp(-0.25 * n * n - n * Im - R2) + exp(-0.25 * n * n + n * Im - R2)) / (0.25 * n * n + R2) * Re;
        double HG_I = exp_minus(-0.25 * n * n - R2, n * Im) / (0.25 * n * n + R2) * (-0.5 * n);

        HG.real(HG.real() + HG_R);
        HG.imag(HG.imag() + HG_I);
    }

    // for G
    for (int n = M + N; n > max(M_N - 1, N); --n) {
        G += exp(-0.25 * n * n + n * Im - R2) / (0.25 * n * n + R2) * Complex(Re, -0.5 * n);
    }

    H *= 1. / (2 * Pi);
    G *= 1. / (2 * Pi);
    HG *= 1. / (2 * Pi);

    Complex a(cos(-2 * Re * Im), sin(-2 * Re * Im));
    Complex b = a * (H + G + HG);
    double real = erfR + E_.real() + F_R - b.real();
    double imag = 0 + E_.imag() + 0 - b.imag();

    if (conj) return Complex(real, -imag);
    return Complex(real, imag);
}