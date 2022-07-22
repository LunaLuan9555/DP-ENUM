#ifndef TOOLS__H
#define TOOLS__H

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <iostream>

#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/RR.h>
#include <NTL/vec_RR.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>

#include <fplll.h>
#include <fplll/bkz.h>

#define MPI_BUFFERSIZE 524288 //512kb
#define PRUN_PATH "default.json"

NTL_CLIENT

using namespace std;
using namespace fplll;

// -------------- read arguments --------------------

template<typename T>
std::istream& operator>>(const std::string& in, T& v)
{
    static std::istringstream* iss = 0;
    if (iss) delete iss;
    iss = new std::istringstream(in);
    *iss>>v;
    return *iss;
}

#define PARSE_MAIN_ARGS \
    std::ostringstream message; \
    message << argv[0];  \
    for (int i=0; i<argc; ++i)

#define MATCH_MAIN_ARGID(arg,dest) \
    if (i==0)  {\
	message << " [" << arg << " " << dest << "]"; \
    } else { \
	if (std::string(argv[i]) == std::string(arg))  {\
	    std::string(argv[i+1]) >> dest; \
	    i++; continue; \
	} \
    }

#define DETECT_MAIN_ARGID(arg,dest,value) \
    if (i==0)  {\
	message << "[" << arg << "]"; \
    } else { \
	if (std::string(argv[i]) == std::string(arg))  {\
	    dest = value; \
	    continue; \
	} \
    }

#define SYNTAX() \
    if (i==0) continue; \
    std::cerr << "Syntax: " << message.str() << std::endl; \
    exit(1);


// ================== gdb debug and tools in main.cpp ==================
void debug();

double timeDiff(timespec start, timespec end);

void genSalt(char *salt);

// ------------------------mpi communication--------------------------

template <class Ty>
void conv_str(string &result, Ty A){
	ostringstream oss;
	oss.clear();
	oss<<A;
	result = oss.str();
}

template <class Ty>
void conv_NTLtype(Ty &A, string mystr){
	stringstream ss;
	ss.clear();
	ss.str(mystr);
	ss >> A;
}

template <class Ty>
void conv_fptype(Ty &A, string mystr){
	stringstream ss;
	ss.clear();
	ss.str(mystr);
	ss >> A;
}

template <class Ta, class Tb>
void convtpye_A2B(const Ta& A, Tb & B){
	ostringstream sA;
	sA.clear();
	stringstream sB;
	sB.clear();
	string tmp;
	sA << A;
	tmp = sA.str();
	sB.str(tmp);
	sB >> B;
}

// ------------------------ tools --------------------------
template <class T>
bool isZero(const Vec<T>& v){
	long n = v.length();
	for(long i = 0; i < n; i++){
		if(v[i])
		{
			return false;
		}
	}
	return true;
}

template <class T>
void clear(Vec<T>& v)
{
	long n = v.length();
	for(long i=0; i<n; i++)
		v[i] = 0;
}

// ------------------------ arithmetic --------------------------
// transpose matrix
template <class T>
void TRANSPOSE(Mat<T>& Bt, const Mat<T>& B)
{
	long n = B.NumRows();
	long m = B.NumCols();
	Bt.kill();
	Bt.SetDims(m, n);
	for(long i = 0; i<m; i++)
		for(long j = 0; j<n; j++)
			Bt[i, j] = B[j, i];
}

template <class T>
Mat<T> TRANSPOSE(const Mat<T>& B)
{
	Mat<T> Bt;
	long n = B.NumRows();
	long m = B.NumCols();
	Bt.SetDims(m, n);
	for(long i = 0; i<m; i++)
		for(long j = 0; j<n; j++)
			Bt[i, j] = B[j, i];
	return Bt;
}

template <class T>
T SUM(const Vec<T>& v)
{
	long n = v.length();
	T sum;
	sum = 0;
	for(long i=0; i<n; i++)
		sum += v[i];
	return sum;
}

// c = a + b
template <class T>
void ADD(Vec<T>& c, const Vec<T>& a, const Vec<T>& b)
{
	long n = a.length();
	if(b.length() != n)
	{
		cerr<<"ERROR: dimension dismatch!"<<endl;
		c.SetLength(0);
		return; 
	}
	c.SetLength(n);
	clear(c);
	for(long i=0; i<n; i++)
		c[i] = a[i] + b[i];
}

// c = a - b
template <class T>
void SUB(Vec<T>& c, const Vec<T>& a, const Vec<T>& b)
{
	long n = a.length();
	if(b.length() != n)
	{
		cerr<<"ERROR: dimension dismatch!"<<endl;
		c.SetLength(0);
		return; 
	}
	c.SetLength(n);
	clear(c);
	for(long i=0; i<n; i++)
		c[i] = a[i] - b[i];
}

// vector = vector * matrix -> c=a*B
template <class T>
void MUL(Vec<T>& c, const Vec<T>& a, const Mat<T>& B)
{
	long n = B.NumRows();
	long m = B.NumCols();
	if(a.length() != n)
	{
		cerr<<"ERROR: dimension dismatch!"<<endl;
		c.SetLength(0);
		return; 
	}
	c.SetLength(m);
	clear(c);
	for(long i=0; i<m; i++)
	{
		for(long j = 0; j<n; j++)
		{
			c[i] = c[i] + a[j] * B[j][i];
		}
	}
}

// innerproduct
template <class T>
void MUL(T& c, const Vec<T>& a, const Vec<T>& b)
{
	long n = a.length();
	c = 0;
	if(b.length() != n)
	{
		cerr<<"ERROR: dimension dismatch!"<<endl;
		
		return; 
	}
	for(long i=0; i<n; i++)
	{
		c += a[i] * b[i];
	}
}

template <class T>
void MUL(ZZ& c, const Vec<T>& a, const Vec<T>& b)
{
	long n = a.length();
	clear(c);
	if(b.length() != n)
	{
		cerr<<"ERROR: dimension dismatch!"<<endl;
		
		return; 
	}
	for(long i=0; i<n; i++)
	{
		c = c + conv<ZZ>(a[i] * b[i]);
	}
}

template <class T>
void SWAP(T& a, T& b)
{
	T tmp;
	tmp = a;
	a = b;
	b = tmp;
}
// ------------------------UniMat.cpp--------------------------
const double pi = 3.1415926897932384;

double Box_Muller();

void Gen_SparseTriangle(mat_ZZ &A, long m);

void Gen_Triangle(mat_ZZ &A, long m);

void Gen_UniMat(mat_ZZ &M, long m);

void rePermutation(mat_ZZ &A);

void rePermutation_keephead(mat_ZZ &A, long head);




// ------------------------ latticeINFO.cpp --------------------------

void generate_random_HNF(vec_ZZ& out,long n,long bit, ZZ seed);

void GenerateBasis(mat_ZZ &B, long n, long seed=0);

void G_S_Orthogonal(const mat_ZZ& B, mat_RR &B1, mat_RR &U);

ZZ Factorial(long k);

ZZ Factorial(const ZZ& k);

ZZ Binomial(const ZZ& m, const ZZ& n);

RR ballVol(long n, const RR& radius);

RR fullEnumCost(long n, const RR& radius, const vec_RR& B1_norm2);

RR gaussianHeuristic(const RR& R, const mat_ZZ& B);

RR GHL(const mat_ZZ& B);


// ------------------------ sampling.cpp --------------------------
void randomVector(vec_ZZ& v, long n, const ZZ& bound);

void randomVector(vec_ZZ& v, long n, int bound);

void randomVector(Vec<long>& v, long n, int bound);

void randomVector(Vec<int>& v, long n, int bound);



#endif

