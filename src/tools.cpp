#include "tools.h"

NTL_CLIENT
using namespace std;
using namespace fplll;

void debug(){}


void genSalt(char *salt)
{
	struct timespec timeseed;
	clock_gettime(CLOCK_MONOTONIC, &timeseed);
	sprintf(salt, "%ld%ld%02d", timeseed.tv_sec, timeseed.tv_nsec,rand()%100);
}


double timeDiff(timespec start, timespec end)
{
	double t = 0;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		t = end.tv_sec-start.tv_sec-1;
		t += 1.0 + (double)(end.tv_nsec-start.tv_nsec)/1e9;
	} else {
		t = end.tv_sec-start.tv_sec;
		t += (double)(end.tv_nsec-start.tv_nsec)/1e9;
	}
	return t;
}

//  =========================== generate unimodular matrix ==========================
double Box_Muller()
{ //standard Gaussian distribution
	double x1, x2, y;
	x1 = rand()%RAND_MAX / (double)RAND_MAX;
	x2 = rand()%RAND_MAX / (double)RAND_MAX;
	y = sqrt(-2*log(x1)) * cos(2*pi*x2);
	return y;
}

void Gen_SparseTriangle(mat_ZZ &A, long m)
{	//VERY sparse! upper triangle
	long i, j;
	A.SetDims(m, m);
	for(i=0; i<m; i++){
		set(A[i][i]);
	}
	for(i=0; i<m; i++){
		for(j=i+1; j<m; j++){
			// A[i][j] = (int)round(Box_Muller()); // Gaussian in [-MAT_MAX, MAT_MAX]
			A[i][j] = (int)(Box_Muller());
		}
	}
}


// void Gen_Sparse_Triangle(mat_ZZ &A, long m)
// {	//upper triangular matrix A_[m,m] with about m positions setting to 1
// 	long i, j;
// 	for(i=0; i<m; i++){
// 		set(A[i][i]);
// 	}
// 	for(i=0; i<m; i++){
// 		for(j=i+1; j<m; j++){
// 			A[i][j] = 1; 
// 		}
// 	}
// }


void rePermutation(mat_ZZ &A)
{	//repermutation Ai
	long m = A.NumRows();
	long t = m<64 ? 8 : (m/8);
	t = m>160 ? 20 : t;
	long a, b;
	vec_ZZ temp;
	for(long i=0; i<t; i++)
	{
		a = rand()%m;
		b = rand()%m;
		temp = A[a];
		A[a] = A[b];
		A[b] = temp;
	}
}

void rePermutation_keephead(mat_ZZ &A, long head)
{	//repermutation Ai, avoiding changing the top [1:head] vectors
	long m = A.NumRows();
	long t = m<64 ? 8 : (m/8); 
	t = m>160 ? 20 : t;
	long a, b;
	vec_ZZ temp;
	for(long i=0; i<t; i++)
	{
		a = rand() % (m - head) + head;
		b = rand() % (m - head) + head;
		temp = A[a];
		A[a] = A[b];
		A[b] = temp;
	}
}

void Gen_Triangle(mat_ZZ &A, long m)
{	//upper triangle
	long i, j;
	A.SetDims(m, m);
	for(i=0; i<m; i++){
		set(A[i][i]);
	}
	for(i=0; i<m; i++){
		for(j=i+1; j<m; j++){
			A[i][j] = rand() % (2 * m + 1) - m; //ZZ <- int  uniformly in [-m, m]
		}
	}
}

void Gen_UniMat(mat_ZZ &A, long m)
{
	long i, j;
	mat_ZZ L, U;
	A.SetDims(m, m);
	Gen_Triangle(U, m);
	Gen_Triangle(L, m);
	L = transpose(L);//lower 
	A = L * U;
}

// ======================== calculate lattice information =======================
void generate_random_HNF(vec_ZZ& out,long n,long bit, ZZ seed)
{
    SetSeed(seed);
    out.SetLength(n); clear(out);
    ZZ p; GenPrime(p,bit*n);
    out(1) = p;
    for (int i=2; i<=n; i++)
    {
	RandomBnd(out(i),p);
    }
}

void GenerateBasis(mat_ZZ &B, long n, long seed)
	// WARNING : the SVP challenge offcial generator uses an NTL older than NTL 9.4
	 // (since NTL 9.4 and later versions use a different pseudorandom number generator) 
	// to create challenges on your local machine
{
	long bit  = 10; //do not change this
	vec_ZZ v; 
	generate_random_HNF(v, n, bit, conv<ZZ>(seed));
    B.SetDims(n,n);
    clear(B);
    B(1,1) = v(1);
    for (int i=2; i<=n; i++)
    {
		B(i,1)=v(i);
		B(i,i)=1;
    }
}

void G_S_Orthogonal(const mat_ZZ& B, mat_RR &B1, mat_RR &U)
//B1 is orthogonal basis, U =μ(i,j)
{
	long i,j;
	RR t3, t4;
	vec_RR t;
	long n;
	n = B.NumRows();
	t.SetLength(n);
	U.SetDims(n,n);
	B1.SetDims(n,n);
	for(i=1; i<=n; i++)
	{
		B1(i)=conv<vec_RR>(B(i));
		for(j=1; j<i; j++)
		{
			InnerProduct(t3, conv<vec_RR>(B(i)), B1(j));
			InnerProduct(t4, B1(j), B1(j));
			U(i,j) = t3/t4;
			//B1(i) -= U(i,j) * B1(j);
			mul(t, U(i,j), B1(j));
			sub(B1(i), B1(i), t);
		}
	}
	for(i=1;i<=n;i++)
		U(i,i)=1;
}

ZZ Factorial(long k)
 // return k!
{
	ZZ one;
	set(one);
	if(k == 0) 
		return one;
	else 
		return k * Factorial(k-1);
}

ZZ Factorial(const ZZ& k)
 // return k!
{
	ZZ one;
	set(one);
	if(k == 0) 
		return one;
	else 
		return k * Factorial(k-1);
}

ZZ Binomial(const ZZ& m, const ZZ& n) //C_m ^n //choose n in m
{
	ZZ res;
	res  = 1;
	ZZ i;
	for(i=m; i>m-n; i--)
		res *= i;
	res = res/Factorial(n);
	return res;
}

RR ballVol(long n, const RR& radius)
// volume of n-dimensional sphere with radius R
{
	RR PI = conv<RR>(pi);
	RR res, Cn;
	ZZ tmp;
	clear(res);
	res = conv<RR>(1);
	long k;
	if( n % 2 == 1 )
	//n = 2k+1 , C_n = C_(2k+1) = 2^{2k+1}k!Pi^k / (2k+1)!
	{ 
		k = (n-1)/2;
		Cn = power(conv<RR>(2), n) * conv<RR>(Factorial(k)) * power(PI, k);
		Cn = Cn / conv<RR>(Factorial(n));
	}
	else 
	//n =2k , Cn = pi^k / k!
	{
		k = n/2;
		Cn = power(PI, k) / conv<RR>(Factorial(k));
	}
	res = Cn * power(radius, n);
	return res;
}

RR fullEnumCost(long n, const RR& radius, const vec_RR& B1_norm2)
// number of nodes in Enum
{
	long i,k;
	RR N, Ni;
	clear(N);
	for(k=1; k<=n; k++)
	{
		RR prod; set(prod);
		// calculate denominator
		for(i=n-k+1; i<=n; i++)
		{
			mul(prod, prod, SqrRoot(B1_norm2[i-1]));
		}
		Ni = ballVol(k, radius) / prod;
		N += Ni;
	}
	return (N / conv<RR>(2));
}

RR gaussianHeuristic(const RR& R, const mat_ZZ& B)
// number of lattice points in a ball with radius R ~= B_n(R)/vol(L) 
{
	ZZ vol_L;
	RR vol_Ball;
	vol_L = SqrRoot(determinant( B*transpose(B) ));
	vol_Ball = ballVol(B.NumCols(), R);
	return vol_Ball/conv<RR>(vol_L);
}

RR GHL(const mat_ZZ& B)
{
	ZZ vol_L;
	RR vol_Ball;

	long m, n;
	m = B.NumCols();
	n = B.NumRows();

	vol_L = SqrRoot(determinant( B*transpose(B) ));
	vol_Ball = ballVol(n, conv<RR>(1));
	return pow(conv<RR>(vol_L)/vol_Ball, conv<RR>(1)/conv<RR>(n) );
}


// ====================== some random sampling functions =====================
void randomVector(vec_ZZ& v, long n, const ZZ& bound)
{ //n-dim vector v with range [0, bound)
	SetSeed(conv<ZZ>(rand()));
	v.SetLength(n);
	for(long i=0; i<n; i++)
	{
		RandomBnd(v[i], bound);
	}
}

void randomVector(vec_ZZ& v, long n, int bound)
{ //n-dim vector v with range [0, bound)
	v.SetLength(n);
	for(long i=0; i<n; i++)
	{
		v[i] = conv<ZZ>(rand()%bound);
	}
}

void randomVector(Vec<long>& v, long n, int bound)
{ //n-dim vector v with range [0, bound)
	v.SetLength(n);
	for(long i=0; i<n; i++)
	{
		v[i] = rand()%bound;
	}
}

void randomVector(Vec<int>& v, long n, int bound)
{ //n-dim vector v with range [0, bound)
	v.SetLength(n);
	for(long i=0; i<n; i++)
	{
		v[i] = rand()%bound;
	}
}

