#include "tools.h"
#include "cellenum.h"
#include "reduction.h"

// define candidate vector (squared) length
#define RADIUS_RELAXATION_FACTOR2	1.2
#define INSERTION_FACTOR_DELTA		0.99
#define BIN_THRESHOLD 0.005

NTL_CLIENT
using namespace std;
using namespace fplll;

bool lastOdd(const Vec<int>& vout, long n)
{
	long k = n-1;
	while((!vout[k]) && k) // vout[k]==0 AND k>0
	{
		k--;
	}
	if(vout[k]%2) // last non-zero component of vout[] is odd
		return true;
	else // even
		return false;
}

void updateBasis(mat_ZZ &B, const cnode& node, long n)
{
	vec_ZZ v;
	mat_ZZ Bnew;
	long i, j, h_v;
	h_v = node.h_v;
	conv(v, node.vec_v);

	Bnew.SetDims(n+1,n);
	Bnew(h_v) = v;
	for(i=1; i<=n+1; i++)
	{
		if(i < h_v)
			Bnew(i) = B(i);
		else if(i > h_v)
			Bnew(i) = B(i-1);
	}
	//reduce basis B
	LLL_default(Bnew, 0.99, FT_DEFAULT, 100);
	j=1;
	for(i=1; i<=n+1; i++)
	{
		if( !IsZero(Bnew(i)) )
		{
			B(j) = Bnew(i);
			j++;
		}
		else continue;
	}
}

// === decoding with insertion ===
void ENUM_discrete_strategy(Vec<long>& sol, Vec<int>& sol_z, Vec<cnode>& Candidate, const Vec<Vec<int>> List, const Mat<long>& B, const Vec<double>& B1norm2_d, \
	const Mat<double>& U, const double& target2, const long & n)
{
	long N = List.length();
	Candidate.SetLength(0);
	Vec<int> z;
	Vec<long> sol_temp;
	int i, j, d;
	double t;
	double y;
	Vec<long> u;
	double partial = 0;
	u.SetLength(n);
	int findflag = 1; //TRUE
	long h_v;

	double target2_relax = target2 * RADIUS_RELAXATION_FACTOR2;
	Vec<double> B1norm2_d_relax;
	B1norm2_d_relax = B1norm2_d; // 0.99 ||b_i^*||^2
	for(i=0; i<n; i++)
	{
		B1norm2_d_relax[i] = B1norm2_d[i] * INSERTION_FACTOR_DELTA;
	}

	for(int k = 0; k < N; k++)
	{
		z = List[k];
		partial = 0;
		h_v = n+1; 
		findflag = true;
		for(i = 0; i < n; i++)
		{
			u[i] = 0;
		}
		for(i=n-1; i>=0; i--) // μ(i,j) u[i] is the cordinate for v=\sum {u_i \dot b_i}
		{
			y = 0;
			for(j=i+1; j<n; j++) // negative partial sum:  - \sum_{j=i+1}^{n}{u_j \dot b_j^*}
			{
				y  -=  u[j] * U[j][i];
			}
			u[i] = floor(y + 0.5);
			// 
			if(u[i] <= y)
				d = -1; // (-1)^zi = -1
			else
				d = 1;
			// 
			long tmp; 
			if(z[i]%2 ==0)
			{
				tmp = z[i]/2;
				u[i] = u[i] + d * tmp;
			}
			else
			{
				tmp = (z[i] + 1) / 2;
				u[i] = u[i] - d * tmp;
			}
			// after determin u[i], calculating partial sum of ||v||^2
			partial += (u[i] - y) * (u[i] - y) * B1norm2_d[i];
			if(partial > target2_relax)
			{
				findflag = 0;
				break; //this z is longer than target length
			}
			else if (partial > target2) // candidate vector, calculating h(v)?
			{
				findflag = -1;
				if(partial < B1norm2_d_relax[i] ) // ||pi_i(v)||^2 < 0.99 ||b_i^*||^2, update h(v)
					h_v = i+1;
			}
			// else: partial<target2, continue;
		}
		
		// deal with coeff u[] and h_v;
		if(findflag == 1)
		{
			sol_z = z;
			MUL(sol, u, B); //v = u * B
			break;
		}
		else if (findflag == -1)
		{
			cnode node_temp;
			MUL(node_temp.vec_v, u, B); //v = u * B
			node_temp.vec_z = z;
			node_temp.h_v = h_v;
			node_temp.v_norm = (long)partial; // not accuracy, just an estimation
			Candidate.append(node_temp);
		}
	}
}

// === cellENUM + decoding with insertion ===
int DP_ENUM_ins_decode(Vec<long>& sol, Vec<int>& sol_z, Vec<cnode>& Candidate, const Vec<int>& z, const Mat<long>& B, \
	const Vec<double>& B1norm2_d, const Vec<double>& B1norm2_d_relax,\
	const Mat<double>& U, const double & target2, const double & target2_relax, const long & n)
{
	int i, j, d;
	double t;
	double y;
	Vec<long> u;
	double partial;
	u.SetLength(n);
	int findflag = 1;
	// === decode ===
	partial = 0;
	long h_v = n+1; 
	findflag = 1;
	for(i = 0; i < n; i++)
	{
		u[i] = 0;
	}
	for(i=n-1; i>=0; i--) //μ(i,j) u[i] is the cordinate for v=\sum {u_i \dot b_i}
	{
		y = 0;
		for(j=i+1; j<n; j++) // negative partial sum:  - \sum_{j=i+1}^{n}{u_j \dot b_j^*}
		{
			y  -=  u[j] * U[j][i];
		}
		u[i] = floor(y + 0.5);
		// 
		if(u[i] <= y)
			d = -1; //(-1)^zi = -1
		else
			d = 1;
		// 
		long tmp; 
		if(z[i]%2 ==0)
		{
			tmp = z[i]/2;
			u[i] = u[i] + d * tmp;
		}
		else
		{
			tmp = (z[i] + 1) / 2;
			u[i] = u[i] - d * tmp;
		}
		// after determin u[i], calculating partial sum of ||v||^2
		partial += (u[i] - y) * (u[i] - y) * B1norm2_d[i];
		if(partial > target2_relax)
		{
			findflag = 0;
			break; //this z is longer than target length
		}
		else if (partial > target2) // candidate vector, calculating h(v)?
		{
			findflag = -1;
			if(partial < B1norm2_d_relax[i] ) // ||pi_i(v)||^2 < 0.99 ||b_i^*||^2, update h(v)
				h_v = i+1;
		}
		// else: partial<target2, continue;
	}
	if(findflag == 1)
	{
		sol_z = z;
		MUL(sol, u, B); //v = u * B
		return 1;
	}
	else if(findflag == -1)
	{
		cnode node_temp;
		MUL(node_temp.vec_v, u, B); //v = u * B
		node_temp.vec_z = z;
		node_temp.h_v = h_v;
		node_temp.v_norm = (long)partial; // not accuracy, just an estimation
		Candidate.append(node_temp);
		return -1;
	}
	else
		return 0;
}

// === cellENUM + decoding with insertion ===
void DP_ENUM_ins_Wrap(long & counter, const Vec<int>& pi, double Bound, const Vec<double>& r, long n, \
	Vec<long>& sol, Vec<int>& sol_z, Vec<cnode>& Candidate, const Mat<long>& B, const Vec<double>& B1norm2_d, \
	const Mat<double>& U, const double& target2)
{
	Candidate.SetLength(0);
	counter = 0;
	long i, j, k;
	int findsol;
	// 
	Vec<int> v, vout; // current node, output pi(v[i])
	Vec<double> c; // partial length
	v.SetLength(n);
	clear(v); 
	// v[0] = 1;
	vout.SetLength(n);
	clear(vout); 
	c.SetLength(n+1);
	clear(c);

	double target2_relax = target2 * RADIUS_RELAXATION_FACTOR2;
	Vec<double> B1norm2_d_relax;
	B1norm2_d_relax = B1norm2_d; // 0.99 ||b_i^*||^2
	for(i=0; i<n; i++)
	{
		B1norm2_d_relax[i] = B1norm2_d[i] * INSERTION_FACTOR_DELTA;
	}

	// 
	for(k=n-1; k>=1; k--)
	{
		c[k] = c[k+1] + (v[k]+1)*(v[k]) * r[k]; // using f(t) defined in [ANS18]
	}
	// 
	k = 0;
	while(1)
	{
		c[k] = c[k+1] + (v[k]+1)*(v[k]) * r[k];
		if(c[k] < Bound)
		{
			if(k==0) //find a cell
			{
				for(i=0; i<n; i++)
				{
					vout[pi[i]] = v[i];
				}
				if(lastOdd(vout, n)) // the last non-zero component of vout[] must be odd
				{
					// GOTO decoding 
					counter++;
					findsol = DP_ENUM_ins_decode(sol, sol_z, Candidate, vout, B, B1norm2_d, B1norm2_d_relax, U, \
						target2, target2_relax, n);
					if(findsol==1){
						break;
					}
					else if(findsol == -1)
					{
						// append to candidate list
					}
				}
				v[k] = v[k] + 1; 
			}
			else //go down the tree
			{
				k--;
				v[k] = 0;
			}
		}
		else
		{
			k++;
			if(k == n)
				break; // no more vectors
			//else k < n
			v[k] = v[k] + 1;
		}
	}
}


// === decoding with out insertion ===
void ENUM_discrete_postprocessing(Vec<long>& sol, Vec<int>& sol_z, const Vec<Vec<int>> List, const Mat<long>& B, \
	const Vec<double>& B1norm2_d, const Mat<double>& U, const double& target2, const long & n)
{
	long N = List.length();
	Vec<int> z;
	int i, j, d;
	double t;
	double y;
	Vec<long> u;
	double partial = 0;
	u.SetLength(n);
	bool findflag = true;
	for(int k = 0; k < N; k++)
	{
		z = List[k];
		partial = 0;
		findflag = true;
		for(i = 0; i < n; i++)
		{
			u[i] = 0;
		}
		for(i=n-1; i>=0; i--) //μ(i,j) u[i] is the cordinate for v=\sum {u_i \dot b_i}
		{
			y = 0;
			for(j=i+1; j<n; j++) // negative partial sum:  - \sum_{j=i+1}^{n}{u_j \dot b_j^*}
			{
				y  -=  u[j] * U[j][i];
			}
			u[i] = floor(y + 0.5);
			// 
			if(u[i] <= y)
				d = -1; //(-1)^zi = -1
			else
				d = 1;
			// 
			long tmp; 
			if(z[i]%2 ==0)
			{
				tmp = z[i]/2;
				u[i] = u[i] + d * tmp;
			}
			else
			{
				tmp = (z[i] + 1) / 2;
				u[i] = u[i] - d * tmp;
			}
			// after determin u[i], calculating partial sum of ||v||^2
			partial += (u[i] - y) * (u[i] - y) * B1norm2_d[i];
			if(partial > target2)
			{
				findflag = false;
				break; //this z is longer than target length
			}
		}
		if(findflag)
		{
			sol_z = z;
			MUL(sol, u, B); //v = u * B
			break;
		}
	}
}

// === cell ENUM + decoding with out insertion ===
int DP_ENUM_decode(Vec<long>& sol, Vec<int>& sol_z, const Vec<int>& z, const Mat<long>& B, const Vec<double>& B1norm2_d, \
	const Mat<double>& U, const double& target2, const long & n)
{
	int i, j, d;
	double t;
	double y;
	Vec<long> u;
	double partial;
	u.SetLength(n);
	bool findflag = true;

	// === decode ===
	partial = 0;
	findflag = true;
	for(i = 0; i < n; i++)
	{
		u[i] = 0;
	}
	for(i=n-1; i>=0; i--) //μ(i,j) u[i] is the cordinate for v=\sum {u_i \dot b_i}
	{
		y = 0;
		for(j=i+1; j<n; j++) // negative partial sum:  - \sum_{j=i+1}^{n}{u_j \dot b_j^*}
		{
			y  -=  u[j] * U[j][i];
		}
		u[i] = floor(y + 0.5);
		// 
		if(u[i] <= y)
			d = -1; //(-1)^zi = -1
		else
			d = 1;
		// 
		long tmp; 
		if(z[i]%2 ==0)
		{
			tmp = z[i]/2;
			u[i] = u[i] + d * tmp;
		}
		else
		{
			tmp = (z[i] + 1) / 2;
			u[i] = u[i] - d * tmp;
		}
		// after determin u[i], calculating partial sum of ||v||^2
		partial += (u[i] - y) * (u[i] - y) * B1norm2_d[i];
		if(partial > target2)
		{
			findflag = false;
			break; //this z is longer than target length
		}
	}
	if(findflag)
	{
		sol_z = z;
		MUL(sol, u, B); //v = u * B
		return 1;
	}
	else
		return 0;
}

// === cell ENUM + decoding with out insertion ===
void DP_ENUM_Wrap(long & counter, const Vec<int>& pi, double Bound, const Vec<double>& r, long n, \
	Vec<long>& sol, Vec<int>& sol_z, const Mat<long>& B, const Vec<double>& B1norm2_d, \
	const Mat<double>& U, const double target2)
{
	counter = 0;
	long i, j, k;
	int findsol;
	// 
	Vec<int> v, vout; // current node, output pi(v[i])
	Vec<double> c; // partial length
	v.SetLength(n);
	clear(v); 
	// v[0] = 1;
	vout.SetLength(n);
	clear(vout); 
	c.SetLength(n+1);
	clear(c);
	// 
	for(k=n-1; k>=1; k--)
	{
		c[k] = c[k+1] + (v[k]+1)*(v[k]) * r[k]; // using f(t) defined in [ANS18]
	}
	// 
	k = 0;
	while(1)
	{
		c[k] = c[k+1] + (v[k]+1)*(v[k]) * r[k];
		if(c[k] < Bound)
		{
			if(k==0) //find a cell
			{
				for(i=0; i<n; i++)
				{
					vout[pi[i]] = v[i];
				}
				if(lastOdd(vout, n)) // the last non-zero component of vout[] must be odd
				{
					// GOTO decoding 
					counter++;
					findsol = DP_ENUM_decode(sol, sol_z, vout, B, B1norm2_d, U, target2, n);
					if(findsol==1){
						break;
					}
				}
				v[k] = v[k] + 1; 
			}
			else //go down the tree
			{
				k--;
				v[k] = 0;
			}
		}
		else
		{
			k++;
			if(k == n)
				break; // no more vectors
			//else k < n
			v[k] = v[k] + 1;
		}
	}
}

// === cell enumeration ===
void cellENUM(Vec<Vec<int>>& List, long & counter, const Vec<int>& pi, double Bound, const Vec<double>& r, long n)
{
	counter = 0;
	long i, j, k;
	// 
	Vec<int> v, vout; // current node, output pi(v[i])
	Vec<double> c; // partial length
	v.SetLength(n);
	clear(v); 
	// v[0] = 1;
	vout.SetLength(n);
	clear(vout); 
	c.SetLength(n+1);
	clear(c);
	// 
	for(k=n-1; k>=1; k--)
	{
		c[k] = c[k+1] + (v[k]+1)*(v[k]) * r[k]; // using f(t) defined in [ANS18]
	}
	// 
	k = 0;
	while(1)
	{
		c[k] = c[k+1] + (v[k]+1)*(v[k]) * r[k];
		if(c[k] < Bound)
		{
			if(k==0) //find a cell
			{
				for(i=0; i<n; i++)
				{
					vout[pi[i]] = v[i];
				}
				if(lastOdd(vout, n)) // the last non-zero component of vout[] must be odd
				{
					List.append(vout);
					counter++;
				}
				v[k] = v[k] + 1; 
			}
			else //go down the tree
			{
				k--;
				v[k] = 0;
			}
		}
		else
		{
			k++;
			if(k == n)
				break; // no more vectors
			//else k < n
			v[k] = v[k] + 1;
		}
	}
}

long cellENUM_probe(const Vec<int>& pi, double Bound, const Vec<double>& r, long n, long M)
{
	long counter = 0;
	long i, j, k;
	// 
	Vec<int> v, vout; // current node, output pi(v[i])
	Vec<double> c; // partial length
	v.SetLength(n);
	clear(v);
	// v[0] = 1;
	vout.SetLength(n);
	clear(vout);
	c.SetLength(n+1);
	clear(c);
	
	// === define bounds for counter ===
	// double THRESHOLD = 0.05;
	double rightM = M * (1 + BIN_THRESHOLD);
	double leftM = M * (1 - BIN_THRESHOLD);
	// 
	for(k=n-1; k>=0; k--)
	{
		// c[k] = c[k+1] + (v[k]+0.5)*(v[k]+0.5) * r[k];
		c[k] = c[k+1] + (v[k]+1)*(v[k]) * r[k]; // using f(t) defined in [ANS18]
	}
	// 
	k = 0;
	while(1)
	{
		c[k] = c[k+1] + (v[k]+1)*(v[k]) * r[k];
		if(c[k] < Bound)
		{
			if(k==0) //find a cell
			{
				for(i=0; i<n; i++)
				{
					vout[pi[i]] = v[i];
				}
				if(lastOdd(vout, n)) // the last non-zero component of vout[] must be odd
				{
					counter++;
					if(counter > rightM) //last non-zero: both even and odd are counted, counter/2
						return 1; // Bound too large
				}
				v[k] = v[k] + 1; 
			}
			else //go down the tree
			{
				k--;
				v[k] = 0;
			}
		}
		else
		{
			k++;
			if(k == n)
				break; // no more vectors
			//else k < n
			v[k] = v[k] + 1;
		}
	}
	if(counter > rightM) 
		return 1;
	else if(counter < leftM)
		return -1; // Bound too small
	else
		return 0; // Bound acceptable
}

double computeBound(double Bleft, double Bright, const Vec<int>& pi, const Vec<double>& r, long n, long M) //binary search Bound r
{
	double Bmid;
	long flag;
	while(Bleft < Bright)
	{
		Bmid = (Bleft + Bright)/2;
		flag = cellENUM_probe(pi, Bmid, r, n, M);
		if(flag > 0.5) // Bmid too large
		{
			Bright = Bmid;
		}
		else if(flag < -0.5) // Bmid too small
		{
			Bleft = Bmid;
		}
		else
		{
			return Bmid;
		}
	}
	return -1; //seraching failed!
}

void computePrunPara(Vec<int>& pi, Vec<double>& r, const Vec<double>& Bnorm2, long n)
// compute permutation and radius r[i]
{
	long i, j;
	pi.SetLength(n);
	r.SetLength(n);
	for(i=0; i<n; i++)
	{
		pi[i] = i;
		r[i] = Bnorm2[i];
	}
	for(i=0; i<n-1; i++)
	{
		for(j=0; j<n-1-i; j++)
		{
			if(r[j] > r[j+1])
			{
				swap(r[j], r[j+1]);
				swap(pi[j], pi[j+1]);//saving permutation
			}
		}
	}
}

void ENUM_discrete_strategy_M(long M, Vec<long>& sol, Vec<int>& sol_z, Vec<cnode>& Candidate, const Vec<Vec<int>> List, const Mat<long>& B, const Vec<double>& B1norm2_d, \
	const Mat<double>& U, const double target2, long n)
{
	long N = List.length();
	Candidate.SetLength(0);
	Vec<int> z;
	Vec<long> sol_temp;
	int i, j, d;
	double t;
	double y;
	Vec<long> u;
	double partial = 0;
	u.SetLength(n);
	int findflag = 1; //TRUE
	long h_v;

	double target2_relax = target2 * RADIUS_RELAXATION_FACTOR2;
	Vec<double> B1norm2_d_relax;
	B1norm2_d_relax = B1norm2_d; // 0.99 ||b_i^*||^2

	for(i=0; i<n; i++)
	{
		B1norm2_d_relax[i] = B1norm2_d[i] * INSERTION_FACTOR_DELTA;
	}
	N = M<N ? M:N;
	for(int k = 0; k < N; k++)
	{
		z = List[k];
		partial = 0;
		h_v = n+1; 
		findflag = true;
		for(i = 0; i < n; i++)
		{
			u[i] = 0;
		}
		for(i=n-1; i>=0; i--) // μ(i,j) u[i] is the cordinate for v=\sum {u_i \dot b_i}
		{
			y = 0;
			for(j=i+1; j<n; j++) // negative partial sum:  - \sum_{j=i+1}^{n}{u_j \dot b_j^*}
			{
				y  -=  u[j] * U[j][i];
			}
			u[i] = floor(y + 0.5);
			// 
			if(u[i] <= y)
				d = -1; // (-1)^zi = -1
			else
				d = 1;
			// 
			long tmp; 
			if(z[i]%2 ==0)
			{
				tmp = z[i]/2;
				u[i] = u[i] + d * tmp;
			}
			else
			{
				tmp = (z[i] + 1) / 2;
				u[i] = u[i] - d * tmp;
			}
			// after determin u[i], calculating partial sum of ||v||^2
			partial += (u[i] - y) * (u[i] - y) * B1norm2_d[i];
			if(partial > target2_relax)
			{
				findflag = 0;
				break; //this z is longer than target length
			}
			else if (partial > target2) // candidate vector, calculating h(v)?
			{
				findflag = -1;
				if(partial < B1norm2_d_relax[i] ) // ||pi_i(v)||^2 < 0.99 ||b_i^*||^2, update h(v)
					h_v = i+1;
			}
			// else: partial<target2, continue;
		}
		
		// deal with coeff u[] and h_v;
		if(findflag == 1)
		{
			sol_z = z;
			MUL(sol, u, B); //v = u * B
			break;
		}
		else if (findflag == -1)
		{
			cnode node_temp;
			MUL(node_temp.vec_v, u, B); //v = u * B
			node_temp.vec_z = z;
			node_temp.h_v = h_v;
			node_temp.v_norm = (long)partial; // not accuracy, just an estimation
			Candidate.append(node_temp);
		}
	}
}

void ENUM_discrete_postprocessing_M(long M, Vec<long>& sol, Vec<int>& sol_z, const Vec<Vec<int>> List, const Mat<long>& B, const Vec<double>& B1norm2_d, \
	const Mat<double>& U, const double target2, long n)
{
	long N = List.length();
	Vec<int> z;
	int i, j, d;
	double t;
	double y;
	Vec<long> u;
	double partial = 0;
	u.SetLength(n);
	bool findflag = true;
	N = M<N ? M:N;
	for(int k = 0; k < N; k++)
	{
		z = List[k];
		partial = 0;
		findflag = true;
		for(i = 0; i < n; i++)
		{
			u[i] = 0;
		}
		for(i=n-1; i>=0; i--) //μ(i,j) u[i] is the cordinate for v=\sum {u_i \dot b_i}
		{
			y = 0;
			for(j=i+1; j<n; j++) // negative partial sum:  - \sum_{j=i+1}^{n}{u_j \dot b_j^*}
			{
				y  -=  u[j] * U[j][i];
			}
			u[i] = floor(y + 0.5);
			// 
			if(u[i] <= y)
				d = -1; //(-1)^zi = -1
			else
				d = 1;
			// 
			long tmp; 
			if(z[i]%2 ==0)
			{
				tmp = z[i]/2;
				u[i] = u[i] + d * tmp;
			}
			else
			{
				tmp = (z[i] + 1) / 2;
				u[i] = u[i] - d * tmp;
			}
			// after determin u[i], calculating partial sum of ||v||^2
			partial += (u[i] - y) * (u[i] - y) * B1norm2_d[i];
			if(partial > target2)
			{
				findflag = false;
				break; //this z is longer than target length
			}
		}
		if(findflag)
		{
			sol_z = z;
			MUL(sol, u, B); //v = u * B
			break;
		}
	}
}

// void cellENUM_main()
// {
// 	// variants
// 	mat_ZZ B;
// 	Vec<int> pi;
// 	Vec<double> r, Bnorm2;
// 	vec_RR Bnorm2_RR;
// 	mat_RR U;
// 	long n, i, j, M;
// 	double Bound;
// 	// input
// 	ifstream finB("basis.txt");
// 	finB >> B;
// 	finB.close();
// 	n = B.NumRows();
// 	// LLL and BKZ
// 	ZZ det0;
// 	LLL(det0, B, 99, 100);
// 	BKZ_XD(B, 0.99, 20); 
// 	// GS info
// 	ComputeGS(B, U, Bnorm2_RR);
// 	Bnorm2.SetLength(n);
// 	for(i=0; i<n; i++)
// 	{
// 		Bnorm2[i] = conv<double>(Bnorm2_RR[i]);
// 	}
//  // cout<<Bnorm2_RR<<endl;
//  // cout<<Bnorm2<<endl;
// 	// get M
// 	M = 100000;
// 	// parameters
// 	double SS = SUM(Bnorm2);
// 	double minBound, maxBound;
// 	minBound = SS/4; //when tag is zero vector, E(Cell(t)) = sum{0.5^2*Bi}
// 	maxBound = 10.0 * SS; //heuristically set. I donnot know how to set the upperbound.
// 	computePrunPara(pi, r, Bnorm2, n); //r[i] = Bnorm2[pi[i]]
// 	Bound = computeBound(minBound, maxBound, pi, r, n , M);
//  // cout<<pi<<endl<<r<<endl;
// 	cout<<"Find bound R^2 = "<<Bound<<endl;
// 	if(Bound < 0)
// 	{
// 		cout<<"Upperbound is too small or there is anything wrong."<<endl;
// 		return -1;
// 	}
// 	// ENUM
// 	Vec<Vec<int>> List;
// 	List.SetLength(0);
// 	long counter;
// 	cellENUM(List, counter, pi, Bound, r, n);
// 	cout<<"Find vectors amount: "<<counter<<endl;
// 	// output
// 	ofstream fout("Tag.txt");
// 	for(i=0; i<counter; i++)
// 		fout<<List[i]<<endl;
// 	fout.close();
// }