#include "psucc.h"

NTL_CLIENT
using namespace std;


void makeListOdd(Vec<Vec<int>> & List, const long & n)
{
    long j, k;
    long K = List.length();
    for(k=0; k<K; k++)
    {
        // List[k]
        j=n-1;
        while((j) && (!List[k][j])) // j>0 and List[k][j]==0
        {
            j--;
        }
        if(List[k][j]%2 == 0)
            List[k][j] -= 1;
    }
}

void makeListEven(Vec<Vec<int>> & List, const long & n)
{
    long j, k;
    long K = List.length();
    for(k=0; k<K; k++)
    {
        // List[k]
        j=n-1;
        while((j) && (!List[k][j])) // j>0 and List[k][j]==0
        {
            j--;
        }
        if(List[k][j]%2 == 1)
            List[k][j] += 1;
    }
}

double P_succ_wrap(Vec<Vec<int>> & List, const long & n, const Vec<double> & B1norm2_d, \
	const double & target_d, const long & M, const int sound)
{
	double a = max(50.0, 30.0 + 3*sqrt((double)n));
	long K = 40, J = 30;
	Vec<int> cell_max;
	Vec<Vec<Vec<complex<double>> > > ERFtable;
	double psucc;
	Vec<double> B1norm_d;
	B1norm_d.SetLength(n);
	for(long i=0; i<n; i++)
		B1norm_d[i]= sqrt(B1norm2_d[i]);
	// ===========================
	if(List.length()==0)
	{
		cerr<<"ERROR in p_succ: No candidate vector!"<<endl;
		return 0;
	}
	if(sound)
	{
		psucc = estimatePsuccEnsemble_sound(List, B1norm_d, B1norm2_d, M, target_d, a, n, K, J);
	}
	else // use method in [AN17] but there's some flaw in it, so we use harmonic mean value to heuristically repair
	{
		double psucc1, psucc2;
		makeListOdd(List, n); 
		findMaxCell(cell_max, List, n);
		ERF_buildTable(ERFtable, B1norm_d, target_d, cell_max, a, K+J+1, n);
		psucc1 = estimatePsuccEnsemble(ERFtable, List, B1norm_d, B1norm2_d, M, target_d, a, n, K, J);
		makeListEven(List, n);
		findMaxCell(cell_max, List, n);
		ERF_buildTable(ERFtable, B1norm_d, target_d, cell_max, a, K+J+1, n);
		psucc2 = estimatePsuccEnsemble(ERFtable, List, B1norm_d, B1norm2_d, M, target_d, a, n, K, J);
		psucc = 2/(1/psucc1+1/psucc2);
		makeListOdd(List, n);
	}
	return psucc;
}

//  ======== quick estimation =======
void findMaxCell(Vec<int>& cell_max, const Vec<Vec<int>>& List, const long & n)
{
	long N = List.length();
	cell_max = List[0];
	for(long i = 1; i<N; i++)
	{
		for(long j =0; j<n; j++)
		{
			if(cell_max[j]<List[i][j])
				cell_max[j] = List[i][j];
		}
	}
}

void ERF_buildTable(Vec<Vec<Vec<complex<double>> > >& ERFtable, const Vec<double>& GSnorm, const double& R, const Vec<int>& max_t, \
	const double& a, const long & max_m, const long& n)
// erf(alpha_i * sqrtz(a+(m-0.5)PI*I))
{
	long m, i, ti;
	complex<double> root_s, z;
	ERFtable.SetLength(max_m+1);
	for(m = 1; m <= max_m; m++)
	{
		ERFtable[m].SetLength(n);
		for(i = 0; i < n; i++)
		{
			ERFtable[m][i].SetLength(max_t[i]+3); //0,1,...,max_t[i]+1
		}
	}

	for (long m=1; m <= max_m; m++)
	{
		root_s = sqrtz_sm(a, m);
		for(i=0; i<n; i++)
		{
			for(ti=0; ti<=(max_t[i]+2); ti++)
			{
				if(ti==0)
				{
					ERFtable[m][i][ti] = complex<double>(0, 0);
				}
				else
				{
					ERFtable[m][i][ti] = erfz(mul_db_cplx((double)ti*GSnorm[i]/(2*R), root_s));
				}
			}
		}
	}
}

double DPE_ImF_sm(const Vec<Vec<Vec<complex<double>> > >& ERFtable, const Vec<double>& GSnorm, const double& R, const Vec<int>& t, 
	const double& a, const long& m, const long& n)
{
	complex<double> F(1, 0);
	complex<double> Fi(1, 1);
	complex<double> coeff;
	double alpha_i, beta_i;
	for(long i=0; i<n; i++)
	{
		// 存储alpha_i, beta_i, 还可以优化
		// alpha_i = (double)t[i]*GSnorm[i]/(2*R);
		// beta_i = (double)(t[i]+1)*GSnorm[i]/(2*R); 
		// 2 (beta[i] - alpha[i]) = GSnorm[i]/(R)
		Fi = div_cplx_db(ERFtable[m][i][(t[i]+1)] - ERFtable[m][i][(t[i])], GSnorm[i]/R);
 // printf("m=%ld, i=%ld, erf(beta_i√s)=%.8g+%.8gI, erf(alpha_i√s)=%.8g+%.8gI, prod increment = %.8g+%.8gI\n",\
 	m, i, ERFtable[m][i][(t[i]+1)].real(), ERFtable[m][i][(t[i]+1)].imag(), ERFtable[m][i][(t[i])].real(),\
 	ERFtable[m][i][(t[i])].imag(), Fi.real(), Fi.imag());
		F = F*Fi;
 // printf("partial prod = %.8g+%.8gI\n", F.real(), F.imag());
	}
	coeff = EQ1_cplx_Fm(a, m, n);
	F = F * coeff;
 // printf("coeff = %.8g+%.8gI, F=%.8g+%.8gI\n", coeff.real(), coeff.imag(), F.real(), F.imag());
	return F.imag();
}

double DPE_psucc(Vec<Vec<Vec<complex<double>> > >& ERFtable, const Vec<double>& GSnorm, const double& R, const Vec<int>& t, const double& a, const long& n, 
	const long& K, const long& J)
{
	double Dmnt = 0; //dominant term
	double Tail = 0;
	long m, j, l;
	Vec<double> Fm;
	Fm.SetLength(K + J + 1 + 1); // ignore Fm[0], contains (K+J+1) terms in total
	for(m=1; m<=K+J+1; m++)
	{
		Fm[m] = DPE_ImF_sm(ERFtable, GSnorm, R, t, a, m, n);
	}

	//The first K terms
	for(m = 1; m<=K; m++)
	{
		if(m%2 == 1)
			Dmnt -= Fm[m];
		else
			Dmnt += Fm[m];
 // printf("partial sum =%.16g, increment = %.16g\n", Dmnt, Fm[m]);
	}

	// van-Wijingaarden transformation using Fm[K+1+0] to Fm[K+1+J} terms
	Mat<double> S;
	S.SetDims(J+1, J+1);
	for(j = 0; j<=J; j++)
	{
		S[0][j] = 0; 
		for(m = 0; m<= j; m++)
		{
			if(m%2==1)
				S[0][j] -= Fm[m+K+1];
			else
				S[0][j] += Fm[m+K+1];
		}
	}
	for(l=1; l<=J; l++)
	{
		for(j=0; j<=J-l; j++)
		{
			S[l][j] = (S[l-1][j] + S[l-1][j+1])/2.0;
		}
	}
	if((K+1)%2==1)
		Tail = -S[J][0];
	else
		Tail = S[J][0];
 // printf("Tail sum =%.16g\n", Tail);
	return exp(a) * (Dmnt + Tail);

}

double estimatePsuccEnsemble(Vec<Vec<Vec<complex<double>> > >& ERFtable, const Vec<Vec<int>>& List, const Vec<double>& GSnorm, const Vec<double>& GSnorm2_d,\
	const long& M, const double& R, const double& a, const long& n, const long& K, const long& J)
{
	long i, j;
	double psucc_cell, psucc;
	long N = List.length();
	// quick sort List according to Expectation
	Vec<Vec<int>> List_sort;
	Vec<double> Ev;
	List_sort = List;
	Ev.SetLength(List.length());
	for(i=0; i<N; i++)
	{
		Ev[i] = 0;
		for(j=0; j<n; j++)
		{
			Ev[i] += ((double)List[i][j]+0.5)*((double)List[i][j]+0.5)*GSnorm2_d[j];
		}
	}
	quickSort(List_sort, Ev, 0, N-1);
	
	//stratified sampling
	N = M < N? M : N;
	long setNum, sampleNum;
	long start, end, split;
	if(N/500 < 2000)
		sampleNum = 2000;
	else if((N/500) > 8000)
		sampleNum = 8000;
	else
		sampleNum = (N/500);
	//mid(10000, M/500, 2000)
	split = N / sampleNum;
	start = 0; 
	end = split;
	sampleNum = 0;
	psucc=0;
	while(start<N)
	{
		// sample one tag in List[start, end)
		sampleNum += 1;
		i = start + rand()%split -1;
		if(i<0) i=0;
		if(i>=N) i=N-1;
 // double psucct;
 // psucct = DPE_psucc(ERFtable, GSnorm, R, List[i],  a, n, K, J);
 // psucc += (double)(end-start) * psucct;	
		psucc += (double)(end-start) * DPE_psucc(ERFtable, GSnorm, R, List_sort[i],  a, n, K, J);
 // printf("i=%ld, psucc=%.8g\t",i, psucct);
		start += split;
		end += split;
		end = end > N ? (N-1): end;
	}
 // double psucct;
 // psucct = 0;
 // for(i=0; i<N; i++)
 // 	psucct+=DPE_psucc(ERFtable, GSnorm, R, List[i],  a, n, K, J);
 // printf("the real psucc=%.8g\n", psucct);
	// return (psucc)/sampleNum;
	return psucc;
}


// ==== concrete psucc calculation ======
double DPE_ImF_sm_sound(const Vec<double>& GSnorm, const double& R, const Vec<int>& t, 
	const double& a, const long& m, const long& n)
{
	complex<double> F(1, 0);
	complex<double> Fi(1, 1);
	complex<double> coeff, root_s, A, B;
	double alpha_i, beta_i;
	root_s = sqrtz_sm(a, m);
	for(long i=0; i<n; i++)
	{
		// --- 存储alpha_i, beta_i, 还可以优化 ---
		alpha_i = (double)t[i]*GSnorm[i]/(2*R);
		beta_i = (double)(t[i]+1)*GSnorm[i]/(2*R); 
		// --- 2 (beta[i] - alpha[i]) = GSnorm[i]/(R) ---
 // printf("i= %ld partial prod = %.8g+%.8gI\n", i, F.real(), F.imag());
		if(t[i]==0)
			A = complex<double>(0, 0);
		else
			A = erfz(mul_db_cplx((double)t[i]*GSnorm[i]/(2*R), root_s));
		B = erfz(mul_db_cplx((double)(t[i]+1)*GSnorm[i]/(2*R), root_s));
		Fi = div_cplx_db(B-A, GSnorm[i]/R);
 // printf("m=%ld, i=%ld, erf(beta_i√s)=%.8g+%.8gI, erf(alpha_i√s)=%.8g+%.8gI, prod increment = %.8g+%.8gI\n",\
 	m, i, ERFtable[m][i][(t[i]+1)].real(), ERFtable[m][i][(t[i]+1)].imag(), ERFtable[m][i][(t[i])].real(),\
 	ERFtable[m][i][(t[i])].imag(), Fi.real(), Fi.imag());
		F = F*Fi;
 // printf("i= %d partial prod = %.8g+%.8gI\n", i, F.real(), F.imag());
	}
	coeff = EQ1_cplx_Fm(a, m, n);
	F = F * coeff;
 // printf("coeff = %.8g+%.8gI, F=%.8g+%.8gI\n", coeff.real(), coeff.imag(), F.real(), F.imag());
	return F.imag();
}


double DPE_psucc_sound(const Vec<double>& GSnorm, const double& R, const Vec<int>& t, const double& a, const long& n, 
	const long& K, const long& J)
{
	double Dmnt = 0; //dominant term
	double Tail = 0;
	long m, j, l;
	Vec<double> Fm;
	Fm.SetLength(K + J + 1 + 1); // ignore Fm[0], contains (K+J+1) terms in total

	// === modify R and t ===
	double Rtrunc, utrunc, atrunc;
	Vec<int> ttrunc;
	long trunc=n-1;
	while(t[trunc]==0 && trunc>0)
		trunc --;
	// --- now t[trunc] is the last non-zero component of t ---
	if(t[trunc]%2==0)
		utrunc= -t[trunc]/2;
	else
		utrunc = (t[trunc]+1)/2;
	Rtrunc = R*R - utrunc*GSnorm[trunc]*utrunc*GSnorm[trunc];
	Rtrunc = sqrt(Rtrunc);
	atrunc = max(50.0, 30.0 + 3*sqrt((double)trunc));

	for(m=1; m<=K+J+1; m++)
	{
		Fm[m] = DPE_ImF_sm_sound(GSnorm, Rtrunc, t, atrunc, m, trunc);
	}

	//The first K terms
	for(m = 1; m<=K; m++)
	{
		if(m%2 == 1)
			Dmnt -= Fm[m];
		else
			Dmnt += Fm[m];
 // printf("partial sum =%.16g, increment = %.16g\n", Dmnt, Fm[m]);
	}

	// van-Wijingaarden transformation using Fm[K+1+0] to Fm[K+1+J} terms
	Mat<double> S;
	S.SetDims(J+1, J+1);
	for(j = 0; j<=J; j++)
	{
		S[0][j] = 0; 
		for(m = 0; m<= j; m++)
		{
			if(m%2==1)
				S[0][j] -= Fm[m+K+1];
			else
				S[0][j] += Fm[m+K+1];
		}
	}
	for(l=1; l<=J; l++)
	{
		for(j=0; j<=J-l; j++)
		{
			S[l][j] = (S[l-1][j] + S[l-1][j+1])/2.0;
		}
	}
	if((K+1)%2==1)
		Tail = -S[J][0];
	else
		Tail = S[J][0];
 // printf("Tail sum =%.16g\n", Tail);
	return exp(atrunc) * (Dmnt + Tail);
}


double estimatePsuccEnsemble_sound(const Vec<Vec<int>>& List, const Vec<double>& GSnorm, const Vec<double>& GSnorm2_d,\
	const long& M, const double& R, const double& a, const long& n, const long& K, const long& J)
{
	long i, j;
	double psucc_cell, psucc;
	long N = List.length();
	// quick sort List according to Expectation
	Vec<Vec<int>> List_sort;
	Vec<double> Ev;
	List_sort = List;
	Ev.SetLength(List.length());
	for(i=0; i<N; i++)
	{
		Ev[i] = 0;
		for(j=0; j<n; j++)
		{
			Ev[i] += ((double)List[i][j]+0.5)*((double)List[i][j]+0.5)*GSnorm2_d[j];
		}
	}
	quickSort(List_sort, Ev, 0, N-1);
	
	// // --- stratefied sampling ------
	// N = M < N? M : N;
	// long setNum, sampleNum;
	// long start, end, split;
	// if(N/500 < 2000)
	// 	sampleNum = 2000;
	// else if((N/500) > 8000)
	// 	sampleNum = 8000;
	// else
	// 	sampleNum = (N/500);
	// split = N / sampleNum;
	// start = 0; 
	// end = split;
	// sampleNum = 0;
	// psucc=0;
	// while(start<N)
	// {
	// 	// sample one tag in List[start, end)
	// 	sampleNum += 1;
	// 	i = start + rand()%split -1;
	// 	if(i<0) i=0;
	// 	if(i>=N) i=N-1;
	// 	psucc += (double)(end-start) * DPE_psucc_sound(GSnorm, R, List_sort[i],  a, n, K, J);
	// 	start += split;
	// 	end += split;
	// 	end = end > N ? (N-1): end;
	// }
	// // ------ counting one by one ------
	psucc = 0;
	for(i=0; i<N; i++)
		psucc += DPE_psucc_sound(GSnorm, R, List[i],  a, n, K, J);
	return psucc;
}





//  =========== functional tools ================
complex<double> EQ1_cplx_Fm(const double& a, const long& m, const long& n) 
// compute PI^{n/2}/s^{n/2+1} = sqrt(PI/s)^n / s
//in the DPE FILT problem, s = (a , (m-0.5)PI), sqrt(s) should be in [0, PI/4]
{
	double nd = (double)n;
	double x, y, rs, rt, rz, thetas, thetat, thetaz;
	complex<double> s(a, ((double)m - 0.5)*PI);
	rs = sqrt(s.real()*s.real() + s.imag()*s.imag());
	thetas = arg(s);

	complex<double> root_s; //sqrtz(s)
	root_s = sqrtz_sm(a, m);
	
	x = root_s.real();
	y = root_s.imag();
	complex<double> t( PIroot * x / (x*x + y*y), -PIroot * y / (x*x + y*y)); //t = sqrt(pi)/sqrt(s)
	rt = sqrt(t.real()*t.real() + t.imag()*t.imag());
	thetat = atan2(t.imag(), t.real()); //arg(t)
	rz = pow(rt, nd) / rs;
	thetaz = thetat * nd - thetas; //z = t^n / s
// printf("pi^{n/2}/s^{n/2-1} = %.16g + %.16g*I, arg=%.10g\n", rz*cos(thetaz), rz*sin(thetaz), atan2(rz*sin(thetaz),rz*cos(thetaz)));
	return complex<double>(rz*cos(thetaz), rz*sin(thetaz));
}

complex<double> sqrtz_sm(double a, long m)
//|arg(z)| < pi/4
// z = (a, (m-0.5)Pi)
{
	// complex<double> z(a, (m - 0.5)*PI); 
	//arg(z) is in first quadrant
	double b, x, y, r;
	b = ((double)m - 0.5)*PI;
	r = sqrt(a*a+b*b);
	x = sqrt((r+a)/2);
	y = sqrt((r-a)/2); //sgn(b) = '+'
	return complex<double>(x, y);
}

void quickSort(Vec<Vec<int>>& List, Vec<double>& Ev, long low, long high)
{
	if(low<high)
	{
		// pivot
		long pivot, start, end, mid;
		double tmpEv;
		Vec<int> tmpVec;
		start = low;
		end = high;
		// median-of-three
		mid = low + (high - low)/2;
		if(Ev[mid] > Ev[high]){
			SWAP(Ev[mid], Ev[high]);
			SWAP(List[mid], List[high]);
		}
		if(Ev[low] > Ev[high]){
			SWAP(Ev[low], Ev[high]);
			SWAP(List[low], List[high]);
		}
		if(Ev[mid] > Ev[low]){
			SWAP(Ev[mid], Ev[low]);
			SWAP(List[mid], List[low]);
		}
		tmpEv = Ev[low];
		tmpVec = List[low];
		// 
		while(start<end)
		{
			while(start < end && Ev[end] >= tmpEv){
				end--;
			}
			if(start < end){
				Ev[start] = Ev[end];
				List[start] = List[end];
			}
			while(start < end && Ev[start] <= tmpEv){
				start++;
			}
			if(start < end){
				Ev[end] = Ev[start];
				List[end] = List[start];
			}
		}
		Ev[start] = tmpEv;
		List[start] = tmpVec;
		pivot = start;
		// iteration
		quickSort(List, Ev, low, pivot-1);
		quickSort(List, Ev, pivot+1, high);
	}
}
