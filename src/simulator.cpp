#include "gsa.h"
#include "simulator.h"
#include "tools.h"
#include "cellenum.h"
#include "psucc.h"
#include "reduction.h"
#include <iostream>
#include <cmath>


NTL_CLIENT
using namespace std;

double qmean[35] = {0.9683072236392071, 0.9693395227829056, 0.9692174139255462, 0.9701276187796496, 0.969944832086452,\
 0.9707132398274835, 0.9706377331212223, 0.9713153410601959, 0.9712211497785079, 0.9718226874663856, 0.9717787383713438, \
 0.9723172742374745, 0.9722572424929301, 0.9727658561266128, 0.9727289528585941, 0.9731841268230101, 0.9732056023699784, \
 0.9736233234456441, 0.973622007955314, 0.9740366998244139, 0.9739786951748196, 0.9743989519412324, 0.9744600769743728, \
 0.9748890613469775, 0.9749352666407893, 0.9753210213436875, 0.9753748136338195, 0.9756242402865503, 0.97556320325273, \
 0.9758327541208127, 0.975863034206352, 0.9760323481264208, 0.9760542998392626, 0.9762285674169696, 0.9763178533202796};

Vec<mat_ZZ> BKZList_Basis;
Vec<double> BKZList_Time;
bool NM_SAVEBKZ; //should be set by DP_NelderMead()

// ============= Main simulator ==================
costAll DP_NelderMead(const mat_ZZ& B0, const long & n, const ZZ& det, const double &target_d, long &M_op,  int& bkz_op, \
    int DP_SIM_FLAG, int k, bool NM_VERBOSE)
{
    //return time
    long i, j;
    int bkz[3]; // blocksize 
    long M[3];
    int bkz_mid, bkz_rfl, bkz_e, bkz_c;
    long M_mid, M_rfl, M_e, M_c;
    double Time_rfl, Time_e, Time_c;
    double Time[3];
    costAll Cost[3];
    costAll Cost_temp, Cost_rfl, Cost_e, Cost_c;
    mat_ZZ B;
    if(DP_SIM_FLAG==1 || DP_SIM_FLAG ==2)
    {
        if(B0.NumRows()!= n)
        {
            cerr << "DP_NelderMead: Basis dimension dismatch!" <<endl;
            Cost_temp.rounds=0;
            Cost_temp.oneround=0;
            return Cost_temp;
        }
        B = B0;
    }
    // for(i=0; i<3; i++)
    // {
    //     bkz[i] = rand()%60+20; // beta \in [20, 60] at current version for n<200
    //     if(bkz[i]%2==0)
    //         bkz[i] -= 1;       // beta is forced to be odd
    //     M[i] = rand()%1e7 + 1e4; // M \in [10000, 10000000]
    // }
    // --- set initial points far away -------
    BKZList_Time.SetLength(n);
    BKZList_Basis.SetLength(n);
    NM_SAVEBKZ = true;
    clear(BKZList_Time);
    bkz[0] = 25; 
    bkz[1] = 31;
    bkz[2] = 40;
    M[0]   = 10000;
    M[1]   = 1e5;
    M[2]   = 5e5; // the largest $M$ should not be too large in case of the memory explosion
    while(1)
    {
        for(i=0; i<3; i++)
        {
            Cost[i] = DP_simulator_warp(B, n, det, target_d, M[i], bkz[i], DP_SIM_FLAG, k);
            Time[i] = Cost[i].rounds * Cost[i].oneround + Cost[i].preproc;
        }
        // sorting
        for(i=0; i<2; i++)
        {
            for(j=0; j<2-i; j++)
            {
                if(Time[j]>Time[j+1])
                {
                    swap(bkz[j], bkz[j+1]);
                    swap(M[j],M[j+1]);
                    swap(Time[j], Time[j+1]);
                }
            }
        }
        if(NM_VERBOSE)
        {
            cout << "+++++ New iteration : initial points +++++" <<endl;
            cout <<"bkz:  "<< bkz[0]<<'\t'<<bkz[1]<<'\t'<<bkz[2] <<endl;
            cout <<"M:    "<< M[0]<<'\t'<<M[1]<<'\t'<<M[2] <<endl;
            cout <<"Cost: "<< Time[0]<<'\t'<<Time[1]<<'\t'<<Time[2] <<endl;
        }
        //
        if( (abs(Time[0]-Time[2])<0.01*Time[0]) || (abs(bkz[2]-bkz[0])<=3 && abs(M[2]-M[0])<0.05*M[0])\
            || (abs(bkz[2]-bkz[0])<=2 && abs(M[2]-M[0])<10000) )
        // return the optimal values
        {
            bkz_op = bkz[0];
            M_op = M[0];
            // Cost_c = DP_sim_pure(det, target_d, M_op, n, k, bkz_op);
            return Cost[0];
        }
        bkz_mid = round((bkz[0] + bkz[1])/2);
        bkz_rfl = 2*bkz_mid - bkz[2];
        M_mid   = round((M[0] + M[1])/2);
        M_rfl   = 2*M_mid - M[2];
        if(NM_VERBOSE)
        {
            cout << "----- middle results -----" <<endl;
            cout <<"middle point:     "<< bkz_mid<<'\t'<<M_mid<<endl;
            cout <<"reflection point: "<< bkz_rfl<<'\t'<<M_rfl<<'\t'<<endl;
        }
        Cost_rfl = DP_simulator_warp(B, n, det, target_d, M_rfl, bkz_rfl, DP_SIM_FLAG, k);
        Time_rfl = Cost_rfl.oneround * Cost_rfl.rounds + Cost_rfl.preproc;
        if(NM_VERBOSE)
        {
            cout <<"reflection point Cost: "<< Time_rfl<<endl;
        }
        // 
        if(Time[0] <= Time_rfl && Time_rfl < Time[1]) //Reflection
        {
            if(NM_VERBOSE) {cout<<"Go to reflection."<<endl;}
            bkz[2] = bkz_rfl;
            M[2]   = M_rfl;
            continue;
        }
        else if (Time_rfl < Time[0]) // Expasion
        {
            if(NM_VERBOSE) {cout<<"Go to expasion."<<endl;}
            bkz_e  = bkz_mid + 2*(bkz_rfl - bkz_mid);
            M_e    = M_mid + 2*(M_rfl - M_mid);
            Cost_e = DP_simulator_warp(B, n, det, target_d, M_e, bkz_e, DP_SIM_FLAG, k);
            Time_e = Cost_e.oneround * Cost_e.rounds + Cost_e.preproc;
            if(Time_e < Time_rfl)
            {
                bkz[2] = bkz_e;
                M[2]   = M_e;
            }
            else
            {
                bkz[2] = bkz_rfl;
                M[2]   = M_rfl;
            }
            continue;
        }
        else if (Time[1] <= Time_rfl && Time_rfl<= Time[2]) //Contraction
        {
            bkz_c = round(bkz_mid + (bkz_rfl - bkz_mid)/2);
            M_c   = round(M_mid + (M_rfl - M_mid)/2);
            Cost_c = DP_simulator_warp(B, n, det, target_d, M_c, bkz_c, DP_SIM_FLAG, k);
            Time_c = Cost_c.oneround * Cost_c.rounds + Cost_c.preproc;
            if(Time_c < Time[2])
            {
                if(NM_VERBOSE) {cout<<"Go to contraction."<<endl;}
                bkz[2] = bkz_c;
                M[2]   = M_c;
                continue;
            }
            // else , goto shrink
        }
        else if (Time_rfl > Time[2]) // Contraction
        {
            bkz_c = round(bkz_mid + (bkz[2] - bkz_mid)/2);
            M_c   = round(M_mid + (M[2] - M_mid)/2);
            Cost_c = DP_simulator_warp(B, n, det, target_d, M_c, bkz_c, DP_SIM_FLAG, k);
            Time_c = Cost_c.oneround * Cost_c.rounds + Cost_c.preproc;
            if(Time_c < Time[2])
            {
                if(NM_VERBOSE) {cout<<"Go to contraction."<<endl;}
                bkz[2] = bkz_c;
                M[2]   = M_c;
                continue;
            }
        }
        if(NM_VERBOSE) {cout<<"Go to shrink."<<endl;}
        for(i=1; i<3; i++) //Shrink
        {
            bkz[i] = round(bkz[0] + (bkz[i]-bkz[0])/2);
            M[i] = round(M[0] + (M[i]-M[0])/2);
        }
        if(NM_VERBOSE)
        {
            cout << "----- iterated points -----" <<endl;
            cout <<"bkz:  "<< bkz[0]<<'\t'<<bkz[1]<<'\t'<<bkz[2] <<endl;
            cout <<"M:    "<< M[0]<<'\t'<<M[1]<<'\t'<<M[2] <<endl;
        }
    } //END WHILE
}

costAll DP_simulator_warp(const mat_ZZ& B_instance, const long& n, const ZZ& det, const double& target_d, \
    long M, int bkz,
    int DP_SIM_FLAG, int TOURS)
{
    /*
    There are three methods provided: 
    DP_SIM_FLAG = 0 -> only use simulated GS sequnce to estimate the running time;
    DP_SIM_FLAG = 1 -> use the BKZ_beta-reduced(2.0 ver) input basis `B_test` to estimate p_succ;
    DP_SIM_FLAG = 2 -> do one round of DP enum to record the time, and use the reduced `B_test` to estimate p_succ;
    */
    costAll re;
    if(DP_SIM_FLAG==0)
    {
        re = DP_sim_pure(det, target_d, M, n, TOURS, bkz);
    }
    else if(DP_SIM_FLAG ==1)
    {
        re = DP_sim_hybrid(B_instance, target_d, M, n, TOURS, bkz, FT_DEFAULT);
    }
    else if(DP_SIM_FLAG == 2)
    {
        re = DP_sim_basis(B_instance, target_d, M, n, TOURS, bkz, FT_DEFAULT);
    }
    return re;
}

costAll DP_sim_pure(const ZZ& det, const double& target_d, const long& M, const long& n, int k, int bkz)
// only use fitting equations
{
    if( (M<5000) || (bkz<0) || (bkz>(n/2)) ) 
    //we restrict beta < n^{1-1/4} for now. if BKZ_simulator is available in the later version, we will change it to $beta<n$
    {
        costAll penalty;
        // penalty for NM method
        penalty.rounds   = 1;
        penalty.oneround = NM_MAX;
        return penalty;
    }
    costAll Time;
    long i, j;
    //  === simulate GS ===
    Vec<double> GS2_fitted;
    double q, logvol, logb0;
    GS2_fitted.SetLength(n);
    if(bkz>=10 && bkz <45)
        q = qmean[bkz-10];
    else
        q = 1-exp(-0.008668967071869243*bkz-3.405597365733596);
    logvol = conv<double>(log(det));
    logb0 = (logvol - (double)n*(n-1)/2*log(q)) / (double)n;

    // === padding GS2 ======
    GS2_fitted[0] = exp(2*logb0);
    for(i=1; i<n; i++)
    {
        GS2_fitted[i] = exp(2*(log(q)*i + logb0));
    }
    //  === calculate T_reduce ===
    RR r_block, block_RR;
    double T_reduce = 0;
    RR nodes;
    clear(nodes);
    for(i=0; i<n-1; i++)
    {
        long block = bkz<(n-i)?bkz:(n-i);
        vec_RR GS2_local;
        GS2_local.SetLength(block);
        for(long t=0; t<block; t++)
            GS2_local[t] = conv<RR>(GS2_fitted[i+t]);
        block_RR = conv<RR>(block);
        r_block = exp(0.354605*block_RR*log(block_RR) -1.533128*block_RR  + 4.898155*log(block_RR) - 2.908438); 
            //if you use full-enum, r_block=1
        nodes += fullEnumCost(block, SqrRoot(GS2_local[0]), GS2_local)/r_block;
           // + 0.000904381*pow(block,3)*n*n + 28752188;
    }

    T_reduce = conv<double>(nodes)*k*205.45;
    //  === calculate T_bin ===
    double T_bin;
    T_bin = 0.11341*M*n*n + 13.155*M*n*log(n) + 265.65*M - 84679*n + 15455380;
    //  === calculate T_cell ===
    double T_cell;
    T_cell = 2.4339*n*M + 108.74*M-17455*n + 1334139;
    //  === calculate T_recover ===
    double T_recover;
    T_recover = 0.39045*n*n + 167.06*n - 4350.4;
    T_recover = T_recover * M;
    // === calculate psucc ===
    double p_succ;
    p_succ = 1/SimulateCost_q(det, bkz, n, target_d, M); 
    if(p_succ<1) {p_succ=1;}
    Time.preproc = 0;
    Time.rounds = p_succ;
    Time.oneround = (T_reduce + T_bin + T_cell + T_recover);
    return Time;
}

costAll DP_sim_hybrid(const mat_ZZ& B0, const double& target_d, const long& M, const long& n, int k, int bkz, FloatType float_type)
// given a target (and reduced) lattice only for calculating p_succ
{
    if( (M<5000) || (bkz<0) || (bkz>(n/2)) )
    {
        costAll penalty;
        // penalty for NM method
        penalty.rounds   = 1;
        penalty.oneround = NM_MAX;
        return penalty;
    }
    if(B0.NumRows()!= n)
    {
        cerr << "DP_sim_hybrid: Basis dimension dismatch!" <<endl;
        costAll Time;
        Time.rounds = 0;
        Time.oneround = 0;
        return Time;
    }
    long i, j;
    int PREC = RR::precision();
    struct timespec time1, time2;
    // === processing lattice basis ===
    mat_ZZ B;
    B = B0;
    mat_RR U;
    vec_RR GS2;
    Vec<double> GS2_d;
    GS2.SetLength(n);
    GS2_d.SetLength(n);
    if(BKZList_Time[bkz]>0 && NM_SAVEBKZ)
    {
        B = BKZList_Basis[bkz];
    }
    else
    {
        clock_gettime(CLOCK_MONOTONIC, &time1);
        if(bkz>n/2)
        {
            LLL_default(B, 0.99, FT_DEFAULT, PREC);
            BKZ_quick(B, bkz, FT_DEFAULT);
        }
        else
        {
            LLL_default(B, 0.99, FT_DEFAULT, PREC);
            BKZ_default(B, bkz, FT_DEFAULT);
        }
        clock_gettime(CLOCK_MONOTONIC, &time2);
        if(NM_SAVEBKZ)
        {
            BKZList_Time[bkz] = timeDiff(time1, time2);
            BKZList_Basis[bkz] = B;
        }
    }
    ComputeGS(B, U, GS2);
    conv(GS2_d, GS2);

    //  === calculate T_reduce ===
    RR r_block, block_RR;
    double T_reduce = 0;
    RR nodes;
    clear(nodes);
    for(i=0; i<n-1; i++)
    {
        long block = bkz<(n-i)?bkz:(n-i);
        vec_RR GS2_local;
        GS2_local.SetLength(block);
        for(long t=0; t<block; t++)
            GS2_local[t] = GS2[i+t];
        block_RR = conv<RR>(block);
        r_block = exp(0.354605*block_RR*log(block_RR) -1.533128*block_RR  + 4.898155*log(block_RR) - 2.908438); 
            //if you use full-enum, r_block=1
        nodes += fullEnumCost(block, SqrRoot(GS2_local[0]), GS2_local)/r_block;
    }
    T_reduce = conv<double>(nodes)*k*205.45;
    //  === calculate T_bin ===
    double T_bin;
    T_bin = 0.11341*M*n*n + 13.155*M*n*log(n) + 265.65*M - 84679*n + 15455380;
    //  === calculate T_cell ===
    double T_cell;
    T_cell = 2.4339*n*M + 108.74*M-17455*n + 1334139;
    //  === calculate T_recover ===
    double T_recover;
    T_recover = 0.39045*n*n + 167.06*n - 4350.4;
    T_recover = T_recover * M;
    // === calculate psucc ===
    double p_succ = 0;
    long counter = 0;
    double NNRBound; // NNR expected value's upperbound
    double minBound, maxBound;
    Vec<int> pi; //coordinate permutation
    Vec<double> r;

    Vec<Vec<int>> List;
    minBound = 0;
    maxBound = 0;
    // --- compute maxbound ----
    for(i=0; i<n; i++)
    {
        double zi = pow(M, 1/(double)n);
        zi = ceil(zi);
        maxBound += (zi+1) * zi * GS2_d[i];
    }
    // ------
    computePrunPara(pi, r, GS2_d, n);
    NNRBound = computeBound(minBound, maxBound, pi, r, n, M);
    List.SetLength(0);
    cellENUM(List, counter, pi, NNRBound, r, n);
    p_succ += 1/P_succ_wrap(List, n, GS2_d, target_d, counter);
    //
    if(p_succ<1) {p_succ=1;}
    costAll Time;
    if(NM_SAVEBKZ)
        Time.preproc = BKZList_Time[bkz];
    else
        Time.preproc = timeDiff(time1, time2);
    Time.rounds = p_succ;
    Time.oneround = (T_reduce + T_bin + T_cell + T_recover);
    return Time;
}

costAll DP_sim_basis(const mat_ZZ& B0, const double& target_d, const long& M, const long& n, int k, int bkz, FloatType float_type)
// given a target (and reduced) lattice and conduct a round of real calculation
{
    if( (M<5000) || (bkz<0) || (bkz>(n/2)) )
    {
        costAll penalty;
        // penalty for NM method
        penalty.rounds   = 1;
        penalty.oneround = NM_MAX;
        return penalty;
    }
    if(B0.NumRows()!= n)
    {
        cerr << "DP_sim_basis: Basis dimension dismatch!" <<endl;
        costAll Time;
        Time.rounds = 0;
        Time.oneround = 0;
        return Time;
    }
    long i, j;
    int PREC = RR::precision();
    struct timespec timestart, timeend;
    struct timespec time1, time2;
    // === processing lattice basis ===
    mat_ZZ B;
    B = B0;
    mat_RR U;
    Mat<double> U_d;
    Mat<long> B_l; 
    vec_RR GS2;
    Vec<double> GS2_d;
    GS2.SetLength(n);
    GS2_d.SetLength(n);
    if(BKZList_Time[bkz]>0 && NM_SAVEBKZ)
    {
        B = BKZList_Basis[bkz];
    }
    else
    {
        clock_gettime(CLOCK_MONOTONIC, &time1);
        if(bkz>n/2)
        {
            LLL_default(B, 0.99, FT_DEFAULT, PREC);
            BKZ_quick(B, bkz, FT_DEFAULT);
        }
        else
        {
            LLL_default(B, 0.99, FT_DEFAULT, PREC);
            BKZ_default(B, bkz, FT_DEFAULT);
        }
        clock_gettime(CLOCK_MONOTONIC, &time2);
        if(NM_SAVEBKZ)
        {
            BKZList_Time[bkz] = timeDiff(time1, time2);
            BKZList_Basis[bkz] = B;
        }
    }
    ComputeGS(B, U, GS2);
    conv(GS2_d, GS2);
    conv(U_d, U);
    conv(B_l, B);
    double target2_d = target_d*target_d;

    //  === calculate one round of DP-ENUM ===
    clock_gettime(CLOCK_MONOTONIC, &timestart);
    long counter = 0;
    double NNRBound; // NNR expected value's upperbound
    double minBound, maxBound;
    Vec<int> pi; //coordinate permutation
    Vec<double> r;
    Vec<Vec<int>> List;
    Vec<long> sol;
    Vec<int> sol_z;
    // --- re-preprocessing ---
    rePermutation_keephead(B, 0);
    BKZ_k_tours(B, bkz, k, float_type);
    ComputeGS(B, U, GS2);
    conv(GS2_d, GS2);
    conv(U_d, U);
    conv(B_l, B);
    // --- compute maxbound ----
    minBound = 0;
    maxBound = 0;
    for(i=0; i<n; i++)
    {
        double zi = pow(M, 1/(double)n);
        zi = ceil(zi);
        maxBound += (zi+1) * zi * GS2_d[i];
    }
    computePrunPara(pi, r, GS2_d, n);   //r[i] = Bnorm2[pi[i]]
    NNRBound = computeBound(minBound, maxBound, pi, r, n , M); //binary search try how many times
    DP_ENUM_Wrap_sim(counter, pi, NNRBound, r, n, sol, sol_z, B_l, GS2_d, U_d, target2_d);

    clock_gettime(CLOCK_MONOTONIC, &timeend);

    // === calculate psucc ===
    double p_succ = 0;
    // ------
    // computePrunPara(pi, r, GS2_d, n);
    // NNRBound = computeBound(minBound, maxBound, pi, r, n, M);
    List.SetLength(0);
    cellENUM(List, counter, pi, NNRBound, r, n);
    // harmonic mean value
    p_succ += 1/P_succ_wrap(List, n, GS2_d, target_d, counter);
    //
    if(p_succ<1) {p_succ=1;}
    costAll Time;
    if(NM_SAVEBKZ)
        Time.preproc = BKZList_Time[bkz];
    else
        Time.preproc = timeDiff(time1, time2);
    Time.rounds = p_succ;
    Time.oneround = timeDiff(timestart, timeend);
    return Time;
}

//  === modified version of DP_ENUM_Wrap as a subalgorithm of DP_sim_basis() =====
// === cell ENUM + decoding with out insertion ===
void DP_ENUM_Wrap_sim(long & counter, const Vec<int>& pi, double Bound, const Vec<double>& r, long n, \
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
                        // break;
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


// costAll DP_sim_pure_fullenum(const ZZ& det, const double& target_d, const long& M, const long& n, int k, int bkz)
// // only use fitting equations
// {
//     costAll Time;
//     double cycles = 0;
//     long i, j;
//     //  === simulate GS ===
//     Vec<double> GS2_fitted;
//     double q, logvol, logb0;
//     GS2_fitted.SetLength(n);
//     if(bkz>=10 && bkz <45)
//         q = qmean[bkz-10];
//     else
//         q = 1-exp(-0.0046334765*bkz-3.49696026590);
//     logvol = conv<double>(log(det));
//     logb0 = (logvol - (double)n*(n-1)/2*log(q)) / (double)n;
//     // === padding GS2 ======
//     GS2_fitted[0] = exp(2*logb0);
//     for(i=1; i<n; i++)
//     {
//         GS2_fitted[i] = exp(2*(log(q)*i + logb0));
//     }
//     //  === calculate T_reduce ===
//     double r_block;
//     double T_reduce = 0;
//     RR nodes;
//     clear(nodes);
//     for(i=0; i<n-1; i++)
//     {
//         long block = bkz<(n-i)?bkz:(n-i);
//         vec_RR GS2_local;
//         GS2_local.SetLength(block);
//         for(long t=0; t<block; t++)
//             GS2_local[t] = conv<RR>(GS2_fitted[i+t]);
//         //if you use full-enum, r_block=1
//         nodes += fullEnumCost(block, SqrRoot(GS2_local[0]), GS2_local);
//             //+ 0.000904381*pow(block,3)*n*n + 28752188; THERE ARE some problems in this part
//     }
//     T_reduce = conv<double>(nodes)*k*205.45;
//     //  === calculate T_bin ===
//     double T_bin;
//     T_bin = 0.11341*M*n*n + 13.155*M*n*log(n) + 265.65*M - 84679*n + 15455380;
//     //  === calculate T_cell ===
//     double T_cell;
//     T_cell = 2.4339*n*M + 108.74*M-17455*n + 1334139;
//     //  === calculate T_recover ===
//     double T_recover;
//     T_recover = 0.39045*n*n + 167.06*n - 4350.4;
//     T_recover = T_recover * M;
//     // === calculate psucc ===
//     double p_succ;
//     p_succ = 1/SimulateCost_q(det, bkz, n, target_d, M) + 1/SimulateCost_q_even(det, bkz, n, target_d, M);
//     p_succ = p_succ/2;
//     if(p_succ<1) {p_succ=1;}
//     Time.rounds = p_succ;
//     Time.oneround = T_reduce + T_bin + T_cell + T_recover;
//     return Time;
// }

// costAll findMinCost_pure(const ZZ& det, const double& target_d, const long& n, \
//     long& M, int &bkz, int& k)
// {
//     costAll minCycles;
//     costAll tempCycles;
//     long Mt;
//     int bkzt;

//     minCycles = DP_sim_pure(det, target_d, 10000, n, k, 10);
//     for(bkzt=10; bkzt<n/2; bkzt+=3)
//     {
//         for(Mt=10000; Mt<200000; Mt+=5000)
//         {
//             tempCycles = DP_sim_pure(det, target_d, Mt, n, k, bkzt);
//         // cout<<Mt<<'\t'<<bkzt<<'\t'<<tempCycles<<endl;
//             if(tempCycles.oneround*tempCycles.rounds < minCycles.oneround*minCycles.rounds)
//             {
//                 M = Mt;
//                 bkz= bkzt;
//                 minCycles = tempCycles;
//             }
//         }
//     }
//     return minCycles;
// }