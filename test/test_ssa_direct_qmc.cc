#include <string>
#include <cstring>
#include <cstddef>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <random>

#include "rdmini/vandercorput.h"
#include "rdmini/ssa_direct.h"

#include "gtest/gtest.h"

//class ssa_distribution : public ::testing:Test {
//protected: 
//    virtual void SetUp () {


namespace KHineq {
    // Compute theoretical bound on the discrepancy of the VdC sequence:
    // see P.Kritzer, "A new upper bound on the star discrepancy of 
    // (0,1)-sequences", Electronic journal of combinatorial number theory 2005:
    //      Dstar(x1,...,xN) <= f_b * logN/N + c_b * 1/N
    // This value is then used in Koksma-Hlawka inequality.
    constexpr size_t base = 10;
    constexpr double a_b  = (base%2)?((base-1.)/4.):((base*base)/(4.*(base+1.)));
    constexpr double f_b  = a_b/std::log(base);
    constexpr double c_b  = (2.0>1.0+1.0/base+a_b)?(2.0):(1.0+1.0/base+a_b);
}


TEST(SsaDistribution, MomentTest) {
    
    // Random generator and distributions
    std::minstd_rand R;
    rdmini::Linear_RNG Rlin;
    std::uniform_real_distribution<double> U(0.5,1);
    rdmini::vdc_uniform_real_distribution<double> U_vdc(0.,1.);

    // Populate vector of propensities and compute the sum
    constexpr size_t n_proc = 20;
    std::vector<double> prop(n_proc);
    for (size_t i=0; i<n_proc; ++i) 
        prop[i] = i+1;//std::ldexp(U(R),-i);
    std::shuffle(prop.begin(),prop.end(),R);
    double total = std::accumulate(prop.begin(),prop.end(),0.0);

    // Compute exact moments of the event indexes distribution
    double exact_mu1 = 0.;
    double exact_mu2 = 0.;
    for (size_t j=0; j<n_proc; ++j) {
        exact_mu1 += double(j) * prop[j] / total ;
        exact_mu2 += double(j*j) * prop[j] / total ;
    }

    // Instantiate a ssa_direct object and add all propensities
    rdmini::ssa_direct<size_t,double> ssa(prop.size());
    for (size_t i=0; i<prop.size(); ++i) ssa.update(i,prop[i]);
    
    // Generate several events and compute approximations of mu1 and mu2. Then,
    // compute the error with respect to the exact value and check whether this
    // error satisfies the theoretical error estimate by Koksma-Hlawka
    //      err_N <= V_fun * Dstar(x1,...,xN)
    // where V_fun is the total variation of the integrand.
    constexpr double V_f_mu1 = n_proc;          // total variation for mu1
    constexpr double V_f_mu2 = n_proc*n_proc;   // total variation for mu2
    constexpr size_t n_events = 1000*1000;
    std::vector<double> appr_mu1(n_events);
    size_t idx_ssa = 0.;
    double approx_mu1 = 0.;
    double approx_mu2 = 0.;
    idx_ssa = ssa.inverse_cdf(U_vdc(Rlin));
    approx_mu1 += idx_ssa;
    approx_mu2 += idx_ssa*idx_ssa;
    bool KH_mu1 = true, KH_mu2 = true;
    std::ostringstream os1, os2;
    for (size_t N=1; N<n_events && (KH_mu1 || KH_mu2); ++N) {
        idx_ssa = ssa.inverse_cdf(U_vdc(Rlin));
        if (KH_mu1) {
            approx_mu1 = (approx_mu1*N + idx_ssa)/(N+1); 
            if (std::abs(approx_mu1-exact_mu1)>V_f_mu1*(KHineq::f_b*std::log(N+1)/(N+1)+KHineq::c_b/(N+1))) {
                os1 << "After " << N << "iterations:"
                    << " expected max error " << KHineq::f_b*std::log(N+1)/(N+1)+KHineq::c_b/(N+1)
                    << " but observed error " << std::abs(approx_mu1-exact_mu1) << std::endl;
                KH_mu1 = false;
            }
        }
        if (KH_mu2) {
            approx_mu2 = (approx_mu2*N + idx_ssa)/(N+1); 
            if (std::abs(approx_mu2-exact_mu1)>V_f_mu2*(KHineq::f_b*std::log(N+1)/(N+1)+KHineq::c_b/(N+1))) {
                os2 << "After " << N << "iterations:"
                    << " expected max error " << KHineq::f_b*std::log(N+1)/(N+1)+KHineq::c_b/(N+1)
                    << " but observed error " << std::abs(approx_mu2-exact_mu2) << std::endl;
                KH_mu2 = false;
            }
        }
    }
            
    EXPECT_TRUE(KH_mu1) << os1.str();
    EXPECT_TRUE(KH_mu2) << os2.str();
}

