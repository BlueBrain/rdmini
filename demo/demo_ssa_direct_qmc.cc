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


//-----------------------------------------------------------------------------
int main(int argc, char **argv) {
   
    // Random generator and distributions
    std::minstd_rand R;
    rdmini::Linear_RNG Rlin;
    std::uniform_real_distribution<double> U(0.5,1);
    rdmini::vc_uniform_real_distribution<double> U_vdc(0.,1.);

    // Generate random propensities of the form a_j = u_j * 2^-j
    constexpr size_t n_proc=10;
    std::vector<double> prop(n_proc);
    for (size_t i=0; i<n_proc; ++i) {
        prop[i]=i+1;//std::ldexp(U(R),-i);
    }
    std::shuffle(prop.begin(),prop.end(),R);

    // Sum of all propensities, a_0 = sum_j a_j
    double total=std::accumulate(prop.begin(),prop.end(),0.0);

    // Instantiate a ssa_direct object and add all propensities
    typedef rdmini::ssa_direct<size_t,double> S;
    S ssa(prop.size());
    for (size_t i=0; i<prop.size(); ++i)
        ssa.update(i,prop[i]);
    
    // Compute exact mean of the event indexes
    double exact_mean = 0.0;
    for (size_t j=0; j<prop.size(); ++j)
        exact_mean += double(j) * double(prop[j]) / total ;

    // Generate several events e_1,...,e_N and compute partial sums 
    constexpr size_t n_events = 1e6;
    std::vector<size_t> partial_sums(n_events);
    partial_sums[0] = ssa.inverse_cdf(U_vdc(Rlin));
    for (size_t l=1; l<n_events; ++l) 
        partial_sums[l] = partial_sums[l-1] + ssa.inverse_cdf(U_vdc(Rlin));
    
    // Compute approx means
    std::vector<double> appr_mean(n_events);
    for (size_t l=0; l<n_events; ++l) {
        appr_mean[l] = double(partial_sums[l])/(l+1) ;
    }

    // Compute theoretical bound on the discrepancy of the VdC sequence:
    // see P.Kritzer, "A new upper bound on the star discrepancy of 
    // (0,1)-sequences", Electronic journal of combinatorial number theory 5(3),
    // 2005.
    // Dstar(x1,...,xN) <= f_b * logN/N + c_b * 1/N
    constexpr size_t base = 10;
    constexpr double a_b  = (base%2)?((base-1.)/4.):((base*base)/(4.*(base+1.)));
    constexpr double f_b  = a_b/std::log(base);
    constexpr double c_b  = (2.0>1.0+1.0/base+a_b)?(2.0):(1.0+1.0/base+a_b);
    
    // Compute total variation of the integrand F^{-1}(u), to obtain an error
    // estimate using Koksma-Hlawka inequality
    constexpr double V_f = n_proc; // total variation

    // Check if the theoretical error estimate hold in this case
    // err_N <= V_f * Dstar
    bool flag_success = true;
    for (size_t N=0; N<n_events; ++N) {
            if (std::abs(appr_mean[N]-exact_mean)>V_f*(f_b*std::log(N+1)/(N+1)+c_b/(N+1))) {
            std::cout << "After " << N << " iterations" 
                << " expected max error " << f_b*std::log(N+1)/(N+1)+c_b/(N+1)
                << " but observed error " << std::abs(appr_mean[N]-exact_mean) << std::endl;
            break;
            flag_success = false; 
        }
    }

    std::cout << "Test " << (flag_success?"passed!":"failed!") << std::endl;

    return 0;
}

