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

    // Generate random propensities of the form a_j = u_j * 2^-j
    size_t n_proc=10;
    std::vector<double> prop(n_proc);
    for (size_t i=0; i<n_proc; ++i) {
        prop[i]=i+1;//std::ldexp(U(R),-i);
    }
    std::shuffle(prop.begin(),prop.end(),R);

    // Sum of all propensities, a_0 = sum_j a_j
    double total=std::accumulate(prop.begin(),prop.end(),0.0);

    // Compute cdf of indexes
    std::vector<double> cdf(n_proc);
    cdf[0] = prop[0]/total;
    for (size_t i=1; i<prop.size(); ++i)
        cdf[i] = cdf[i-1] + prop[i]/total;

    // Instantiate a ssa_direct object and add all propensities
    typedef rdmini::ssa_direct<size_t,double, rdmini::vc_uniform_real_distribution<double> > S;
    S ssa(prop.size());
    for (size_t i=0; i<prop.size(); ++i)
        ssa.update(i,prop[i]);
    
    // Compute exact mean of the event indexes
    double exact_mean = 0.0;
    for (size_t j=0; j<prop.size(); ++j)
        exact_mean += double(j) * double(prop[j]) / total ;

    // Generate several events e_1,...,e_N and compute partial sums 
    size_t n_events = 1e8;
    std::vector<size_t> partial_sums(n_events);
    partial_sums[0] = (ssa.next(Rlin)).key();
    for (size_t l=1; l<n_events; ++l)
        partial_sums[l] = partial_sums[l-1] + (ssa.next(Rlin)).key();
    
    // Compute approx means
    std::vector<double> appr_mean(n_events);
    for (size_t l=0; l<n_events; ++l) {
        appr_mean[l] = double(partial_sums[l])/(l+1) ;
    }

    for (size_t N=n_events-100; N<n_events; ++N)
        std::cout << "(Measured error) / (logN/N) = " 
            << std::abs(appr_mean[N]-exact_mean) / ( std::log(N+1) / (N+1)) 
            << std::endl;



    
    return 0;
}

