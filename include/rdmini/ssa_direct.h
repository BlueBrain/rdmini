#ifndef SSA_DIRECT_H_
#define SSA_DIRECT_H_

#include <random>

#include "rdmini/exceptions.h"

/** Implementation of 'direct' SSA method. */

namespace rdmini {

// KeyType              must be unsigned integral
// ValueType            must be floating point
// RealDistribution     must be a real distribution over an interval [a,b]

template <typename KeyType,
         typename ValueType> 
struct ssa_direct {
    typedef KeyType key_type;
    typedef ValueType value_type;
    typedef value_type& reference;

private:
    size_t n_key;
    std::uniform_real_distribution<value_type> U;
    std::exponential_distribution<value_type> E;
    std::vector<value_type> propensities;
    value_type total;

public:
    // Nested class Event_type: an Event_type is a pair (idx, dt)
    class event_type: std::pair<key_type,value_type> {
        typedef std::pair<key_type,value_type> pair_type;
    public:
        using typename pair_type::first_type;
        using typename pair_type::second_type;

        event_type() =default;
        event_type(key_type k_,value_type dt_): pair_type(k_,dt_) {}

        key_type key() const { return this->first; }
        value_type dt() const { return this->second; }
    };

    // Constructor    
    explicit ssa_direct(size_t n_key_=0) : U(0.,1.), E(1.) { reset(n_key_); }

    // Getter for number of keys
    size_t size() const { return n_key; }

    // Computes inverse CDF
    key_type inverse_cdf(value_type u) const {
        value_type x = u*total;    
        key_type i=0;
        for (i=0; i<n_key; ++i) {
            x-=propensities[i];
            if (x<0) break;
        }
        if (i>=n_key) throw rdmini::ssa_error("fell off propensity ladder (rounding?)");
        return i;
    }

    // Computes next event: which one (idx) and when (dt) 
    template <typename R>
    event_type next(R &g) {
        return event_type{inverse_cdf(U(g)), E(g)/total};
    }
        

    // Setter for number of keys
    void reset(size_t n_key_) {
        n_key=n_key_;
        propensities.assign(n_key,0.0);
        total=0.0;
    }

    // Setter for propensity with index k 
    void update(key_type k,value_type r) {
        value_type &p=propensities[k];
        total+=r-p;
        p=r;
    }

    //Getter for propensity with index k
    value_type propensity(key_type k) const { return propensities[k]; }

    // Getter for total of propensities 
    value_type total_propensity() const { return total; };

};

} // namespace rdmini

#endif // ndef SSA_DIRECT_H_
