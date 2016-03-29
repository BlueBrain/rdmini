#ifndef SAMPLE_H_
#define SAMPLE_H_

#include <cstddef>
#include <random>
#include <vector>
#include <cmath>

#include "rdmini/categorical.h"

/** Weightied random sampling algorithms
 *
 * The samplers below have an API loosely modelled on the C++ random
 * distribution concept.
 *
 * Any particular sampler S is parameterized by a type S::param_type.
 * For a param_type p, and instance s of S:
 *
 *     S::param_type s.param() const        Retrieve parameter
 *     void s.param(p)                      Set parameter
 *
 * Conventionally, if the parameter type can be constructed with
 * arguments args..., then the sampler can be constructed with
 * the same arguments.
 *
 * Samplers expose and may possibly be parameterized by two types:
 *
 *     S::size_type        Represents sample and population counts
 *     S::real_type        Represents probabilities and weight parameters
 *
 * Samples are drawn from a collection described by an iterator
 * range and written to an output iterator. Depending on the algorithm,
 * the requirements on the iterators may be more stringent, for example,
 * a reservoir algorithm might demand that the output iterator is also
 * a random access iterator.
 *
 * If [b,e) is the iterator range, o is the output iterator, and g is a
 * uniform random number generator:
 *
 *     s.sample(b,e,o,g)   Sample from [b,e) and writey to o, using g.
 *
 * The size of the sample (and the number of items thus written) may be
 * fixed or variable, depending on the sampler.
 *
 *     size_type s.min()   Minimum sample size.
 *     size_type s.max()   Maximum sample size.
 *
 * The algorithm may also require a fixed size population; if the
 * range [b,e) to the sample method is smaller than the required
 * population size, the sampler should throw a std::out_of_range
 * exception; elements in the range that are beyond than the population
 * size are ignored.
 *
 *     size_type s.size9)  Minimum population size to draw from.
 */

namespace rdmini {

/** Ordered systematic sampler
 *
 * The parameter is constructed from a sequence of inclusion probabilities
 * each in the interval [0,1].
 *
 * Each item in the range [b,e) will be drawn zero or more times, according
 * to the inclusion probability; probabilities past the end of the supplied
 * parameter are taken to be zero. The total number in the sample will be
 * equal to the sum of the weights over the range, rounded up or down.
 *
 * No particular care is taken to eliminate the effects of numerical round off
 * error, and as such even if the weights over the range sum exactly to n,
 * there is the possibility that fewer than n items will be sampled.
 */

struct ordered_systematic_sampler {
    using size_type=std::size_t;
    using real_type=double;

    struct param_type {
        param_type() {}

        template <typename Iter>
        param_type(Iter pi_begin,Iter pi_end): pi_psum(pi_begin,pi_end) {
            std::partial_sum(pi_psum.begin(),pi_psum.end(),pi_psum.begin(),
                [](real_type s,real_type x) -> real_type {
                    if (x<0.0 || x>1.0) throw std::out_of_range("invalid inclusion probability");
                    return s+x;
                });
        }

    private:
        std::vector<real_type> pi_psum;
        friend class ordered_systematic_sampler;
    };

    ordered_systematic_sampler() {}

    template <typename Iter>
    ordered_systematic_sampler(Iter pi_begin,Iter pi_end): P(pi_begin,pi_end) {}

    param_type P;

    void param(const param_type &P_) { P=P_; } 
    const param_type &param() const { return P; }

    void reset() {} // nop

    size_type min() const { return 0; }
    size_type max() const {
        if (P.pi_psum.empty()) return 0;
        else return (size_type)std::ceil(P.pi_psum.back());
    }

    size_type size() const { return 0; }

    template <typename InIter,typename OutIter,typename Rng>
    size_type sample(InIter b,InIter e,OutIter o,Rng &g) {
        static std::uniform_real_distribution<real_type> U;

        real_type u=U(g);
        size_type n=0;
        size_type n_max=max();

        for (real_type v: P.pi_psum) {
            if (b==e || n==n_max) break;

            if (u<v) {
                *o++=*b;
                u+=1;
            }
            ++b;
        }
        return n;
    }
};

/** Multinomial draw sampler
 *
 * Sample n items with replacement with unequal weights.
 *
 * The draw sampler uses a categorical distribution to draw n items
 * from the population of N. It performs well when n < N; when n is
 * significantly more than N, a sequential algorithm that draws
 * from a sequence of Binomial distributions would be preferable.
 */

struct multinomial_draw_sampler {
    using size_type=std::size_t;
    using real_type=double;

    using categorical=rdmini::categorical_distribution<size_type,real_type>;
    categorical cat;
    size_type n=0;

    multinomial_draw_sampler() {}

    template <typename Iter>
    multinomial_draw_sampler(size_type n_,Iter b,Iter e): n(n_),cat(b,e) {}

    struct param_type {
        size_type n=0;
        categorical::param_type cat_param;

        param_type() {}
        param_type(size_type n_,const categorical::param_type &p): n(n_), cat_param(p) {}

        template <typename Iter>
        param_type(size_type n_,Iter mu_begin,Iter mu_end): n(n_), cat_param(mu_begin,mu_end) {}
    };

    void param(const param_type &P) {
        n=P.n;
        cat.param(P.cat_param);
    } 

    param_type param() const { return param_type{n,cat.param()}; }

    void reset() {} // nop

    size_type min() const { return n; }
    size_type max() const { return n; }
    size_type size() const { return 1+cat.max(); }

    template <typename FwdIter,typename OutIter,typename Rng>
    size_type sample(FwdIter b,FwdIter e,OutIter o,Rng &g) {
        // draw n items from categorical distribution
        if (n==0) return 0;

        if (std::distance(b,e)<size()) throw std::out_of_range("population range too small");

        for (size_type i=0; i<n; ++i) {
            FwdIter x=std::next(b,cat(g));
            *o++=*x;
        }
        return n;
    }
};

} // namespace rdmini

#endif // ndef SAMPLE_H_
