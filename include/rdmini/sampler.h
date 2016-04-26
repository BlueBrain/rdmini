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
 *     size_type s.size()  Minimum population size to draw from.
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

namespace impl {
    /** Generic order reservoir sampling implementation
     *
     * \param n    Reservoir size (unsigned integral type)
     * \param b    Beginning of population range (input iterator)
     * \param e    End of population range (input iterator)
     * \param o    Reservoir (random access output iterator)
     * \param f    Order-generating functor
     * \return     Number or items stored in reservoir
     */

    template <typename size_type,typename InIter,typename OutIter,typename F>
    size_type order_reservoir_sample(size_type n,InIter b,InIter e,OutIter o,F f) {
        using order_type=typename std::result_of<F()>::type;

        if (n==0) return 0;

        using key=std::pair<order_type,size_type>;

        std::vector<key> heap;
        heap.reserve(n);

        size_type i=0;
        for (; i<n && b!=e; ++i) {
            heap.emplace_back(f(),i);
            o[i]=*b++;
        }
        
        if (i<n) return i;

        // heapify
        std::make_heap(heap.begin(),heap.end());

        // keep least P.n elements...
        for (; b!=e; ++b) {
            order_type q=f();
            if (q<heap.front().first) {
                key k{q,heap.front().second};
                std::pop_heap(heap.begin(),heap.end());
                heap.back()=k;
                o[k.second]=*b;
                std::push_heap(heap.begin(),heap.end());
            }
        }

        return n;
    }
}

/** Adjusted Pareto sampler
 *
 * Sample n items without replacement using an order method
 * with ranking variables Q_i = U_i/(1-U_i)·(1-p_i)/p_i·a_i
 * where the p_i are the desired inclusion probabilities and
 * a_i is an adjustment factor:
 *     a_i = exp(p_i(1-p_i)(p_i-1/2)/d^2)
 *     d = sum p_i(1-p_i)
 *
 * The actual inclusion probabilities of this sampling method
 * are only approximated by p_i, but approach asymptotically
 * as d approaches infinity.
 *
 * The supplied inclusion probabilities should sum to n.
 */

struct adjusted_pareto_sampler {
    using size_type=std::size_t;
    using real_type=double;

    adjusted_pareto_sampler() {}

    template <typename Iter>
    adjusted_pareto_sampler(size_type n_,Iter b,Iter e): P(n_,b,e) {}

    struct param_type {
        size_type n=0;
        std::vector<real_type> qcoef;

        param_type() {}

        template <typename Iter>
        param_type(size_type n_,Iter pi_begin,Iter pi_end): n(n_), qcoef(pi_begin,pi_end) {
            real_type d=0;
            for (real_type q: qcoef) d+=q*(1-q);

            real_type ood2=1/(d*d);
            for (real_type &q: qcoef) {
                double loga=q*(1-q)*(q-real_type(0.5))*ood2;
                double a=1+loga+0.5*loga*loga; // approximates exp(loga) as loga is small.
                q=(1-q)/q*a;
            }
        }
    } P;

    template <typename Rng>
    struct next_order {
        std::uniform_real_distribution<real_type> U;
        const std::vector<real_type> &q;
        size_type i;
        Rng &g;
        
        explicit next_order(const param_type &P,Rng &g_): q(P.qcoef),i(0),g(g_) {}

        real_type operator()() {
            if (i>=q.size()) return std::numeric_limits<real_type>::max();
            real_type r=q[i++];
            real_type u=U(g);
            return u*r/(1-u);
        }
    };
 
    void param(const param_type &P_) { P=P_; }
    param_type param() const { return P; }

    void reset() {} // nop

    size_type min() const { return P.n; }
    size_type max() const { return P.n; }
    size_type size() const { return P.n; }

    template <typename InIter,typename OutIter,typename Rng>
    size_type sample(InIter b,InIter e,OutIter o,Rng &g) {
        next_order<Rng> f(P,g);
        return impl::order_reservoir_sample(P.n,b,e,o,f);
    }
};

/* Efraimidis and Spirakis sampler
 *
 * Reservoir implementation of Efraimidis and Spirakis (2006)
 * method for weighted sampling without replacement
 * (doi: 10.1016/j.ipl.2005.11.003).
 *
 * Parameters are not the inclusion probabilities, but the proportional
 * weights to draw an item each round. For weights not too distant
 * from n/N (n being sample size, N the population size), this approximates
 * the inclusion probabilities.
 */

struct efraimidis_spirakis_sampler {
    using size_type=std::size_t;
    using real_type=double;

    efraimidis_spirakis_sampler() {}

    template <typename Iter>
    efraimidis_spirakis_sampler(size_type n_,Iter b,Iter e): P(n_,b,e) {}

    struct param_type {
        size_type n=0;
        std::vector<real_type> oolambda;

        param_type() {}

        template <typename Iter>
        param_type(size_type n_,Iter pi_begin,Iter pi_end): n(n_), oolambda(pi_begin,pi_end) {
            for (auto &l: oolambda) l=1/l;
        }
    } P;

    template <typename Rng>
    struct next_order {
        std::exponential_distribution<real_type> E;
        const std::vector<real_type> &q;
        size_type i;
        Rng &g;
        
        explicit next_order(const param_type &P,Rng &g_): q(P.oolambda),i(0),g(g_) {}

        real_type operator()() {
            if (i>=q.size()) return std::numeric_limits<real_type>::max();
            return E(g)*q[i++];
        }
    };
 
    void param(const param_type &P_) { P=P_; }
    param_type param() const { return P; }

    void reset() {} // nop

    size_type min() const { return P.n; }
    size_type max() const { return P.n; }
    size_type size() const { return P.n; }

    template <typename InIter,typename OutIter,typename Rng>
    size_type sample(InIter b,InIter e,OutIter o,Rng &g) {
        next_order<Rng> f(P,g);
        return impl::order_reservoir_sample(P.n,b,e,o,f);
    }
};


} // namespace rdmini

#endif // ndef SAMPLE_H_
