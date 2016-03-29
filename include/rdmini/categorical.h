#ifndef CATEGORICAL_H_
#define CATEGORICAL_H_

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cstddef>
#include <sstream>

#include "rdmini/util/ios_util.h"

/** Categorical distribution implemented with alias method.
 *
 * ref: Vose (1991), A linear algorithm for generating random numbers
 *      with a given distribution. IEEE Transactions on Software Engineering 17(9),
 *      972–975. doi:10.1109/32.92917
 */

namespace rdmini {

template <typename IntType,typename RealType=double>
struct categorical_distribution {
    using value_type=IntType;
    using real_type=RealType;

    struct param_type {
        param_type() {}

        template <typename Iter>
        param_type(Iter b,Iter e): tbl(b,e) {
            value_type n=static_cast<value_type>(tbl.size());

            // normalise q values
            real_type sum=0;
            for (value_type i=0; i<n; ++i) sum+=q(i);

            real_type scale=real_type(n)/sum;
            for (value_type i=0; i<n; ++i) q(i)*=scale;

            // build alias table ...
            value_type i_small=0,i_big=0;
            while (i_small<n && !(q(i_small)<=1)) ++i_small;
            while (i_big<n && q(i_big)<=1) ++i_big;

            value_type i_small_top=i_small;
            
            while (i_small<n && i_big<n) {
                alias(i_small)=i_big;
                q(i_big)=(q(i_big)+q(i_small))-1;

                // advance i_small from i_small_top unless we've just
                // made a new small behind i_small_top:

                bool new_small=(q(i_big)<=1);
                if (new_small && i_big<i_small_top) i_small=i_big;
                else {
                    i_small=i_small_top+1;
                    while (i_small<n && !(q(i_small)<=1)) ++i_small;
                    i_small_top=i_small;
                }

                // advance i_big if we've made a new small
                if (new_small) {
                    while (i_big<n && q(i_big)<=1) ++i_big;
                }
            }

            // anything left over should be given probability 1

            if (i_small<n) { 
                q(i_small)=1;
                i_small=i_small_top+1;
                while (i_small<n) q(i_small++)=1;
            }

            while (i_big<n) q(i_big++)=1;
        }

        bool operator==(const param_type &P) const { return tbl==P.tbl; }
        bool operator!=(const param_type &P) const { return tbl!=P.tbl; }

        friend std::ostream &operator<<(std::ostream &O,const param_type &p) {
            scoped_ios_format save(O);
            
            O << std::setprecision(std::numeric_limits<real_type>::max_digits10);

            O << p.size();
            for (auto t: p.tbl) { O << ' ' << t.first << ' ' << t.second; }
            return O;
        }

        friend std::istream &operator>>(std::istream &I,param_type &p) {
            scoped_ios_format save(I);
            value_type n=0;

            I >> n;
            if (!I) return I;
            if (n<0) { I.setstate(std::ios::failbit); return I; }

            std::vector<table_entry> tbl(n);
            for (value_type i=0; i<n; ++i) {
                table_entry t;
                //I >> std::skipws >> t.first >> std::skipws >> t.second;
                I >> t.first >> t.second;
                if (!I) return I;
                if (t.first<0 || t.first>1 || t.second<0 || t.second>=n) {
                    I.setstate(std::ios::failbit);
                    return I;
                }
                tbl[i]=t;
            }
            p.tbl=tbl;

            return I;
        }

        const real_type &q(value_type i) const { return tbl[i].first; }
        const value_type &alias(value_type i) const { return tbl[i].second; }
        value_type size() const { return static_cast<value_type>(tbl.size()); }

    private:
        // tbl[i].first = probability bin i should give i instead of alias
        // tbl[i].second = alias for bin i
        struct table_entry: std::pair<real_type,value_type> {
            table_entry(real_type q_=0,value_type a_=0):
                std::pair<real_type,value_type>(q_,a_) {}
        };
        std::vector<table_entry> tbl;

        real_type &q(value_type i) { return tbl[i].first; }
        value_type &alias(value_type i) { return tbl[i].second; }
    };

    const param_type &param() const { return P; }
    void param(const param_type &P_) { P=P_; }

    value_type min() const { return 0; }
    value_type max() const { return P.size()-1; }

    void reset() {}

    bool operator==(const categorical_distribution &them) const { return P==them.P; }
    bool operator!=(const categorical_distribution &them) const { return !(P==them.P); }

    categorical_distribution() {}
    explicit categorical_distribution(const param_type &P_): P(P_) {}

    template <typename Iter>
    categorical_distribution(Iter b,Iter e): P(b,e) {}

    template <typename Rng>
    value_type operator()(Rng &g,const param_type &p) {
        value_type n(p.size());
        if (n==0) return 0;

        std::uniform_real_distribution<real_type> draw(0,n);
        real_type d=draw(g);

        value_type bin=value_type(d); // integer part
        real_type u=d-bin;

        return u<p.q(bin)?bin:p.alias(bin);
    }

    template <typename Rng>
    value_type operator()(Rng &g) { return (*this)(g,P); }

    friend std::ostream &operator<<(std::ostream &O,const categorical_distribution &C) {
        return O << C.P;
    }

    friend std::istream &operator>>(std::istream &I,categorical_distribution &C) {
        return I >> C.P;
    }

private:
    param_type P;
};

} // namespace rdmini

#endif // ndef CATEGORICAL_H_
