#ifndef SSA_PP_PROCSYS_H_
#define SSA_PP_PROCSYS_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "rdmini/iterspan.h"
#include "rdmini/rdmodel.h"
#include "rdmini/tiny_multiset.h"
#include "rdmini/exceptions.h"

/** SSA process system that maintains process dependencies
 * factored through populations, and computes propensities
 * on demand from cached factors. */

template <unsigned MaxOrder=3>
struct ssa_pp_procsys {
    typedef uint32_t key_type;
    typedef double value_type;
    typedef uint32_t pop_type;
    typedef int32_t count_type;

    static constexpr size_t max_process_order=MaxOrder;
    static constexpr size_t max_population_index=std::numeric_limits<pop_index>::max()-1;
    static constexpr size_t max_count=std::numeric_limits<count_type>::max();
    static constexpr size_t max_participants=max_population_index;
    static constexpr size_t max_instances=std::numeric_limits<pop_type>::max()-1;

private:
    /** Fundamental data structures:
     *
     * pop_count:
     *     pop_count[j][p] holds population count for population index p in instance j
     *
     * pop_to_pc_tbl:
     *     pop_to_pc_tbl[p] references a collection of propensity contributions for
     *     population index p. Each contribution comprises a process index and a
     *     a slot index, which tells us how to update the propensity calculation tables.
     *     The implementation requires that a population's contributions to the same
     *     process are stored contiguously.
     *
     * rate:
     *     rate[k] is the rate constant for the process k
     *
     * propensity_tbl:
     *     propensity_tbl[j][k] is a (short) sequence of count_type values used to
     *     compute (together with rate[k]) the propensity of process k in instance j
     *
     * proc_delta_tbl:
     *     proc_delta_tbl[k] is a (short) sequence of pairs (p,d) that describe
     *     which populations p should be adjusted by a delta d when the process
     *     k is applied.
     */
    
    size_t n_pop;            // number of populations
    size_t n_proc;           // number of processes 
    size_t n_instance;       // number of instances

    std::vector<std::vector<pop_type>> pop_count;
    std::vector<value_type> rate;

    typedef std::array<count_type,max_process_order> propensity_tbl_entry;
    std::vector<std::vector<propensity_tbl_entry>> propensity_tbl;

    struct pc_entry {
        key_type k;     // process number
        unsigned index; // in range [0,MaxOrder)
    };
    std::vector<std::vector<pc_entry>> pop_to_pc_tbl;

    struct pd_entry {
        pop_type p;     // population index
        int delta;      // change to apply

        pd_entry() {}
        pd_entry(std::pair<pop_type,int> pd): p(pd.first), delta(pd.second) {}
    };
    std::vector<std::vector<pd_entry>> proc_delta_tbl;


    void initialise(size_t n) {
        n_instance=n;
        n_pop=0;
        n_proc=0;

        pop_count=std::vector<std::vector<pop_type>>(n_instance);
        propensity_tbl=std::vector<std::vector<propensity_tbl_entry>>(n_instance);

        rate.clear();
        pop_to_pc_tbl.clear();
        proc_delta_tbl.clear();
    }

public:
    explicit ssa_pp_procsys(size_t n_instance_=1) {
        initialise(n_instance_);
    }

    void reset() {
        #pragma omp parallel for
        for (size_t j=0;j<n_instance;++j) {
            for (size_t p=0;p<n_pop;++p) set_count(p,0,j);
        }
    }

    template <typename ProcDesc>
    void add(const ProcDesc &q) {
        if (n_proc>=std::numeric_limits<key_type>::max())
            throw rdmini::invalid_value("process index out of bounds");

        key_type key=n_proc++;
        // check for unsigned overflow ...
        if (n_proc<key) throw std::overflow_error("number of processes overflow");

        // in minimizing assumptions about q.left() and q.right()
        // datastructures, aim for one-pass through the supplied info.

        small_map<pop_type,int> proc_delta_entry;
        std::array<count_type,max_process_order> left_sorted;
        size_t nleft=0;

        pop_type max_pop=n_pop;
        for (auto p: q.left()) {
            if (nleft>max_process_order)
                throw rdmini_invalid_value("too many reactants");

            --proc_delta_entry[p];
            if (proc_delta_entry.size()>max_participants)
                throw rdmini_invalid_value("too many participants");

            left_sorted[nleft++]=p;
            if (p>max_pop) max_pop=p;
        }
        std::sort(left_sorted.data(),left_sorted.data()+nleft);

        for (auto p: q.right()) {
            ++proc_delta_entry[p];
            if (proc_delta_entry.size()>max_participants)
                throw rdmini_invalid_value("too many participants");

            if (p>max_pop) max_pop=p;
        }

        // extend population-indexed data structures if required
        if (max_pop>n_pop) {
            #pragma omp parallel for
            for (size_t j=0;j<n_instance;++j) {
                pop_count[j].resize(max_pop);
            }

            pop_to_pc_tbl.resize(max_pop);
            n_pop=max_pop;
        }

        // update proc_delta_tbl:
        proc_delta_tbl.emplace_back(proc_delta_entry.begin(),proc_delta_entry.end());

        // update pop_to_pc_tbl and propensity_tbl
        propensity_tbl_entry prop_entry;
        std::fill(prop_entry.begin(),prop_entry.end(),1);

        count_type c=0; // population contribution to propensity
        pop_type p_prev=0;
        for (size_t i=0;i<nleft;++i) {
            pop_type p=left_sorted[i];
            pop_to_pc_tbl[p].emplace_back(key,i);

            if (i==0 || p!=p_prev) c=pop_count[p];
            else --c;

            prop_entry[i]=c;
        }

        #pragma omp parallel for
        for (size_t j=0;j<n_instance;++j) 
            propensity_tbl[j].push_back(prop_entry);

        rate.push_back(q.rate());

    }

    template <typename In>
    void add(In b,In e) {
        while (b!=e) add(*b++);
    }

    /** Remove all processes, population counts */
    void clear() {
        proc_propensity_tbl.clear();
        for (auto &entry: pop_contribs_tbl) entry.clear();

        pop_count.assign(n_pop,0);
        rate.clear();
    }

    size_t size() const { return n_proc; }
    
    count_type count(size_t p) const { return pop_count[p]; }

    template <typename F>
    void set_count(size_t p,count_type c,F update_notify,size_t j=0) {
        for (auto kci: pop_to_pc_tbl[p])
            apply_contrib_update(kci,c-pop_count[p],update_notify,j);
        pop_count[p]=c;
    }

    void set_count(size_t p,count_type c,size_t j=0) { set_count(p,c,[](key_type) {},j); }

    template <typename F>
    void apply(key_type k,F update_notify,size_t j=0) {
        for (auto pd: proc_delta_tbl[k]) {
            for (auto kci: pop_to_pc_tbl[pd.p])
                apply_contrib_update(kci,pd.delta,update_notify,j);
            pop_count[pd.p]+=pd.delta;
        }
    }

    void apply(key_type k,size_t j=0) { apply(k,[](key_type) {},j); }

    value_type propensity(key_type k,size_t j=0) {
        const proc_propensity_entry &kp=proc_propensity_tbl[j][k];
        value_type r=rate[k];
        for (auto c: kp.counts) r*=c;
        return r;
    }

    friend std::ostream& operator<<(std::ostream &O,const ssa_pp_procsys &sys) {
        // primarily for debugging
        O << "ssa_pp_procsys: n_pop=" << sys.n_pop << ", n_proc=" << sys.n_proc << "\n";
        O << "pop_to_pc_tbl:\n";
        size_t idx=0;
        for (const auto &e: sys.pop_to_pc_tbl) {
            O << "    " << std::setw(6) << std::right << idx++ << ":";
            for (const auto &kci: e) 
                O << ' ' << kci.k  << ':'
                  << std::showpos << kci.i << std::noshowpos;
            O << "\n";
        }
        O << "proc_delta_tbl:\n";
        idx=0;
        for (const auto &e: sys.proc_delta_tbl) {
            O << "    " << std::setw(6) << std::right << idx++ << ":";
            for (const auto &pd: e) 
                O << ' ' << pd.p  << ':'
                  << std::showpos << pd.delta << std::noshowpos;
            O << "\n";
        }
        O << "rate:\n";
        idx=0;
        for (const auto &r: sys.rate) {
            O << "    " << std::setw(6) << std::right << idx++ << ":"
              << ' ' << r << '\n';
        }
        return O;
    }
};


#endif // ndef  SSA_PP_PROCSYS_H_
