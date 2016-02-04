#ifndef SSA_PP_PROCSYS_H_
#define SSA_PP_PROCSYS_H_

#include <cstdint>
#include <cmath>
#include <vector>
#include <array>
#include <map>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <iostream>
#include <iomanip>

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

    void extend_n_pop(n) {
        // increase table sizes to accomodate a (larger) number of populations n.
        if (n<=n_pop) return;

        #pragma omp parallel for
        for (size_t j=0;j<n_instance;++j) {
            pop_count[j].resize(n);
        }

        pop_to_pc_tbl.resize(n);
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
        const auto &left=q.left();
        const auto &right=q.right();

        if (left.size()>max_process_order)
            throw rdmini_invalid_value("too many reactants");

        if (left.size()+right.size()>max_participants)
            throw rdmini_invalid_value("too many participants");

        if (n_proc>=std::numeric_limits<key_type>::max())
            throw rdmini::invalid_value("process index out of bounds");

        key_type key=n_proc++;
        // check for unsigned overflow ...
        if (n_proc<key) throw std::overflow_error("number of processes overflow");

        tiny_multiset<pop_type,max_process_order> lpop;
        for (auto p: left) lpop.insert(p);

        std::array<count_type,max_process_order> propensity_template;
        propensity_template.fill(1);


        // extend pop_to_pc_tbl, pop_count to account for new populations
        ...
        /*
        pop_type pop_max=std::max_element(left.begin(),left.end());
        pop_max=std::max(std::max_element(right.begin(),rig
        if (lpop_max>=n_pop) extend_n_pop(lpop_max+1);
        */

        


        proc_propensity_tbl.resize(n_proc);
        proc_propensity_tbl[key].rate=proc.rate();

        std::vector<pop_index> left_sorted;
        std::map<pop_index,count_type> delta_map;
        for (auto p: proc.left()) {
            if (p>=n_pop) throw rdmini::invalid_value("population index out of bounds");
            --delta_map[p];
            left_sorted.push_back(p);
        }
        for (auto p: proc.right()) {
            if (p>=n_pop) throw rdmini::invalid_value("population index out of bounds");
            ++delta_map[p];
        }

        proc_delta_tbl.resize(n_proc);
        auto &pd_entry=proc_delta_tbl[key];
        for (auto pd: delta_map)
           if (pd.second) pd_entry.push_back({pd.first,pd.second});

        std::sort(left_sorted.begin(),left_sorted.end());
        uint32_t index=0;
        for (auto p: left_sorted) {
            pop_contribs_tbl[p].emplace_back(key,index++);
        }
        
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
    }

    size_t size() const { return n_proc; }
    
    /** Zero all population counts, reset propensity contributions */
    void zero_populations() {
        // init contribs and pop counts
        pop_index n_pop=pop_count.size();
        for (pop_index p=0; p<n_pop; ++p) {
            pop_count[p]=0;
            auto &pc_entry=pop_contribs_tbl[p];
            if (pc_entry.empty()) continue;
 
            // note that indices in pc are grouped by kproc id.
            auto pc_i=pc_entry.begin();
            auto pc_end=pc_entry.end();

            count_type count=0;
            key_type k=pc_i->k;
            proc_propensity_tbl[k].counts[pc_i->i]=count;
            
            while (++pc_i!=pc_end) {
                if (pc_i->k==k) --count;
                else count=0;

                k=pc_i->k;
                proc_propensity_tbl[k].counts[pc_i->i]=count;
            }
        }
    }

    count_type count(size_t p) const { return pop_count[p]; }

    template <typename F>
    void set_count(size_t p,count_type c,F update_notify) {
        for (auto kci: pop_contribs_tbl[p])
            apply_contrib_update(kci,c-pop_count[p],update_notify);
        pop_count[p]=c;
    }

    void set_count(size_t p,count_type c) { set_count(p,c,[](key_type) {}); }

    template <typename F>
    void apply(key_type k,F update_notify) {
        for (auto pd: proc_delta_tbl[k]) {
            for (auto kci: pop_contribs_tbl[pd.p])
                apply_contrib_update(kci,pd.delta,update_notify);
            pop_count[pd.p]+=pd.delta;
        }
    }

    void apply(key_type k) { apply(k,[](key_type) {}); }

    value_type propensity(key_type k) {
        const proc_propensity_entry &kp=proc_propensity_tbl[k];
        value_type r=kp.rate;
        for (auto c: kp.counts) r*=c;
        return r;
    }

    friend std::ostream& operator<<(std::ostream &O,const ssa_pp_procsys &sys) {
        // primarily for debugging
        O << "ssa_pp_procsys: n_pop=" << sys.n_pop << ", n_proc=" << sys.n_proc << "\n";
        O << "pop_contribs_tbl:\n";
        size_t idx=0;
        for (const auto &e: sys.pop_contribs_tbl) {
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
        O << "pop_count:\n";
        idx=0;
        for (const auto &c: sys.pop_count) {
            O << "    " << std::setw(6) << std::right << idx++ << ":"
              << ' ' << c << '\n';
        }
        O << "proc_propensity_tbl:\n";
        idx=0;
        for (const auto &e: sys.proc_propensity_tbl) {
            O << "    " << std::setw(6) << std::right << idx++ << ":";
            O << " rate=" << std::left << std::setw(10) << e.rate;
            O << " counts:";
            for (const auto &c: e.counts) 
                O << ' ' << c;
            O << "\n";
        }
        return O;
    }

private:
    size_t n_pop;
    size_t n_proc;

    template <typename Proc>
    void add_proc(const Proc &proc) {
        if (n_proc>=std::numeric_limits<key_type>::max())
            throw rdmini::invalid_value("process index out of bounds");

        key_type key=n_proc++;

        proc_propensity_tbl.resize(n_proc);
        proc_propensity_tbl[key].rate=proc.rate();

        std::vector<pop_index> left_sorted;
        std::map<pop_index,count_type> delta_map;
        for (auto p: proc.left()) {
            if (p>=n_pop) throw rdmini::invalid_value("population index out of bounds");
            --delta_map[p];
            left_sorted.push_back(p);
        }
        for (auto p: proc.right()) {
            if (p>=n_pop) throw rdmini::invalid_value("population index out of bounds");
            ++delta_map[p];
        }

        proc_delta_tbl.resize(n_proc);
        auto &pd_entry=proc_delta_tbl[key];
        for (auto pd: delta_map)
           if (pd.second) pd_entry.push_back({pd.first,pd.second});

        std::sort(left_sorted.begin(),left_sorted.end());
        uint32_t index=0;
        for (auto p: left_sorted) {
            pop_contribs_tbl[p].emplace_back(key,index++);
        }
    }

    // population to process rate factor map: constant per model (i.e. shareable)
    struct proc_contrib_index {
        key_type k; // which process
        uint32_t i; // which slot in rate contribs

        proc_contrib_index(key_type k_,uint32_t i_): k(k_),i(i_) {}
    };
    typedef std::vector<proc_contrib_index> pop_contribs_entry;
    std::vector<pop_contribs_entry> pop_contribs_tbl;

    template <typename F>
    void apply_contrib_update(proc_contrib_index kci,count_type delta,F notify) {
        proc_propensity_tbl[kci.k].counts[kci.i]+=delta;
        notify(kci.k);
    }

    // process to population delta map: constant per model (i.e. shareable)
    struct proc_delta {
        pop_index p; // which population
        int32_t delta; // delta to apply e.g. +1, -1, 0 specific to the reaction and speco
    };
    typedef std::vector<proc_delta> proc_delta_entry;
    std::vector<proc_delta_entry> proc_delta_tbl;

    // population counts (mutable, instance-specific)
    std::vector<count_type> pop_count;
    
    // process rate factor table (mutable, instance-specific)
    struct proc_propensity_entry {
        value_type rate;
        std::array<count_type,max_process_order> counts;

        proc_propensity_entry(): rate(0) { counts.fill(1); }
    };
    std::vector<proc_propensity_entry> proc_propensity_tbl;
};


#endif // ndef  SSA_PP_PROCSYS_H_
