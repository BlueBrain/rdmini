#ifndef SSA_PP_PROCSYS_PAR_H_
#define SSA_PP_PROCSYS_PAR_H_

#include <cstdint>
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
#include "rdmini/ssa_common.h"

/** SSA process system that maintains process dependencies
 * factored through populations, and computes propensities
 * on demand from cached factors.
 *
 * Type parameters:
 *   key_type: unsigned integral type for enumerating processses. uint32_t probably good choice.
 *   value_type: floating point type representing propensitites. double or float.
 *   count_type: signed integral type for representing population counts and deltas. uint32_t probably suffices.
 */
 
/** Per-instance process propensity calculation table.
 *
 * For each process index (key_type), have a rate and MaxOrder factors.
 * Population changes update the factors.
 * Propensity is calculated as the product of the rate and the factors.
 */

template <typename value_type, typename count_type, unsigned MaxOrder>
struct proc_propensity_entry {
    value_type rate;
    std::array<count_type,MaxOrder> counts;

    proc_propensity_entry(): rate(0) { counts.fill(1); }
};

template <typename value_type, typename count_type, unsigned MaxOrder>
using proc_propensity_table_view=proc_propensity_entry<value_type, count_type, MaxOrder> *;

/** Per-instance population counts */

template <typename count_type>
using pop_count_view=count_type *;

/** Per-model map from population indices to process factor contributions */

template <typename key_type>
struct proc_contrib_index {
    key_type k; // which process
    uint32_t i; // which slot in rate contribs

    proc_contrib_index(key_type k_,uint32_t i_): k(k_),i(i_) {}
};

template <typename key_type>
using pop_contribs_entry=std::vector<proc_contrib_index<key_type>>;

template <typename key_type>
using pop_contribs_table=std::vector<pop_contribs_entry<key_type>>;

template <typename key_type>
using pop_contribs_table_view=const pop_contribs_entry<key_type> *;

/** Per-model map from process indices to population deltas */

template <typename pop_index>
struct proc_delta {
    pop_index p; // which population
    int32_t delta; // delta to apply
};

template <typename pop_index>
using proc_delta_entry=std::vector<proc_delta<pop_index>>;

template <typename pop_index>
using proc_delta_table=std::vector<proc_delta_entry<pop_index>>;

template <typename pop_index>
using proc_delta_table_view=const proc_delta_entry<pop_index> *;

/** One instance of a process system for a given model.
 *
 * Population counts and propensity tables are per-instance;
 * Population contribution table and process to population delta tables are common
 * accross instances.
 */

template <typename KeyType, typename ValueType, typename CountType, unsigned MaxOrder>
struct ssa_pp_procsys_instance {
    typedef KeyType key_type;
    typedef ValueType value_type;
    typedef CountType count_type;

private:
    typedef uint32_t pop_index;
    typedef proc_propensity_table_view<value_type, count_type, MaxOrder> pp_table;
    typedef pop_count_view<count_type> pop_count_data;

public:
    ssa_pp_procsys_instance(): n_pop(0), n_proc(0) {}

    ssa_pp_procsys_instance(size_t n_pop_,
                            size_t n_proc_,
                            pop_contribs_table_view<key_type> pop_contribs_tbl_,
                            proc_delta_table_view<pop_index>  proc_delta_tbl_,
                            pp_table                          prop_propensity_tbl_,
                            pop_count_data                    pop_count_):
        n_pop(n_pop_),
        n_proc(n_proc_),
        pop_contribs_tbl(pop_contribs_tbl_),       // read-only view
        proc_delta_tbl(proc_delta_tbl_),           // read-only view
        proc_propensity_tbl(prop_propensity_tbl_), // mutable view
        pop_count(pop_count_)                      // mutable view
    {
    }

    /** Zero all population counts, reset propensity contributions */
    void zero_populations() {
        for (pop_index p=0; p<n_pop; ++p) {
            pop_count[p]=0;
            const auto &pc_entry=pop_contribs_tbl[p];
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

    size_t size() const { return n_proc; }
    count_type count(size_t p) const { return pop_count[p]; }

    template <typename F>
    void set_count(size_t p, count_type c, F update_notify) {
        for (auto kci: pop_contribs_tbl[p])
            apply_contrib_update(kci, c-pop_count[p], update_notify);
        pop_count[p]=c;
    }

    void set_count(size_t p,count_type c) { set_count(p,c,[](key_type) {}); }

    // batched interface to population counts
    template <typename Out>
    void copy_counts(Out out) const {
        std::copy(&pop_count[0], &pop_count[0]+n_pop, out);
    }

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
        const auto &kp=proc_propensity_tbl[k];
        value_type r=kp.rate;
        for (auto c: kp.counts) r*=c;
        return r;
    }

    friend std::ostream& operator<<(std::ostream &O,const ssa_pp_procsys_instance &sys) {
        // primarily for debugging
        O << "ssa_pp_procsys_instance:\n";
        O << "pop_count:\n";
        for (size_t idx=0; idx<sys.n_pop; ++idx) {
            O << "    " << std::setw(6) << std::right << idx << ":"
              << ' ' << sys.pop_count[idx] << '\n';
        }
        O << "proc_propensity_tbl:\n";
        for (size_t idx=0; idx<sys.n_proc; ++idx) {
            const auto &e=sys.proc_propensity_tbl[idx];
            O << "    " << std::setw(6) << std::right << idx << ":";
            O << " rate=" << std::left << std::setw(10) << e.rate;
            O << " counts:";
            for (const auto &c: e.counts) 
                O << ' ' << c;
            O << "\n";
        }
        return O;
    }

private:
    template <typename F>
    void apply_contrib_update(proc_contrib_index<key_type> kci,count_type delta,F notify) {
        proc_propensity_tbl[kci.k].counts[kci.i]+=delta;
        notify(kci.k);
    }

    size_t n_proc;
    size_t n_pop;

    // immutable:
    proc_delta_table_view<pop_index> proc_delta_tbl;
    pop_contribs_table_view<key_type> pop_contribs_tbl;

    // mutable:
    pp_table proc_propensity_tbl;
    pop_count_data pop_count;
};


/** Represents set of ssa_pp_procsys_instance objects, sharing the same model. */

template <unsigned MaxOrder=3>
struct ssa_pp_procsys_par {
    typedef uint32_t key_type;
    typedef double value_type;
    typedef int32_t count_type;

    typedef ssa_pp_procsys_instance<key_type, value_type, count_type, MaxOrder> instance_type;

private:
    typedef uint32_t pop_index;
    typedef proc_propensity_entry<value_type, count_type, MaxOrder> pp_entry;
    typedef proc_propensity_table_view<value_type, count_type, MaxOrder> pp_table_view;

public:
    static constexpr size_t max_process_order=MaxOrder;
    static constexpr size_t max_population_index=std::numeric_limits<pop_index>::max()-1;
    static constexpr size_t max_count=std::numeric_limits<count_type>::max();
    static constexpr size_t max_participants=max_population_index;

    ssa_pp_procsys_par(): n_instance(0), n_pop(0) {}

    template <typename In>
    ssa_pp_procsys_par(size_t n_instance_, size_t n_pop_, In b, In e): n_instance(n_instance_), n_pop(n_pop_), n_proc(0) {
        if (n_pop_>0 && n_pop_-1>max_population_index)
            throw ssa_error("population index out of bounds");

        pop_count_instances.assign(n_instance,std::vector<count_type>(n_pop,0));
        pop_contribs_tbl.resize(n_pop);

        std::vector<pp_entry> pp_tbl_template;
        while (b!=e) add_proc(*b++, pp_tbl_template);
        
        proc_propensity_tbl_instances.assign(n_instance,pp_tbl_template);

        // do any packing of built process-population dependency structures here ...

        for (size_t i=0; i<n_instance;++i) instance(i).zero_populations();
    }


    size_t size() const { return n_proc; }
    size_t instances() const { return n_instance; }

    instance_type instance(size_t i) {
        pop_contribs_table_view<key_type> pop_contribs_view=&pop_contribs_tbl[0];
        proc_delta_table_view<pop_index>  proc_delta_view=&proc_delta_tbl[0];
        pp_table_view                     pp_tbl={&proc_propensity_tbl_instances[i][0]};
        pop_count_view<count_type>        pop_counts={&pop_count_instances[i][0]};

        return instance_type(n_pop, n_proc, pop_contribs_view, proc_delta_view, pp_tbl, pop_counts);
    }
    
    friend std::ostream& operator<<(std::ostream &O,const ssa_pp_procsys_par &sys) {
        // primarily for debugging
        O << "ssa_pp_procsys_par: n_instance=" << sys.n_instance << ", n_pop=" << sys.n_pop << ", n_proc=" << sys.n_proc << "\n";
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
        // const cheating for debugging -- TODO: make a better interface!
        O << "instance 0:\n" << const_cast<ssa_pp_procsys_par &>(sys).instance(0) << "\n";
        return O;
    }

private:
    size_t n_instance;
    size_t n_pop;
    size_t n_proc;

    /** Called in setup phase */
    template <typename Proc>
    void add_proc(const Proc &proc, std::vector<pp_entry> &proc_propensity_tbl) {
        if (n_proc>=std::numeric_limits<key_type>::max())
            throw ssa_error("process index out of bounds");

        key_type key=n_proc++;

        proc_propensity_tbl.resize(n_proc);
        proc_propensity_tbl[key].rate=proc.rate();

        std::vector<pop_index> left_sorted;
        std::map<pop_index,count_type> delta_map;
        for (auto p: proc.left()) {
            if (p>=n_pop) throw ssa_error("population index out of bounds");
            --delta_map[p];
            left_sorted.push_back(p);
        }
        for (auto p: proc.right()) {
            if (p>=n_pop) throw ssa_error("population index out of bounds");
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

    // model-specific data, shared across instances
    proc_delta_table<pop_index> proc_delta_tbl;
    pop_contribs_table<key_type> pop_contribs_tbl;

    // instance-specifc data
    std::vector<std::vector<count_type>> pop_count_instances;
    std::vector<std::vector<pp_entry>> proc_propensity_tbl_instances;
};


#endif // ndef  SSA_PP_PROCSYS_PAR_H_
