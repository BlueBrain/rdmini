#ifndef SSA_PP_PROCSYS_H_
#define SSA_PP_PROCSYS_H_

#include <cstdint>
#include <vector>
#include <map>
#include <stdexcept>
#include <limits>

#include "rdmini/iterspan.h"
#include "rdmini/rdmodel.h"
#include "rdmini/ssa_common.h"

/** SSA process system that maintains process dependencies
 * factored through populations, and computes propensities
 * on demand from cached factors. */

template <unsigned MaxOrder=3>
struct ssa_pp_procsys {
    typedef uint32_t key_type;
    typedef double value_type;
    typedef uint32_t count_type;
    typedef const std::vector<key_type> &keyset_type;

    static constexpr size_t max_process_order=MaxOrder;
    static constexpr size_t max_population_index=std::numeric_limits<pop_index>::max()-1;
    static constexpr size_t max_count=std::numeric_limits<count_type>::max();
    static constexpr size_t max_participants=max_populations;

    ssa_pp_procsys(): n_proc(0) {}

    template <typename In>
    void define_processes(In b,In e) {
        clear();
        while (b!=e) add_proc(*b++0;
        zero_populations();
    }

    /** Remove all processes, population counts */
    void clear() {
        proc_propensity_tbl.clear();
        for (auto &entry: pop_contribs_tbl) entry.clear();
    }

    size_t size() const { return n_kproc; }
    
    /** Zero all population counts, reset propensity contributions */
    void zero_populations() {
        // init contribs and pop counts
        pop_index n_pop=pop_count.size();
        for (pop_index p=0; p<n_pop; ++p) {
            pop_count[p]=0;
            auto &pc_entry=pop_contribs_tbl[p];
            if (pc_entry.empty()) continue;
 
            // note that indices in pc are grouped by kproc id.
            auto pc_i=pc.begin();
            auto pc_end=pc.end();

            count_type count=0;
            key_type k=pc_iter->k;
            proc_propensity_tbl[k].counts[pc_i->i]=count;
            
            while (++pc_i!=pc_end) {
                if (pc_i->k==k) --count;
                else count=0;

                k=pc_iter->k;
                proc_propensity_tbl[k].counts[pc_iter->i]=count;
            }
        }
    }

    count_type count(size_t p) { return pop_count[p]; }

    template <typename F>
    void set_count(size_t p,count_type c,F update_notify) {
        for (auto kci: pop_contribs_tbl[p])
            apply_contrib_update(kci,pop_count[p]-c,update_notify);
    }

    void set_count(size_t p,count_type c) { set_count(p,c,[](key_type) {}); }

    template <typename F>
    void apply(key_type k,F update_notify) {
        for (auto pd: proc_delta_tbl[k]) 
            for (auto kci: pop_contribs_tbl[pd.p_id])
                apply_contrib_udpate(kci,pd.delta,update_notify);
    }

    void apply(key_type k) { apply(k,[](key_type) {}); }

    value_type propensity(key_type k) {
        const proc_propensity_entry &kp=proc_propensity_tbl[k];
        value_type r=kp.rate;
        for (auto c: kp.counts) r*=c;
        return r;
    }

private:
    typedef uint32_t pop_index;

    template <typename Proc>
    void add_proc(const Proc &proc) {
        if (n_proc>=std::numeric_limits<key_type>::max())
            throw ssa_error("process index out of bounds");

        key_type key=n_proc++;

        proc_propensity_tbl.resize(key);
        proc_propensity_tbl[key].rate=proc.rate();

        std::vector<pop_index> left_sorted;
        std::map<pop_index,ssa_count_type> delta_map;
        for (auto p: proc.left()) {
            if (p>max_population_index) throw ssa_error("population index out of bounds");
            --delta_map[p];
            left_sorted.push_back(p);
        }
        for (auto p: proc.right()) {
            if (p>max_population_index) throw ssa_error("population index out of bounds");
            ++delta_map[p];
        }

        proc_delta_tbl.resize(k_id);
        auto &pd_entry=proc_delta_tbl[k_id];
        for (auto pd: delta_map) pd_entry.push_back({pd.first,pd.second});

        std::sort(left_sorted.begin(),left_sorted.end());
        uint32_t index=0;
        for (auto p: left_sorted) {
            if (p>=pop_contribs_tbl.size()) pop_contribs_tbl.resize(p+1);
            pop_contribs_tbl[p].emplace_back({key,index++});
        }
    }

    template <typename F>
    void apply_contrib_update(proc_contrib_index kci,count_type delta,F notify) {
        proc_contrib_tbl[kci.k_id].counts[kci.i]+=delta;
        notify(kci.k_id);
    }

    size_t n_proc;

    std::vector<count_type> pop_count;
    
    struct proc_contrib_index {
        key_type k; // which process
        uint32_t i; // which slot in rate contribs
    };

    typedef std::vector<proc_contrib_index> pop_contribs_entry;
    std::vector<pop_contribs_entry> pop_contribs_tbl;

    struct proc_propensity_entry {
        value_type rate;
        std::array<count_type,max_process_order> counts;

        proc_propensity_entry(): rate(0) { counts.fill(1); }
    };
    std::vector<prpc_propensity_entry> proc_propensity_tbl;
};


#endif // ndef  SSA_PP_PROCSYS_H_
