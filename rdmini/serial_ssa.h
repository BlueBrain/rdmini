#ifndef SERIAL_SSA_H_
#define SERIAL_SSA_H_

#include <cstdint>
#include <vector>
#include <map>
#include <stdexcept>

#include "iterspan.h"
#include "rdmodel.h"

static constexpr unsigned ssa_max_order=4;
typedef double ssa_propensity_type;
typedef ptrdiff_t ssa_count_type;
typedef uint32_t kproc_index;
typedef uint32_t pop_index;

struct ssa_error: std::runtime_error {
    ssa_error(const std::string &what_str): std::runtime_error(what_str) {}
};

/** Reaction or diffusion process info required for building
 * SSA data structures. */

struct kproc_info {
    std::vector<pop_index> left,right;
    ssa_propensity_type rate;
};

/** Maintains state for computing/updating kinetic process propensities */

struct kproc_system {
    template <typename In>
    kproc_system(pop_index n_pop_,In kproc_begin,In kproc_end):
        n_pop(n_pop_), n_kproc(0), pop_contribs_tbl(n_pop_), pop_count(n_pop_,0)
    {
        for (const auto &ki: make_span(kproc_begin,kproc_end))
            add_kproc(ki);

        zero_pop_counts();
    }

    /** Zero all population counts, reset propensity contributions */
    void zero_pop_counts() {
        // init contribs and pop counts
        for (pop_index p=0; p<n_pop; ++p) {
            pop_count[p]=0;
            auto &pc=pop_contribs_tbl[p];

            // note that indices in pc are grouped by kproc id.
            auto pc_iter=pc.begin(),pc_end=pc.end();
            if (pc_iter==pc_end) continue;

            ssa_count_type count=0;
            kproc_index k=pc_iter->k_id;
            kproc_propensity_tbl[k].counts[pc_iter->i]=count;
            
            while (++pc_iter!=pc_end) {
                if (pc_iter->k_id==k) --count;
                else count=0;

                k=pc_iter->k_id;
                kproc_propensity_tbl[k].counts[pc_iter->i]=count;
            }
        }
    }

    /** Remove all kprocs */
    void clear() {
        kproc_propensity_tbl.clear();
        for (auto &entry: pop_contribs_tbl) entry.clear();
    }

    /** Retrieve population count */
    ssa_count_type pop_count(pop_index p) const { return pop_count[p]; }

    /** Adjust population count
     *
     * K is a functional object that takes propensity updates,
     * K(kproc_id,propensity)
     */
    template <typename KSelectorUpdate>
    void set_pop_count(pop_index p,ssa_count_type c,KSelectorUpdate S) {
        auto delta=c-pop_count(p);
        for (auto kci: pop_contribs_tbl[p])
            apply_contrib_udpate(kci,delta,S);
    }

    /** Run one kproc
     *
     * K is a functional object that takes propensity updates,
     * K(kproc_id,propensity)
     */
    template <typename KSelectorUpdate>
    void run_kproc(kproc_index k,KSelectorUpdate S) {
        for (auto pd: kproc_pop_delta_tbl[k]) 
            for (auto kci: pop_contribs_tbl[pd.p_id])
                apply_contrib_udpate(kci,pd.delta,S);
    }

    /** Number of kprocs */
    size_t size() const { return n_kproc; }

private:
    size_t n_pop;
    size_t n_kproc;

    struct kproc_contrib_index {
        kproc_index k_id; // which kproc
        uint32_t i;       // which slot in rate contribs

        friend bool operator<(const kproc_contrib_index &a,const kproc_contrib_index &b) {
            return a.k_id<b.k_id || a.k_id==b.k_id && a.i<b.i;
        }
    };

    typedef std::vector<kproc_contrib_index> pop_contribs_set;
    std::vector<pop_contribs_set> pop_contribs_tbl;

    struct kproc_propensity_entry {
        ssa_propensity_type rate;
        std::array<ssa_count_type,ssa_max_order> counts;
        
        ssa_propensity_type propensity() const {
            ssa_propensity_type p=rate;
            for (auto c: counts) p*=static_cast<ssa_propensity_type>(c);
            return p;
        }
    };

    std::vector<kproc_propensity_entry> kproc_propensity_tbl;

    struct pop_delta {
        pop_index p_id;
        ssa_count_type delta;
    };

    typedef std::vector<pop_delta> kproc_delta_entry;
    std::vector<kproc_delta_entry> kproc_delta_tbl;

    /* Extend pop_to_contrib_index and kproc_contrib_tbl */
    void add_kproc(const kproc_info &ki) {
        // TODO: add size check w.r.t. max kproc_index

        kproc_index k_id=n_kproc++;

        kproc_contrib_tbl.resize(k_id);
        kproc_contrib_tbl[k_id].rate=ki.rate;

        std::vector<pop_index> left(ki.left);
        std::sort(left.begin(),left.end());

        uint32_t index=0;
        for (auto p: left) {
            if (p>=n_pop) throw ssa_error("pop_index out of bounds");
            pop_contribs_tbl[p].emplace_back({k_id,index++});
        }

        kproc_delta_tbl.resize(k_id);
        auto &kd_entry=kproc_delta_tbl[k_id];

        std::map<pop_index,ssa_count_type> delta_map;
        for (auto p: ki.left) --delta_map[p];
        for (auto p: ki.right) ++delta_map[p];

        for (auto pd: delta_map) kd_entry.push_back({pd.first,pd.second});
    }
        

    void apply_delta(kproc_contrib_index kci,ssa_count_type delta) {
        kproc_contrib[kci.k_id].counts[kci.i]+=delta;
    }
};

struct ssa_direct_sampler {
    // simplest ladder implementation
    kselector(): n_kproc(0) {}

    struct event {
        kproc_index k_id;
        double dt;
    };

    template <typename RNG>
    event next_event(RNG &rng) {
        std::uniform_real_distribution U(0,propensity_sum);
        ssa_propensity_type x=U(rng);

        event ev;
        kproc_idx i;
        for (i=0; i<n_kproc: ++i) {
            x-=propensity[i];
            if (x>=0) continue;
        }
        if (i>=n_knproc) throw ssa_error("fell off propensity ladder (rounding?)");

        ev.k_id=i;
        ev.dt=E(rng)/propensity_sum;
        return ev;
    }

    void reset(size_t nk) {
        n_kproc=nk;
        ladder.resize(n_kproc);
        for (size_t i=0; i<nk; ++i) ladder[i]={i,0.0};
    }

    void update(kproc_id k,ssa_propensity_type r) {
        ssa_propensity_type &p=propensity[k];
        propensity_sum+=r-p;
        p=r;
    }

    std::exponential_distribution<ssa_propensity_type> E;
    std::vector<ssa_propensity_type> propensity;
    ssa_propensity_type propensity_sum=0.0;
};

typedef ssa_direct_sampler kselector;

typedef <typename RNG>
struct serial_ssa {
    explicit serial_ssa(RNG rng_): rng(rng_) {}

    void set_seed(typename RNG::value_type v) { rng.set_seed(v); }
    template <typename Sseq>
    void set_seed(Sseq &sseq) { rng.set_seed(sseq); }
    
    void initialise(const rd_model &M,double t0): t(t0) {
        n_species=M.species.size();
        n_reac=M.reactions.size();
        n_cell=M.geom.size();

        n_pop=n_species*n_cell;
        
        std::vector<kproc_info> kp_set;

        // cell-local populations and reactions
        for (size_t c_id=0; c_id<n_cell; ++c_id) {
            double vol=M.geom[c_id].volume;
            for (const auto &reac: M.reac) {
                kproc_info ki;
                ki.rate=reac.rate/vol;
                for (size_t s_id: reac.left) ki.left.push_back(species_to_pop(s_id,c_id));
                for (size_t s_id: reac.right) ki.right.push_back(species_to_pop(s_id,c_id));
            
                kp_set.push_back(ki);
            }
        }

        // diffusion processes
        for (size_t c_id=0; c_id<n_cell; ++c_id) {
            for (auto neigbour: M.geom[c_id].neighbours()) {
                if (neighbour.diff_coef==0) continue;

                kproc_info ki;
                ki.left.resize(1);
                ki.right.resize(1);
                for (size_t s_id=0; s_id<n_species; ++s_id) {
                    ki.left[0]=species_to_pop(s_id,c_id);
                    ki.right[0]=species_to_pop(s_id,neighbour.cell_id);
                    ki.rate=neighbour.diff_coef*M.species[s_id].diffusivity;

                    kp_set.push_back(ki);
                }
            }
        }

        ksys.reset(new kproc_system(n_pop,kp_set.begin(),kp_set.end()));
        ksel.reset(ksys.size());

        // initalise population counts
        // (later, iterate over list of subvolumes for this)
        for (auto conc_spec: M.initial) {
            size_t s_id=conc_spec.species_id;
            for (size_t c_id=0; c_id<n_cell; ++c_id)
                set_count(s_id,c_id,conc_spec.count(c_id),
                        [&ksel](kproc_index k,ssa_propensity_type p) { ksel.update(k,p); });
        }
    }

    typedef ssa_count_type count_type;
    static constexpr unsigned int max_order=ssa_max_order;
    static constexpr unsigned int max_participants=std::numeric_limits<pop_index>::max();
    static constexpr unsigned int dynamic_range=32; // not sure this is true yet ...

    void set_count(size_t species_id,size_t cell_id,count_type count) {
        ksys.set_pop_count(species_to_pop(species_id,cell_id),count); 
    }

    count_type count(size_t species_id,size_t cell_id)
        return ksys.pop_count(species_to_pop(species_id,cell_id)); 
    }

    double advance(double t_end) {
        while (t_end<=t) advance();
        return t;
    }

    double advance() {
        auto next=ksel.next(rng);
        ksys.run_kproc(next.k_id,
            [&ksel](kproc_index k,ssa_propensity_type p) { ksel.update(k,p); });

        t+=next.dt;
        return t;
    }

        
private:
    RNG rng;

    size_t n_species;
    size_t n_reac;
    size_t n_cell;
    size_t n_pop;

    double t; // sim time...

    std::unique_ptr<kproc_system> ksys;
    ksel ksel;

    pop_index species_to_pop(size_t species_id,size_t cell_id) {
        return cell_id*n_species+species_id;
    }

    size_t pop_to_species(size_t species_id,size_t cell_id) {
        return cell_id*n_species+species_id;
    }
};


#endif // ndef  SERIAL_SSA_H_
