#ifndef SERIAL_SSA_H_
#define SERIAL_SSA_H_

#include <cstdint>
#include <vector>
#include <map>
#include <stdexcept>
#include <limits>

#include "rdmini/iterspan.h"
#include "rdmini/rdmodel.h"
#include "rdmini/ssa_common.h"
#include "rdmini/ssa_pp_procsys.h"
#include "rdmini/ssa_direct.h"

// TODO: split into .h and .cc if we've got fixed template parameters.

struct serial_ssa {
    serial_ssa() {}

    struct kproc_info {
        std::vector<size_t> left_,right_;
        double rate_;

        const std::vector<size_t> &left() { return left_; }
        const std::vector<size_t> &right() { return right_; }
        double rate() { return rate_; }
    };

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

        ksys.define_processes(kp_set.begin(),kp_set.end());
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

    typedef proc_system::count_type count_type;
    static constexpr unsigned int max_order=proc_system::max_order;
    static constexpr unsigned int max_participants=proc_system::max_participants;
    static constexpr unsigned int dynamic_range=32; // TODO: replace with correct value ...

    void set_count(size_t species_id,size_t cell_id,count_type count) {
        ksys.set_pop_count(species_to_pop(species_id,cell_id),count); 
    }

    count_type count(size_t species_id,size_t cell_id)
        return ksys.pop_count(species_to_pop(species_id,cell_id)); 
    }

    template <typename G>
    double advance(double t_end,G &g) {
        while (t_end<=t) advance(g);
        return t;
    }

    template <typename G>
    double advance(G &g) {
        auto next=ksel.next(g);
        ksys.apply(next.k_id,
            [&ksel](kproc_index k) { ksel.update(k,ksys.propensity(k)); });

        t+=next.dt;
        return t;
    }

        
private:
    size_t n_species;
    size_t n_reac;
    size_t n_cell;
    size_t n_pop;

    double t; // sim time...

    typedef ssa_direct ssa_selector;
    typedef ssa_pp_procsys<3> proc_system;

    typedef proc_system::key_type proc_index_type;

    proc_system ksys;
    ssa_selector ksel;

    pop_index species_to_pop(size_t species_id,size_t cell_id) {
        return cell_id*n_species+species_id;
    }

    size_t pop_to_species(size_t species_id,size_t cell_id) {
        return cell_id*n_species+species_id;
    }
};


#endif // ndef  SERIAL_SSA_H_
