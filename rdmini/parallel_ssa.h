#ifndef SERIAL_SSA_H_
#define SERIAL_SSA_H_

#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "rdmini/iterspan.h"
#include "rdmini/rdmodel.h"
#include "rdmini/exceptions.h"
#include "rdmini/ssa_direct.h"

#include "rdmini/ssa_pp_procsys.h"

// TODO: split into .h and .cc if we've got fixed template parameters.

struct parallel_ssa {
private:
    typedef ssa_pp_procsys<3> proc_system;
    typedef proc_system::key_type proc_index_type;

    typedef ssa_direct<proc_index_type,double> ssa_selector;
    typedef ssa_selector::event_type event_type;

    struct ksel_updater_f {
        ksel_updater_f(proc_system &sys_,ssa_selector &sel_,size_t instance_):
            sys(sys_), sel(sel_),instance(instance_) {}
        proc_system &sys;
        ssa_selector &sel;
        size_t instance;

        void operator()(proc_index_type k) { sel.update(k, sys.propensity(k)); }
    };

    ksel_updater_f ksel_update(size_t instance) {
        return ksel_updater_f(ksys,states[instance].ksel,instance);
    }


public:
    typedef proc_system::count_type count_type;
    static constexpr unsigned int max_process_order=proc_system::max_process_order;
    static constexpr unsigned int max_participants=proc_system::max_participants;
    static constexpr unsigned int dynamic_range=32; // TODO: replace with correct value ...

    parallel_ssa() {}

    explicit parallel_ssa(size_t n_instances,const rd_model &M, double t0=0) {
        initialise(n_instances,M,t0);
    }

    struct kproc_info {
        std::vector<size_t> left_,right_;
        double rate_;

        const std::vector<size_t> &left() const { return left_; }
        const std::vector<size_t> &right() const { return right_; }
        double rate() const { return rate_; }
    };

    void initialise(size_t n_instances_,const rd_model &M, double t0) {
        n_instances=n_instances_;

        n_species=M.n_species();
        n_reac=M.n_reactions();
        n_cell=M.n_cells();

        n_pop=n_species*n_cell;
        
        std::vector<kproc_info> kp_set;

        // cell-local populations and reactions
        for (size_t c_id=0; c_id<n_cell; ++c_id) {
            double vol=M.cells[c_id].volume;
            for (const auto &reac: M.reactions) {
                kproc_info ki;
                int order=(int)reac.left.size();
                ki.rate_=reac.rate*std::pow(vol,1-order);
                for (size_t s_id: reac.left) ki.left_.push_back(species_to_pop_id(s_id,c_id));
                for (size_t s_id: reac.right) ki.right_.push_back(species_to_pop_id(s_id,c_id));
            
                kp_set.push_back(ki);
            }
        }

        // diffusion processes
        for (size_t c_id=0; c_id<n_cell; ++c_id) {
            for (auto neighbour: M.cells[c_id].neighbours) {
                if (neighbour.diff_coef==0) continue;

                kproc_info ki;
                ki.left_.resize(1);
                ki.right_.resize(1);
                for (size_t s_id=0; s_id<n_species; ++s_id) {
                    ki.left_[0]=species_to_pop_id(s_id,c_id);
                    ki.right_[0]=species_to_pop_id(s_id,neighbour.cell_id);
                    ki.rate_=neighbour.diff_coef*M.species[s_id].diffusivity;

                    kp_set.push_back(ki);
                }
            }
        }

        ksys=proc_system(n_instances);
        ksys.add(kp_set.begin(),kp_set.end());

        states.resize(n_instances);
        #pragma omp parallel for
        for (size_t i=0; i<n_instances; ++i) {
            auto &state=states[i];

            state.t=t0;
            state.stale=true;

            state.ksel.reset(ksys.size());

            // initalise population counts
            // (later, iterate over list of named cell lists for this)
            for (size_t s_id=0; s_id<n_species; ++s_id) {
                double conc=M.species[s_id].concentration;
                for (size_t c_id=0; c_id<n_cell; ++c_id)
                    ksys.set_count(species_to_pop_id(s_id,c_id),conc*M.cells[c_id].volume,i);
            }

            // initialise selector with propensities
            auto update=ksel_update(i);
            for (proc_index_type k=0; k<ksys.size(); ++k) update(k);
        }
    }

    void set_count(size_t instance,size_t species_id,size_t cell_id,count_type count) {
        auto &state=states[instance];

        ksys.set_count(species_to_pop_id(species_id,cell_id),count,ksel_update(instance),instance);
        state.stale=true;
    }

    count_type count(size_t instance,size_t species_id,size_t cell_id) const {
        return ksys.count(species_to_pop_id(species_id,cell_id)); 
    }

    std::result_of<decltype(&proc_system::counts)(proc_system,size_t)>::type counts(size_t instance) const {
        return ksys.counts(instance);
    }

    template <typename G>
    double advance(size_t instance,double t_end,G &g) {
        auto &state=states[instance];
        auto update=ksel_update(instance);

        for (;;) {
            state.get_next(g);
            if (state.t+state.next_dt>t_end) break;

            ksys.apply(state.next_k_id,update,instance);
            state.t+=state.next_dt;
            state.stale=true;
        }

        state.next_dt-=t_end-state.t;
        state.t=t_end;
        return state.t;
    }

    template <typename G>
    double advance(size_t instance,G &g) {
        auto &state=states[instance];

        state.get_next(g);
        ksys.apply(state.next_k_id,ksel_update(instance),instance);
        state.t+=state.next_dt;
        state.stale=true;

        return state.t;
    }

    size_t population_size() const { return n_pop; }
    size_t instances() const { return n_instances; }

    std::pair<size_t,size_t> pop_to_species_id(size_t pop_id) const {
        return std::pair<size_t,size_t>(pop_id%n_species,pop_id/n_species);
    }

    size_t species_to_pop_id(size_t species_id,size_t cell_id=0) const {
        return cell_id*n_species+species_id;
    }

    friend std::ostream &operator<<(std::ostream &O,const parallel_ssa &S) {
        O << S.ksys;
        return O;
    }
        
private:
    size_t n_instances;
    size_t n_species;
    size_t n_reac;
    size_t n_cell;
    size_t n_pop;


    struct instance_state {
        double t;
        ssa_selector ksel;

        bool stale;
        proc_index_type next_k_id;
        double next_dt;

        template <typename G> 
        void get_next(G &g) {
            if (stale) {
                auto ev=ksel.next(g);
                next_k_id=ev.key();
                next_dt=ev.dt();
            }

            stale=false;
        }
    };

    proc_system ksys;
    std::vector<instance_state> states;
};


#endif // ndef  SERIAL_SSA_H_
