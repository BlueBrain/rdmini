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
private:
    typedef ssa_pp_procsys<3> proc_system;
    typedef proc_system::key_type proc_index_type;

    typedef ssa_direct<proc_index_type,double> ssa_selector;
    typedef ssa_selector::event_type event_type;

    struct ksel_update {
        ksel_update(proc_system &sys_,ssa_selector &sel_): sys(sys_), sel(sel_) {}
        proc_system &sys;
        ssa_selector &sel;
        void operator()(proc_index_type k) { sel.update(k, sys.propensity(k)); }
    };

public:
    typedef proc_system::count_type count_type;
    static constexpr unsigned int max_process_order=proc_system::max_process_order;
    static constexpr unsigned int max_participants=proc_system::max_participants;
    static constexpr unsigned int dynamic_range=32; // TODO: replace with correct value ...

    serial_ssa(): t(0), stale(true) {}

    explicit serial_ssa(const rd_model &M, double t0=0) {
        initialise(M,t0);
    }

    struct kproc_info {
        std::vector<size_t> left_,right_; // lhs and rhs of the chemical recaion A + B = AB
        double rate_;

        const std::vector<size_t> &left() const { return left_; }
        const std::vector<size_t> &right() const { return right_; }
        double rate() const { return rate_; }
    };

    void initialise(const rd_model &M, double t0) {
        t=t0;
        stale=true;

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
                ki.rate_=reac.rate/vol;
                for (size_t s_id: reac.left) ki.left_.push_back(species_to_pop(s_id,c_id));
                for (size_t s_id: reac.right) ki.right_.push_back(species_to_pop(s_id,c_id));
            
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
                    ki.left_[0]=species_to_pop(s_id,c_id);
                    ki.right_[0]=species_to_pop(s_id,neighbour.cell_id);
                    ki.rate_=neighbour.diff_coef*M.species[s_id].diffusivity;

                    kp_set.push_back(ki);
                }
            }
        }

        ksys.reset(n_pop);
        ksys.define_processes(kp_set.begin(),kp_set.end());
        ksel.reset(ksys.size());

        // initalise population counts
        // (later, iterate over list of named cell lists for this)
        for (size_t s_id=0; s_id<n_species; ++s_id) {
            double conc=M.species[s_id].concentration;
            for (size_t c_id=0; c_id<n_cell; ++c_id)
                set_count(s_id,c_id, conc*M.cells[c_id].volume);
        }

        // initialise selector with propensities
        for (proc_index_type k=0; k<ksys.size(); ++k)
            ksel.update(k, ksys.propensity(k));
    }

    void set_count(size_t species_id,size_t cell_id,count_type count) {
        ksel_update U(ksys,ksel);
        ksys.set_count(species_to_pop(species_id,cell_id),count,U);
        stale=true;
    }

    count_type count(size_t species_id,size_t cell_id) const {
        return ksys.count(species_to_pop(species_id,cell_id)); 
    }

    // for time interval
    template <typename G>
    double advance(double t_end,G &g) {
        for (;;) {
            get_next(g);
            if (t+next.dt>t_end) break;

            advance(g);
        }

        next.dt-=t_end-t;
        t=t_end;
        return t;
    }

    // for one exact step
    template <typename G>
    double advance(G &g) {
        get_next(g);
        ksys.apply(next.k_id,
            [&](proc_index_type k) { ksel.update(k,ksys.propensity(k)); });

        t+=next.dt;
        stale=true;

        return t;
    }

    friend std::ostream &operator<<(std::ostream &O,const serial_ssa &S) {
        O << S.ksys;
        // O << S.ksel;
        return O;
    }
        
private:
    template <typename G> 
    void get_next(G &g) {
        if (stale) {
            auto ev=ksel.next(g);
            next={ev.key(), ev.dt()};
        }
        stale=false;
    }

    size_t n_species;
    size_t n_reac;
    size_t n_cell;
    size_t n_pop;

    double t; // sim time...
    struct {
        proc_index_type k_id;
        double dt;
    } next;

    bool stale; // true => need to re-poll ksel

    proc_system ksys;
    ssa_selector ksel;

    size_t species_to_pop(size_t species_id,size_t cell_id) const {
        return cell_id*n_species+species_id;
    }
};


#endif // ndef  SERIAL_SSA_H_
